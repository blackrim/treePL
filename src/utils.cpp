
#include <vector>
#include <map>
#include <stack>
#include <iostream>
#include <algorithm>
#include <fstream>

#include <time.h>
#include <math.h>
#include <limits>
#include <stdlib.h>
#include "node.h"
#include "utils.h"
#include "tree.h"
#include "optimize_nlopt.h"
#include "optimize_tnc.h"
#include "pl_calc_parallel.h"
#include "siman_calc_par.h"
#include "tree_reader.h"
#include "optim_options.h"

using namespace std;
//changed the postorder nodes to the tree, make sure that this works
void set_mins_maxs(Tree * tr,map<Node *,double> * mins, map<Node *,double> * maxs){
    double LARGE = 10000;
    double SMALL = 0.0;
    for (unsigned int i=0;i<tr->getNodeCount();i++){
	Node * nd = tr->getNode(i);
	if(nd->isInternal()){
	    //do the min
	    if((*mins).count(nd) == 0){
		double ymin = SMALL;
		for(int j=0;j<nd->getChildCount();j++){
		    if(nd->getChild(j).getChildCount()>0){//child is internal
			if((*mins)[&nd->getChild(j)] > ymin){
			    ymin = (*mins)[&nd->getChild(j)];
			}
		    }
		}
		(*mins)[nd] = ymin;
		nd->min = ymin;nd->minb = true;
	    }
	    //do the max
	    if((*maxs).count(nd) == 0){
		Node * par = nd;
		double ymax = LARGE;
		while(par->hasParent()){
		    par = par->getParent();
		    double tmax = ymax;
		    if((*maxs).count(par) > 0)
			tmax = (*maxs)[par];
		    if(tmax < ymax)
			ymax = tmax;
		}
		(*maxs)[nd] = ymax;
		nd->max = ymax;nd->maxb = true;
	    }
	    //cout << nd->getName() << " "<< nd->min << " " << (*mins)[nd]<< " " << nd->max << " " << (*maxs)[nd]<< endl;
	}
    }
}

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters){
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void TrimSpaces( string& str)  {
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if(( string::npos == startpos ) || ( string::npos == endpos))
    {
	str = "";
    }
    else
	str = str.substr( startpos, endpos-startpos+1 );

}

double get_start_rate(Tree * tree, vector<double> * durations){
    double trelen = 0;
    double tredur = 0;
    for (int i = 0;i< tree->getNodeCount();i++)
	trelen += tree->getNode(i)->getBL();
    for(int i = 0; i< durations->size(); i++)
	tredur += durations->at(i);
    return trelen/tredur;
}

/**
   need to examine this
   min should probably be much smaller
 **/
void process_initial_branch_lengths(Tree * tree, bool collapse, int numsites){
    double min = 1.;
    for(int i=0;i<tree->getExternalNodeCount();i++){
	if(tree->getExternalNode(i)->getBL()<min/numsites){
	    tree->getExternalNode(i)->setBL(min/numsites);
	    cout << "tiny branch length at "<< tree->getExternalNode(i)->getName()<< ". setting to " << (min/numsites) << endl;
	    //exit(0);
	}
    }
    if (collapse == false){
	for(int i=0;i<tree->getNodeCount();i++){
	    if(tree->getNode(i)->getBL()<min/numsites){
		tree->getNode(i)->setBL(min/numsites);
		cout << "tiny branch length at internal node. setting to " << (min/numsites) << endl;
		 //	 exit(0);
	    }
	}
    }else{
	/*
	 * or collapse
	 */
	bool smallbranches = true;
	while(smallbranches){
	    int icount = 0;
	    for(int i=0;i<tree->getNodeCount();i++){
		if(tree->getNode(i)->getBL()<1./numsites && tree->getNode(i)!= tree->getRoot()){
		    Node * pr = tree->getNode(i)->getParent();
		    vector<Node *> children = tree->getNode(i)->getChildren();
		    for(unsigned int j=0;j<children.size();j++){
			Node * ch = children[j];
			pr->addChild(*ch);
		    }
		    pr->removeChild(*tree->getNode(i));
		    tree->processRoot();
		}else{
		    icount += 1;
		}
	    }
	    if(icount == tree->getNodeCount())
		smallbranches = false;
	}
    }
}

void apply_node_label_preorder(Tree * tr){
    stack<Node *> prestack;
    prestack.push(tr->getRoot());
    int count = 0;
    while(!prestack.empty()){
	Node * tnode = prestack.top();
	prestack.pop();
	tnode->setNumber(count);
	count++;
	for(int i=0;i<tnode->getChildCount();i++){
	    prestack.push(&tnode->getChild(i));
	}
    }
}

void calc_char_durations(Tree * tr, int numsites){
#pragma omp parallel for num_threads(8)
    for(int i=0;i<tr->getNodeCount();i++){
//	cout << i << " " << i/float(tr->getNodeCount()) << endl;
	Node * tnode = tr->getNodeByNodeNumber(i);
	tnode->char_duration = round(tnode->getBL()*numsites);
	tnode->log_fact_char_dur = logFact(tnode->char_duration);
    }
}

void setup_date_constraints(Tree * tr, map<Node *, double> * inmins, map<Node *, double> * inmaxs){
    double maxcon = 0;
    for(int i=0;i<tr->getNodeCount();i++){
	Node * nd = tr->getNodeByNodeNumber(i);
	if(nd->getChildCount()==0){
	    nd->date = 0;
	    nd->free = false;
	    nd->minb = false; nd->maxb = false;
	}else if((*inmins).count(nd) > 0 || (*inmaxs).count(nd) > 0){
	    nd->free=true;
	    if ((*inmins).count(nd) > 0){
		nd->min = (*inmins)[nd];
		if (nd->min > maxcon)
		    maxcon = nd->min;
	    }if ((*inmaxs).count(nd) > 0){
		nd->max = (*inmaxs)[nd];
		if (nd->max > maxcon)
		    maxcon = nd->max;
	    }if ((*inmins).count(nd) > 0 && (*inmaxs).count(nd) > 0){
		if ((*inmaxs)[nd] == (*inmins)[nd]){//fixage
		    nd->date = (*inmins)[nd];
		    nd->free = false;
		}
	    }
	}else{
	    nd->free = true;
	}
    }
}

double logFact(double k){
    return k*(log(k)-1)+log(sqrt(6.28318*k));
}

void extract_tree_info(Tree * tr, vector<int> * free, vector<int> * parent_nds_ints, 
		       vector<int> * child_counts, vector<double> * char_durations, 
		       vector<double> * log_fact_char_durations, vector<double> * min, 
		       vector<double> * max, vector< vector<int> > * children_vec,
		       vector<double> * pen_min, vector<double> * pen_max){
    //should return based on the number in the vector
    //preorder
    free->clear();
    parent_nds_ints->clear();
    char_durations->clear();
    log_fact_char_durations->clear();
    min->clear();
    max->clear();
    child_counts->clear();
    children_vec->clear();
    for(int i=0;i<tr->getNodeCount();i++){
	Node * tnode = tr->getNodeByNodeNumber(i);
	free->push_back(tnode->free);
	if(tnode == tr->getRoot())
	    parent_nds_ints->push_back(-1);
	else
	    parent_nds_ints->push_back(tnode->getParent()->getNumber());
	char_durations->push_back(tnode->char_duration);
	log_fact_char_durations->push_back(tnode->log_fact_char_dur);
	child_counts->push_back(tnode->getChildCount());
	vector<int> childs;
	for (int j=0;j<tnode->getChildCount();j++){
	    childs.push_back(tnode->getChild(j).getNumber());
	}
	children_vec->push_back(childs);
	if(tnode->minb == true)
	    min->push_back(tnode->min);
	else
	    min->push_back(0);//min->push_back(-1);
	if(tnode->maxb == true)
	    max->push_back(tnode->max);
	else
	    max->push_back(1e+20);//max->push_back(-1);
	if(tnode->pen_maxb == true)
	    pen_max->push_back(tnode->pen_max);
	else
	    pen_max->push_back(-1);
	if(tnode->pen_minb == true)
	    pen_min->push_back(tnode->pen_min);
	else
	    pen_min->push_back(-1);
    }
    //ASSERT THAT EACH OF THE VECTORS IS LONG ENOUGH
    /*for(int i=0;i<preorder_nds_ints->size();i++){
 	if(i!=0){
 	    cur_param_ord->push_back(rates->at(i));
 	}
    }
    for(int i=0;i<preorder_nds_ints->size();i++){
 	if(free->at(i)==1)
 	    cur_param_ord->push_back(dates->at(i));
	    }*/
    
}

//get_back_freeparams
int generate_param_order_vector(vector<int> * freeparams, bool lf, 
				 vector<int> * cvnodes, vector<int> * free){
    int icount = 0;
    bool cv = false;
    int numnodes = free->size();
    if(cvnodes != NULL){
	if(cvnodes->size() == numnodes)
	    cv = true;
    }
    //rates
    for(int i=0;i<numnodes;i++){
	if (i == 0){
	    freeparams->push_back(-1);
	}else{
	    if(cv == true){
		if((*cvnodes)[i] == 0){
		    freeparams->push_back(icount);
		    icount += 1;
		}else{
		    freeparams->push_back(-1);
		}
	    }else{
		freeparams->push_back(icount);
		if(lf == false){
		    icount += 1;
		}
	    }
	}
    }
    if (lf == true)
	icount += 1;
    //dates
    for(int i=0;i<numnodes;i++){
	if(free->at(i)==1){
 	    freeparams->push_back(icount);
	    icount += 1;
	}else{
	    freeparams->push_back(-1);
	}
    }
/*    cout << icount << " " << freeparams->size() << endl;
    for(int i=0;i<freeparams->size();i++)
	cout << freeparams->at(i) << " ";
    cout << endl;*/
    return icount;
}

void get_feasible_start_dates(Tree * tr, vector<double> * dates, vector<double> * rates,
			      vector<double> * durations){
    /* Set up some feasible times and stores them in tree. Present-day is time 0.*/
    Node * root = tr->getRoot();
    double minTime=0.0,maxTime=0.0;
    minTime = root->min; maxTime = root->max;
    if (root->free == true){
	if(root->max < 10000)
	    root->date=minTime+(root->max-minTime)*(0.02+myRand()*0.96);
	else
	    root->date=minTime*1.25;
    }
    for(int i=0;i<root->getChildCount();i++){
	a_feasible_time(&root->getChild(i),root->date);
    }
    //setup the vectors
    dates->clear();
    rates->clear();
    durations->clear();
    for(int i=0;i<tr->getNodeCount();i++){
	Node * nd = tr->getNodeByNodeNumber(i);
	dates->push_back(nd->date);
	rates->push_back(nd->rate);
	double duration = 0;
	if (nd->hasParent()){
	    duration = nd->getParent()->date - nd->date;
	}else{//root
	    duration = -1;
	}
	if(duration == 0){
	    duration = numeric_limits<double>::min();
	}
	durations->push_back(duration);
    }
}


int a_feasible_time(Node * node,double timeAnc){
    double minTime=0.0,maxTime=0.0;
    minTime = node->min; maxTime = node->max;
    if (node->free == true){
	if (node->maxb == true){
	    if (node->max < timeAnc){
		timeAnc=node->max;   /* the age of this node must be <= to its maxAge */
	    }
	}
	node->date = timeAnc - (timeAnc-minTime)*(0.02+myRand()*0.96)/log(node->order+3);
    }
    for(int i=0;i<node->getChildCount();i++){
	a_feasible_time(&node->getChild(i),node->date);
    }
    return 1;
}

double round(double x){
    double testx = x-int(x);
    if(testx>=0.5){return ceil(x);}else{return floor(x);}
}

float rand_float_range(float a, float b){
    return ((b-a)*((float)rand()/RAND_MAX))+a;
}

double myRand(void){
    return rand()/(double)RAND_MAX;
}

void set_node_order(Tree * tree){
    for(int i=0;i<tree->getExternalNodeCount();i++){
	Node * nd = tree->getExternalNode(i);
	int order = 0;
	nd->order = order;
	while(nd->hasParent()){
	    nd = nd->getParent();
	    order += 1;
	    if(nd->order < order)
		nd->order = order;
	}
    }
}

/*
 * which one is 
 * 0 tnc
 * 1 lbfgs
 * 2 tnewton
 * 3 mma
 */
int optimize_best(int whichone, bool ad, double * init_x, pl_calc_parallel * plp, int numiter,bool moredetail,double ftol,double xtol){
    if(whichone==0){
	return optimize_plcp_tnc(init_x,plp,numiter,ad,ftol);
    }else if(whichone > 0 && whichone < 6){
	return optimize_plcp_nlopt(init_x,plp,numiter,whichone,ad,moredetail,ftol,xtol);
    } 
    else{
	cout << "BAD OPTIMIZATION NUMBER" << endl;
	exit(0);
    }
    return 0;
}

//THESE ARE ALL AD
int optimize_best_parallel(int whichone, double * init_x, pl_calc_parallel * plp, int numiter, bool moredetail, double ftol,double xtol){
    if(whichone == 0){
	return optimize_plcp_tnc_ad_parallel(init_x,plp,ftol);
    }else if(whichone > 0 && whichone < 6){
	return optimize_plcp_nlopt_ad_parallel(init_x,plp,numiter,whichone, moredetail,ftol,xtol);
    }else{
	cout << "BAD OPTIMIZATION NUMBER" << endl;
	exit(0);
    }
    return 0;
}

//TODO: add timing to this, the faster one may be better too
void prime_optimization(pl_calc_parallel & plp, vector<double> & params){
    bool moredetail, moredetailad,moredetailcv;
    siman_calc_par smp3;
    smp3.set_verbose(false);
    smp3.set_pl(&plp,50,0.999,0.07,0.25,15000,0.001);
    smp3.optimize(params);
    cout << params[0] << endl;
    double exsim = plp.calc_pl(params);
    cout << "exit siman: " << exsim << endl;
    double ftol = 1e-7;
    double xtol = 1e-7;
    
    cout << "now priming regular" << endl;
    double bestsc = 0;
    int bestind = 0;
    for(int i=0;i<6;i++){
	double *x2 = new double[plp.numparams];
	for (unsigned int j = 0; j < params.size(); j++){x2[j]   = params[j];}
	int rc;
	rc = optimize_best(i,false,x2,&plp,10000,false,ftol,xtol);
	vector<double> nparams;
	for (unsigned int j = 0; j < params.size(); j++){nparams.push_back(x2[j]);}
	double out2 = plp.calc_pl(nparams);
	if((exsim-out2) > bestsc){
	    bestind = i;
	    bestsc = exsim-out2;
	    if (rc == 3 && i > 0)
		moredetail = true;
	    else
		moredetail = false;
	}
	cout << rc << " " << out2 << " " << exsim-out2<< endl;
	delete []x2;
    }
    cout << "----"<<endl;
    bestsc = 0;
    int bestindad = 0;
    cout << "now priming AD" << endl;
    for(int i=0;i<6;i++){
	double *x2 = new double[plp.numparams];
	for (unsigned int j = 0; j < params.size(); j++){x2[j]   = params[j];}
	int rc;
	rc = optimize_best(i,true,x2,&plp,10000,false,ftol,xtol);
	vector<double> nparams;
	for (unsigned int j = 0; j < params.size(); j++){nparams.push_back(x2[j]);}
	double out2 = plp.calc_pl(nparams);
	if((exsim-out2) > bestsc){
	    bestindad = i;
	    bestsc = exsim-out2;
	    if (rc == 3 && i > 0)
		moredetailad = true;
	    else
		moredetailad = false;
	}
	cout << rc << " " << out2 << " " << exsim-out2<< endl;
	delete []x2;
    }
    cout << "----"<<endl;
    bestsc = 0;
    int bestindcv = 0;
    cout << "now priming CV (AD)"<< endl;
    //get the one node to CV, last node should always be a tip
    vector<int> cvnodes(plp.free->size(),0);
    cvnodes[cvnodes.size()-1] = 1;
    vector<int> freeparamscv;
    vector<int> free(plp.free->begin(),plp.free->end());
    plp.set_cv_nodes(cvnodes);
    int numparamscv = generate_param_order_vector(&freeparamscv, false, &cvnodes, &free);
    vector<double> cvplparams;
    plp.set_freeparams(numparamscv,false,&freeparamscv, &cvplparams);
    for(int i=0;i<6;i++){
	double *x2 = new double[plp.numparams];
	for (unsigned int j = 0; j < cvplparams.size(); j++){x2[j]   = cvplparams[j];}
	int rc;
	rc = optimize_best_parallel(i,x2,&plp,10000,false,ftol,xtol);
	vector<double> nparams;
	for (unsigned int j = 0; j < cvplparams.size(); j++){nparams.push_back(x2[j]);}
	double out2 = plp.calc_pl(nparams);
	if((exsim-out2) > bestsc){
	    bestindcv = i;
	    bestsc = exsim-out2;
	    if (rc == 3 && i > 0)
		moredetailcv = true;
	    else
		moredetailcv = false;
	}
	cout << rc << " " << out2 << " " << exsim-out2<< endl;
	delete []x2;
    }

    cout << "best: " << bestind  << "(" << moredetail << ") bestad: " << bestindad << " (" << moredetailad << ")"<< 
	" bestcv: " << bestindcv << " (" << moredetailcv << ")" << endl; 
    cout << "PLACE THE LINES BELOW IN THE CONFIGURATION FILE"<< endl;
    cout << "opt = " << bestind << endl; //TODO ITERATIONS
    if (moredetail)
	cout << "moredetail" << endl;
    cout << "optad = " << bestindad << endl;//TODO ITERATIONS
    if (moredetailad)
	cout << "moredetailad" << endl;
    cout << "optcvad = " << bestindcv << endl;
    if (moredetailcv)
	cout << "moredetailcvad" << endl;

}

bool check_possible_cv_nodes(const int cvnode, const vector<int> & samps, const vector<int> & parent_nds_ints){
    int cvpar = parent_nds_ints[cvpar];
    for(int i=0;i<samps.size();i++){
	if(parent_nds_ints[samps[i]] == cvpar)
	    return false;
    }
    return true;
}

double get_median(vector<double> & container){
    std::vector<double>::iterator first = container.begin();
    std::vector<double>::iterator last = container.end();
    std::vector<double>::iterator middle = first + (last - first) / 2;
    std::nth_element(first, middle, last); // can specify comparator as optional 4th arg
    double median = *middle;
    return median;
}

double get_gmean(vector<double> & container){
    double product = 1;
    for(int i=0;i<container.size();i++){
	product *= container[i];
    }
    return pow(product,1/(double)container.size());
}

double get_sum(vector<double> & container){
    double s = 0;
    for(int i=0;i<container.size();i++){
	s += container[i];
    }
    return s;
}

void process_ind8s_inr8s(string ind8s,string inr8s,vector<double> * params, Tree * intree){
    cout << "calculating with input parameters in files " << ind8s << " and " << inr8s << endl;
    ifstream infile(ind8s.c_str());
    ifstream infile2(inr8s.c_str());
    if (!infile || !infile2){
	cerr << "Could not open file." << endl;
	exit(0);
    }
    string line;
    Tree *dtree = NULL;
    Tree *rtree = NULL;
    while (getline(infile, line)){
	if (line.size() > 1){
	    TreeReader tr;
	    dtree = tr.readTree(line);
	    dtree->setHeightFromTipToNodes();
	    cout << "nodes in dates tree: " << dtree->getNodeCount() << endl;
	}
    }while (getline(infile2, line)){
	if (line.size() > 1){
	    TreeReader tr;
	    rtree = tr.readTree(line);
	    cout << "nodes in rates tree: " << rtree->getNodeCount() << endl;
	}
    }
    int pc = 0;
    int numsites = 2744;
    for(int i=0;i<intree->getNodeCount();i++){
	if(i == 0){
	    continue;
	}else{
	    cout << intree->getNodeByNodeNumber(i)->getNumber() << endl;
	    vector<string> names;
	    intree->getNodeByNodeNumber(i)->getExternalNodeNames(&names);
	    params->at(pc) = rtree->getMRCA(names)->getBL()*numsites;
	    pc++;
	}
    }
    for(int i=0;i<intree->getNodeCount();i++){
	if(intree->getNodeByNodeNumber(i)->getChildCount() == 0 ||
	   intree->getNodeByNodeNumber(i)->free == false){
	    continue;
	}else{
	    cout << intree->getNodeByNodeNumber(i)->getNumber() << endl;
	    vector<string> names;
	    intree->getNodeByNodeNumber(i)->getExternalNodeNames(&names);
	    params->at(pc) = dtree->getMRCA(names)->getHeight();
	    pc++;
	}
    }
    cout << params->size() << " " << pc << endl;
}


/*
 * This performs a full optimization round on the tree
 *
 * This includes both analytical deriv and auto deriv as well 
 * as fall back no deriv and simulated annealing
 * 
 * If parallel is set, this will use the adolc instead of the rad
 * for auto diff
 *
 */
//TODO: if exit from optimization step is worse than entrance, switch back
bool optimize_full(pl_calc_parallel & plp, vector<double> * params, const OptimOptions * optims, bool cv){
    bool madeit;
    int numiter = optims->pliter;
    if(plp.lf)
	numiter = optims->lfiter;
    else if (plp.calc_isum((*plp.get_cv_nodes())) > 0)
	numiter = optims->cviter;

    int totaliters = 0;
    for(int i=0;i<numiter;i++){
	madeit = true;
	//start with siman calc
	siman_calc_par smp3;
	smp3.set_verbose(optims->verbose);
	smp3.set_pl(&plp,50/float(i+1),0.999,0.07,0.25,optims->lfsimaniter,0.001);
	smp3.optimize(*params);

	double exsim = plp.calc_pl(*params);
	cout << "exit siman: " << exsim << endl;

	double *x2 = new double[plp.numparams];
	for (unsigned int j = 0; j < params->size(); j++){x2[j]   = (*params)[j];}
	int rc;
	rc = optimize_best(optims->bestopt,false,x2,&plp,optims->maxoptimiters,optims->moredetail,optims->ftol,optims->xtol);
	//check here for success or not
	for (unsigned int j = 0; j < params->size(); j++){(*params)[j] = x2[j];}	
	delete[] x2;
	//need a catch here is value is == LARGE
	double out1 = plp.calc_pl(*params);
	if(out1 == optims->LARGE){
	    cout << "problem with params" << endl;
	    for (unsigned int j = 0; j < params->size(); j++){cout << "p: " << (*params)[j] << endl;}
	    exit(0);
	}
	if((exsim - out1) > 0.0001){
	    madeit = false;
	}else if (exsim == out1 && rc < 1){
	    cout << "calculating without gradient (might want to try a different opt=VALUE)" << endl;
	    x2 = new double[plp.numparams];
	    for (unsigned int j = 0; j < params->size(); j++){x2[j]   = (*params)[j];}
	    rc = optimize_best(5,false,x2,&plp,optims->maxoptimiters,optims->moredetail,optims->ftol,optims->xtol);
	    for (unsigned int j = 0; j < params->size(); j++){(*params)[j] = x2[j];}
	    delete[] x2;
	    out1 = plp.calc_pl(*params);
	}
	cout << "after opt calc1: " << out1 << endl;
	x2 = new double[plp.numparams];
	for (unsigned int j = 0; j < params->size(); j++){x2[j]   = (*params)[j];}
	if (cv == false){
	    if(plp.numparams < 100000)
	        rc = optimize_best(optims->bestadopt,true,x2,&plp,optims->maxoptimiters,optims->moredetailad,optims->ftol,optims->xtol);
            else
		rc = optimize_best(optims->bestadopt,false,x2,&plp,optims->maxoptimiters,optims->moredetailad,optims->ftol,optims->xtol);
	}else{
	    rc = optimize_best_parallel(optims->bestadopt,x2,&plp,10000,optims->moredetailad,optims->ftol,optims->xtol);
	}
	for (unsigned int j = 0; j < params->size(); j++){(*params)[j] = x2[j];}
	delete[] x2;
	double out2 = plp.calc_pl(*params);
	if(out2 == optims->LARGE){
	    cout << "problem with params" << endl;
	    for (unsigned int j = 0; j < params->size(); j++){cout << "p: " << (*params)[j] << endl;}
	    exit(0);
	}
	if((out1-out2) > 0.0001){
	    madeit = false;
	}else if (out2 == out1 && rc < 1){
	    cout << "calculating without gradient (might want to try a different optad=VALUE)" << endl;
	    x2 = new double[plp.numparams];
	    for (unsigned int j = 0; j < params->size(); j++){x2[j]   = (*params)[j];}
	    rc = optimize_best(5,false,x2,&plp,optims->maxoptimiters,optims->moredetail,optims->ftol,optims->xtol);
	    for (unsigned int j = 0; j < params->size(); j++){(*params)[j] = x2[j];}
	    delete[] x2;
	    out2 = plp.calc_pl(*params);
	}
	cout << "after opt calc2: " << out2 << endl;
	totaliters += 1;
	if (optims->thorough == true && (i+1) == optims->lfiter){
	    if(totaliters >= 1000){//just last ditch
		cout << "thorough optimization hit 1000 iterations, breaking" << endl;
		break;
	    }
	    if(madeit == true)
		break;
	    else{
		i=0;
	    }
	}
    }
    return madeit;
}

void mapspace(pl_calc_parallel &plp, vector<double> * params){
    siman_calc_par smp3;
    smp3.mapspace(&plp,*params);
}


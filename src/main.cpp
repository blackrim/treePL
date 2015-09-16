/*
 * main.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <fenv.h>
#include <string.h>
#include <algorithm> 
#include <limits>
#include <omp.h>
#include <adolc/adolc_openmp.h>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "optim_options.h"
#include "pl_calc_parallel.h"
#include "optimize_tnc.h"
#include "optimize_nlopt.h"
#include "siman_calc_par.h"

#define LARGE 1e+15

// line endings
#define LINETERM_UNIX  0
#define LINETERM_MAC   1
#define LINETERM_DOS   2

int main(int argc,char* argv[]) {
    int seed = -1;
    cout.precision(8);
    int retval;
    retval = -1;
    retval = feraiseexcept( FE_ALL_EXCEPT );
    bool setparams = false;
    if(argc != 2 || string(argv[1]) == "-h" || string(argv[1]) == "--help"){
        cout << "treePL version 1.0"<< endl;
        cout << "you need more arguments." << endl;
        cout << "usage: treePL configfile" << endl;
        exit(0);
    }
    string treefile;
    map<string,vector<string> > mrcas;
    map<string,double > mrca_mins;
    map<string,double > mrca_maxs;
    bool cv = false;
    bool collapse = false;
    double smooth = 10;
    double sample = 1.;
    int numsites =0;
    double scale=1;
    bool checkconstraints = false;
    string outfilen = "";
    string cvoutfile = "cv.out";
    OptimOptions oopt;

    double cvstart = 1000;
    double cvstop = 0.1; 
    double cvmultstep = 0.1;
    bool randomcv = false;
    int randomcviter = 10;
    double randomcvsamp = 0.1;
    bool verbose = false;
    bool paramverbose = false;
    bool prime = false;
    bool mapspaceb = false;
    bool log_pen = false;

    string ind8s;
    string inr8s;
    string configFileName = argv[1];
    
    read_config(configFileName, cvstart, cvstop, cvmultstep, randomcv, randomcviter, 
        randomcvsamp, verbose, paramverbose, prime, mapspaceb, log_pen, treefile, 
        cv, collapse, smooth, sample, numsites, scale, checkconstraints, outfilen, 
        cvoutfile, mrcas, mrca_mins, mrca_maxs, oopt, seed, ind8s, inr8s);
    

    /*
     * start the analyses
     */

    //added loop to do multiple trees in the file
    ifstream infile(treefile.c_str());
    if(!infile){
        cerr << "Could not open file." << endl;
        return 1;
    }
    ofstream outFile;
    outFile.open(outfilen.c_str(),ios::out);
    ofstream outFile2;
    string routfilen = outfilen+".r8s";
    outFile2.open(routfilen.c_str(),ios::out);
    ofstream paramverboseFile;

    // this should read each tree in a file and process those
    //while (getline(infile, line)){
    string line;
    while (getlineSafe(infile, line)) { // version of getline where line endings are unknown
        if(line.size() > 1){
            TreeReader tr;
            Tree *tree = NULL;
            tree = tr.readTree(line);
            /*
             * change 0 BL to 1/numsites
             */
            process_initial_branch_lengths(tree,collapse,numsites);
            tree->setHeightFromTipToNodes();
            //TODO: set the max iters higher if bigger tree
            /*
             * setting up min and max constraints
             */
            map<Node *,double> mins;map<Node *,double> maxs;
            map<Node *,string> min_names; map<Node *,string> max_names;
            map<string,double>::iterator it;
            for (it=mrca_mins.begin(); it != mrca_mins.end(); it++ ){
                if(mrcas.count((*it).first) == 0){
                    cout << "problem with mrca " << (*it).first << " (probably no mrca named this)" <<endl;
                    exit(0);
                }
                for(int m=0;m<mrcas[(*it).first].size();m++){
                    Node * tn = tree->getExternalNode(mrcas[(*it).first][m]);
                    if (tn == NULL){
                        cout << "problem with taxa for mrca " << (*it).first << " = " << mrcas[(*it).first][m] << " doesn't exist " << endl;
                        exit(0);
                    }
                }
                Node * tmrca =tree->getMRCA(mrcas[(*it).first]);
                if(tmrca == NULL){
                    cout << "problem with mrca " << (*it).first << " (probably bad names to nodes)" << endl;
                    exit(0);
                }
                mins[tmrca] = (*it).second/scale;
                min_names[tmrca] = (*it).first;
                tmrca->min = (*it).second/scale;
                tmrca->minb = true;
                tmrca->pen_minb = true;
                tmrca->pen_min = (*it).second/scale;
                cout << "setting " << (*it).first << " min: " << (*it).second << endl;
            }
            for (it=mrca_maxs.begin(); it != mrca_maxs.end(); it++ ){
                if(mrcas.count((*it).first) == 0){
                    cout << "problem with mrca " << (*it).first << " (probably no mrca named this)" <<endl;
                    exit(0);
                }
                Node * tmrca =tree->getMRCA(mrcas[(*it).first]);
                if(tmrca == NULL){
                    cout << "problem with mrca " << (*it).first << " (probably bad names to nodes)" << endl;
                    exit(0);
                }
                maxs[tmrca] = (*it).second/scale;
                max_names[tmrca] = (*it).first;
                tmrca->max = (*it).second/scale;
                tmrca->maxb = true;
                if(tmrca->minb == true){
                    if(tmrca->min != tmrca->max){//don't set the penalty if it is a fixed age
                        tmrca->pen_maxb = true;
                        tmrca->pen_max = (*it).second/scale;
                    }else{
                        tmrca->pen_minb = false;
                        tmrca->pen_max = 0;
                    }
                }else{
                    tmrca->pen_maxb = true;
                    tmrca->pen_max = (*it).second/scale;
                }
                cout << "setting " << (*it).first << " max: " << (*it).second << endl;
            }

            /*
             * checking that all the constraints are valid 
             * by asking whether only younger ones are children
             * of older ones
             * can basically just do this with maxes
             */
            if(checkconstraints == true){
                map<Node *,double>::iterator it;
                
            /*  
                bool max = false;
                for (it=maxs.begin(); it != maxs.end(); it++ ){
                    Node * tempn = (*it).first;
                    if(tempn != tree->getRoot()){
                        tempn = tempn->getParent();
                    }else{
                        continue;
                    }
                    double curtest = (*it).second;
                    while (tempn != tree->getRoot()){
                        if(maxs.count(tempn) == 1){
                            if(mins[tempn]  curtest){
                                cout << min_names[tempn] << " in conflict " << min_names[(*it).first] << endl;
                                exit(0);
                            }
                        }
                        tempn = tempn->getParent();
                    }
                }
                if(max == true)
                exit(0);
            */

                bool min = false;
                for (it=mins.begin(); it != mins.end(); it++ ){
                    Node * tempn = (*it).first;
                    if(tempn != tree->getRoot()){
                        tempn = tempn->getParent();
                    }else{
                        continue;
                    }
                    double curtest = (*it).second;
                    while (tempn != tree->getRoot()){
                        if(mins.count(tempn) == 1){
                            if(mins[tempn] < curtest){
                                cout << min_names[tempn] << " in conflict " << min_names[(*it).first] << endl;
                                exit(0);
                            }
                        }
                        tempn = tempn->getParent();
                    }
                }
                if(min == true){
                    exit(0);
                }
            }
            /*
             * END CHECKING THE CONSTRAINTS
             */
            /****************
             * BEGIN REAL ANALYSES
             **************/
            //JUST DO PL_CALC
            //setup the data structures
            //give the tree preorder numbering
            cout << "preorder prep" << endl;
            apply_node_label_preorder(tree);
            cout << "calculating character durations" << endl;
            calc_char_durations(tree, numsites);
            cout << "setting min and max" << endl;
            set_mins_maxs(tree,&mins,&maxs);
            cout << "setting up all constraints" << endl;
            setup_date_constraints(tree,&mins,&maxs);
            vector<int> free;
            vector<double> vmin;
            vector<double> vmax;
            vector<double> penmin;
            vector<double> penmax;
            vector<double> char_durations;
            vector<double> log_fact_char_durations;
            vector<int> parent_nds_ints;
            vector<int> child_counts;
            vector< vector<int> > children_vec;//this is to speed up some operations
            extract_tree_info(tree,&free,&parent_nds_ints,&child_counts, &char_durations,
                      &log_fact_char_durations,&vmin,&vmax,&children_vec,&penmin,&penmax);
            vector<double> start_dates;
            vector<double> start_rates;
            vector<double> start_durations;
            set_node_order(tree);
            cout << "getting feasible start dates" <<endl;
            get_feasible_start_dates(tree,&start_dates,&start_rates,&start_durations);
            double start_rate = get_start_rate(tree,&start_durations);//numsites/20.;
            cout << "start rate " << start_rate << endl;
            start_rates[1] = start_rate;
            double minrate = 0.0;
            /*
             * Start with Langley Fitch (one rate)
             */
            pl_calc_parallel plp;
            plp.setup_starting_bits(&parent_nds_ints,&child_counts, &free, 
                        &char_durations, &log_fact_char_durations, &vmin, 
                        &vmax, &start_dates, &start_rates, &start_durations);
            plp.set_log_pen(log_pen);
            plp.minrate = minrate;
            plp.pen_min = &penmin;
            plp.pen_max = &penmax;
            plp.children_vec = &children_vec;
            //delete current file
            if(paramverbose){
                paramverboseFile.open("paramverboself",ios::out);
            }
            plp.paramverbose = paramverbose;
        
            vector<int> freeparams;
            //lf
            int numparams = generate_param_order_vector(&freeparams, true, NULL, &free);
        
            cout << "numparams:" << numparams << endl;
            vector<double> params;
        
            plp.set_freeparams(numparams, true, &freeparams, &params);
            double initcalc = plp.calc_pl(params);
            cout << "initial calc: " << initcalc << endl;

            // TODO: if initialization fails, try again until no failure.
            if(isnan(initcalc) || isinf(initcalc) || initcalc == LARGE){
                bool good = false;
                int trycount = 0;
                while (!good) {
                    cout << "problem initializing. trying again." << endl;
                    cout << "attempting to get feasible start rates/dates." <<endl;
            // clear vectors of bad values
                    //start_dates.clear();
                    //start_rates.clear();
                    //start_durations.clear();
                    get_feasible_start_dates(tree,&start_dates,&start_rates,&start_durations);
                    start_rate = get_start_rate(tree,&start_durations);
                    cout << "new start rate " << start_rate << endl;
                    start_rates[1] = start_rate;
                    minrate = 0.0;
                    pl_calc_parallel plp;
                    plp.setup_starting_bits(&parent_nds_ints,&child_counts, &free,
                                            &char_durations, &log_fact_char_durations, &vmin,
                                            &vmax, &start_dates, &start_rates, &start_durations);
                    plp.set_log_pen(log_pen);
                    plp.minrate = minrate;
                    plp.pen_min = &penmin;
                    plp.pen_max = &penmax;
                    plp.children_vec = &children_vec;

                    //freeparams.clear();
                    numparams = generate_param_order_vector(&freeparams, true, NULL, &free);
                    plp.set_freeparams(numparams, true, &freeparams, &params);
                    initcalc = plp.calc_pl(params);
                    cout << "initial calc: " << initcalc << endl;

                    if(isnan(initcalc) || isinf(initcalc) || initcalc == LARGE){
                        trycount += 1;
                        if(trycount == 10){ // 10 is a magic number here. just don't want infinite looping if pathological problem.
                            cout << "Failed setting feasible start rates/dates after 10 attempts. Aborting." << endl;
                            exit(0);
                        }
                    }else{
                        cout << "feasible start rates/dates found" << endl;
                        good = true;
                    }
                    //exit(0);
                }
            }
    /*        vector<double> g(params.size());
            cout << "calculating gradient" <<endl;
            g= vector<double>(params.size());
            plp.calc_gradient(params,&g);
            for(int j=0;j<params.size();j++){
                cout << params[j] << "\t" << g[j] << endl;
            }
            cout << "calculating gradient ad" <<endl;
            g=vector<double>(params.size());
            cout << "f: "<< plp.calc_function_gradient(&params,&g) << endl;
            for(int j=0;j<params.size();j++){
                cout << params[j] << "\t" << g[j] << endl;
            }
    */
            //prime optimization
            if(prime == true){
                prime_optimization(plp,params);
                exit(0);
            }

            //optimize
            //this is for testing thorough, has to stay made it for the entire iteration
            bool madeit = optimize_full(plp,&params,&oopt,false);
            cout << "exited lf converged :" << madeit<< endl;
            initcalc = plp.calc_pl(params);
            cout << "lf calc: " << initcalc << endl;
            if(paramverbose){
                paramverboseFile.close();
                paramverboseFile.open("paramverbosepl",ios::out);
            }
    //        exit(0);

            /****
             * CROSS VALIDATION IF REQUIRED
             */
            if(cv == true){
                ofstream cvoutFl;
                cvoutFl.open(cvoutfile.c_str(),ios::out);
                cout << "conducting cross validation analysis" << endl;
                int cvtype = 0; //0 = LOOCV, 1 = RANDOM
                if(randomcv == true){
                    cvtype = 1;
                    cout << "conducting RANDOM SUBSAMPLE cross validation analysis" << endl;
                }
                double curcv = cvstart;
                bool failed = false;
                double lowest_smoothing = 1000;
                double lowest_chi = 100000000;
                //get sampling groups
                vector< vector<int> > samp_groups;
                if(cvtype == 0){
                    for (int i = 0;i < tree->getExternalNodeCount(); i++){
                        vector<int> samps;
                        samps.push_back(tree->getExternalNode(i)->getNumber());
                        samp_groups.push_back(samps);
                    }
                }else if(cvtype == 1){//kfold cv, TODO fix this for new array types
                    //srand(time(NULL)); // this was used twice! don't do that!
                    int sampsize = (tree->getExternalNodeCount()*randomcvsamp)+0.5;//round up
                    /* uncomment for sample without replacement
                      vector<int> poss_nodes;
                    for(int i=0;i<tree->getExternalNodeCount();i++){
                    poss_nodes.push_back(tree->getExternalNode(i)->getNumber());
                    }*/
                    int numtosamp = randomcviter;
                    cout << "sampling " << sampsize << " tips from " << tree->getExternalNodeCount() << " total tips " << randomcviter << " times" << endl;
                    for(int i=0;i<randomcviter;i++){
                        //sampling with replacement
                        vector<int> poss_nodes;
                        for(int j=0;j<tree->getExternalNodeCount();j++){
                            if(tree->getExternalNode(j)->getParent() != tree->getRoot())
                            poss_nodes.push_back(tree->getExternalNode(j)->getNumber());
                        } 
                        vector<int> samps;
                        int value; 
                        int cvfail = 0;
                        bool completefail = false;
                        vector<int>::iterator it;
                        int numtosamp = sampsize;
                        while(numtosamp > 0){
                            value = (int)rand() % poss_nodes.size();
                            //check poss_nodes[value] is not a sister to anything already in samps
                            bool checkcv = check_possible_cv_nodes(poss_nodes[value],samps,parent_nds_ints);
                            if(checkcv == false){
                                cout << "checkcv FALSE" << endl; 
                                cvfail += 1;
                                if(cvfail > 1000){
                                    completefail = true;
                                    break;
                                }else{
                                    continue;
                                }
                            }
                            samps.push_back(poss_nodes[value]);
                            it = poss_nodes.begin()+value;
                            poss_nodes.erase(it);
                            numtosamp -= 1;
                        }
                        samp_groups.push_back(samps);
                        if(completefail == true){
                            cout <<"complete failure at cv for not getting sister species" << endl;
                            exit(0);
                            i--;//make sure this works, should be a singleton category
                        }
                    }
                }
        //      exit(0); // TODO: make sure the groups are broken up correctly

                //start cv
                cout << "--------------"<<endl;
                vector<double> cvstrates(start_rates.size(),params[0]); 
                vector<double> cvstdates(plp.dates);
                vector<double> cvstdur(plp.durations);
                bool startopt = true;
        
                while(curcv >= cvstop){
                    if(startopt){
                        freeparams.clear();
                        numparams = generate_param_order_vector(&freeparams, false, NULL, &free);
                        cout << "numparams:" << numparams << endl;
                        plp.smoothing = curcv;
                        cout << "smoothing:" << plp.smoothing << endl;
                        plp.set_freeparams(numparams, false, &freeparams, &params);
                        cout << plp.calc_pl(params) << endl;
                        optimize_full(plp,&params,&oopt,false);
                        cout << "after opt calc: " << plp.calc_pl(params) << endl;
                        cvstdates = vector<double>(plp.dates);
                        cvstdur = vector<double>(plp.durations);
                        cvstrates = vector<double>(plp.rates);
                    }
                    double chisq = 0;
                    cout << "curcv: "<< curcv << endl;
                    int chisqcount = 0;
                    int i;
                    vector<double> sqerrs;
                //was shared(chisq,chisqcount)
                #pragma omp parallel for ADOLC_OPENMP_NC reduction(+:chisq) reduction(+:chisqcount) num_threads(oopt.nthreads)
                    for(i=0;i<samp_groups.size();i++){
                        //cout << "thread_num: "<< omp_get_thread_num() << endl;
                        //cout << "num threads: " << omp_get_num_threads() << endl;
                        vector<int> cvnodes(free.size(),0);
            //            cout << "cv vec: ";
                        for(int j=0;j<samp_groups[i].size();j++){
                            cvnodes[samp_groups[i][j]]=1;
            //                cout << samp_groups[i][j] << " " ;
            //                cout << tree->getNodeByNodeNumber(samp_groups[i][j])->getName() << " " ;
                        }
            //            cout << endl;
                        pl_calc_parallel plpcv;
                        plpcv.setup_starting_bits(&parent_nds_ints,&child_counts, &free, 
                                &char_durations, &log_fact_char_durations, &vmin, 
                                &vmax, &cvstdates, &cvstrates, &cvstdur);
                        plpcv.set_log_pen(log_pen);
                        plpcv.children_vec = &children_vec;
                        plpcv.pen_min = &penmin;
                        plpcv.pen_max = &penmax;
                        plpcv.smoothing = curcv;
                        plpcv.set_cv_nodes(cvnodes);
                        vector<int> freeparamscv;        
                        int numparamscv = generate_param_order_vector(&freeparamscv, false, &cvnodes, &free);
                        //cout << "numparams:" << numparamscv << endl;
                        vector<double> cvplparams;
                        plpcv.set_freeparams(numparamscv, false, &freeparamscv, &cvplparams);
                        cout << "cv check: " << plpcv.calc_pl(cvplparams) << endl;
                        //cout << "samp_group: "<<samp_groups[i].size() << endl;

                        //add a loop to remove the cv nodes for fold
                        optimize_full(plpcv,&cvplparams,&oopt,true);
                        double poc = plpcv.calc_pl(cvplparams);
                        cout << "post opt cv check: " << poc << endl;
                        int fcount = 0;
                        /*for(int j=0;j<plpcv.rates.size();j++){
                            if(plpcv.freeparams[fcount] != -1){
                                cout << "r:\t" <<cvplparams[plpcv.freeparams[fcount]]/numsites <<endl;
                            }
                            fcount += 1;
                        }
                        for(int j=0;j<plpcv.numnodes;j++){
                            if(plpcv.freeparams[fcount] != -1){
                                cout  <<"d:\t" <<cvplparams[plpcv.freeparams[fcount]]/numsites<< endl;
                            }
                            fcount += 1;
                        }*/
                        //TODO: return 0 if length is 0
                        if(cvtype==0){//loocv chisq
                            double ratees = 0;
                            int par = parent_nds_ints[samp_groups[i][0]];
                            if(par==0){//is the root
                                for(int j=0;j<parent_nds_ints.size();j++){
                                    if(parent_nds_ints[j] == par && j != samp_groups[i][0]){
                                        ratees = plpcv.rates[j];
                                    }
                                }
                            }else{
                                ratees = plpcv.rates[par];
                            }
                            double d = plpcv.durations[samp_groups[i][0]];
                            double expe = ratees*d;
                            double length = char_durations[samp_groups[i][0]];
                            double sq =(length-expe)*(length-expe);// (expe-length)*(expe-length);
                            if(d < 1e-100 || expe < 1e-100){
                                chisq += 0;
                            }else{
                                chisq += sq/expe;//average chisq
                            }
                            //sqerrs.push_back(sq);
                            cout << tree->getNodeByNodeNumber(samp_groups[i][0])->getName()<<  " rate "<< ratees << " expe "<< expe << " dur " << d << " len " << length << " sq " << sq << " chisq " << chisq << endl;                 
                            //name obj 
            //                cout << tree->getNodeByNodeNumber(samp_groups[i][0])->getName() << "\t" << poc << "\t" << sq<< "\t" << chisq << "\t" << ratees << "\t" << d;
            //                for(int j=0;j<plpcv.rates.size();j++){
                //                cout << "\t" << plpcv.rates[j]/float(numsites);
            //                }
            //                for(int j=0;j<plpcv.dates.size();j++){
                //                cout << "\t" << plpcv.dates[j]/float(numsites);
            //                }
            //                cout << endl;
                        }else{//randomcv chisq
                            double tchisq = 0;
                            // for(int k = 0; k < nds.size(); k++){
                            for(int k=0; k < samp_groups[i].size(); k++){
                                double ratees = 0;
                                int par = parent_nds_ints[samp_groups[i][k]];
                                if(par==0){
                                    for(int j=0; j < parent_nds_ints.size(); j++){
                                       if(parent_nds_ints[j] == par && j != samp_groups[i][k]){
                                           ratees = plpcv.rates[j];
                                       }
                                    }
                                }else{
                                    ratees = plpcv.rates[par];
                                }
                                double d = plpcv.durations[samp_groups[i][k]];
                                double expe = ratees*d;
                                double length = char_durations[samp_groups[i][k]];
                                double sq = (expe-length)*(expe-length);
                                if(d < 1e-100 || expe < 1e-100){
                                    tchisq += 0;
                                }else{
                                    tchisq += sq/expe;//average chisq
                //                cout << "chi: " << d << " " << expe << " " << length << " " << sq << endl;
                //                cout << "tchisq: " << tchisq << endl;
                                }
                            }
                            chisq += tchisq/(float)samp_groups[i].size();
                            cout << "(chisq): " << chisq << endl;
                        }
               
                        plpcv.delete_ad_arrays();
                        chisqcount += 1;
                        cout << "cv: " << i+1 << endl;
            //            exit(0);
                    }
            
                    cout << "chisq: (" << curcv << ") "<< chisq << endl;//"\tsq:" << get_sum(sqerrs) << "\tmed: " << get_gmean(sqerrs) << endl;
                    cvoutFl << "chisq: ("<< curcv << ") " << chisq << endl;
                    if(std::numeric_limits<double>::infinity() == chisq){
                        cout << "INF"<<endl;
                        exit(0);
                    }
                    if(chisq < lowest_chi){
                        lowest_chi = chisq;
                        lowest_smoothing = curcv;
                    }
                    curcv = curcv * cvmultstep;
                    cout << "--------------"<< endl;
                    //exit(0);
                }
                cout << "check:"<<plp.calc_pl(params) <<endl;
                plp.smoothing = lowest_smoothing;
                cvoutFl.close();
            }else{
                plp.smoothing = smooth;
            }
            //exit(0);
            /*
             * END CROSS VALIDATION IF REQUIRED
             ****/
            //Now with Penalized Likelihood
            freeparams.clear();
            numparams = generate_param_order_vector(&freeparams, false, NULL, &free);
            cout << "numparams:" << numparams << endl;
            cout << "smoothing:" << plp.smoothing << endl;
            plp.set_freeparams(numparams, false, &freeparams, &params);

            //JUST FOR TESTING R8S things
            if(ind8s.size() > 1 && inr8s.size() > 1){
                process_ind8s_inr8s(ind8s,inr8s,&params, tree);
                cout << plp.calc_pl(params) << endl;
                exit(0);
            }
            cout << plp.calc_pl(params) << endl;
    //        exit(0);
            optimize_full(plp,&params,&oopt,false);
            cout << "after opt calc: " << plp.calc_pl(params) << endl;
            /*
             * write the final dates to a text file
             */
    /*
            for(int i=0;i<params.size();i++){
            cout <<"p: " << params[i]/numsites << endl;
            }
            cout << "calculating gradient ad" <<endl;
            g= vector<double>(params.size());
            cout << "f: "<< plp.calc_function_gradient(&params,&g) << endl;
            for(int j=0;j<params.size();j++)
                cout << params[j]/numsites << "\t" << g[j] << endl;

            cout << "calculating gradient par" <<endl;
            g= vector<double>(params.size());
            plp.calc_pl_function_gradient_adolc(&params,&g);
            for(int j=0;j<params.size();j++)
                cout << params[j] << "\t" << g[j] << endl;
    */

            for(int i=0;i<tree->getNodeCount();i++){
                std::cout.precision(10);
        //        cout << "d " << plp.durations[tree->getNode(i)->getNumber()] << " " <<plp.dates[tree->getNode(i)->getNumber()] << endl;
                if(plp.durations[tree->getNode(i)->getNumber()] != -1){
                    tree->getNode(i)->setBL(plp.durations[tree->getNode(i)->getNumber()]);
                }
                //tree->getNode(i)->setBL(plp.durations[tree->getNode(i)->getNumber()]*scale);
            }
            outFile << tree->getRoot()->getNewick(true) <<";" << endl;

        
            for(int i=0;i<tree->getNodeCount();i++){
                tree->getNode(i)->setBL(plp.rates[tree->getNode(i)->getNumber()]/numsites);
            }
            outFile2 << tree->getRoot()->getNewick(true) << ";"<<endl;

            if(mapspaceb){
                mapspace(plp,&params);
            }

            delete tree;
            if(paramverbose){
                paramverboseFile.close();
            }
        }
    }
    outFile.close();
    outFile2.close();
    infile.close();

    return 0;
}

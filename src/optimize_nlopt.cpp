#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "optimize_nlopt.h"
#include "pl_calc_parallel.h"
#include <math.h>
#include <nlopt.h>
#include <iostream>
#include <stdio.h>
#include <limits>
#include <cmath>

typedef struct{
    vector <double> tncparams;
    vector <double> tncgrads;
    pl_calc_parallel * pl_;
} data_obj_nlopt;

#define LARGE 1e+15
#define SMALL 1e-20

double function_plcp_nlopt_ad(unsigned n, const double *x, double *g, void *state){
    data_obj_nlopt *d = (data_obj_nlopt *) state;
    for(unsigned int i=0;i<d->tncparams.size();i++){
	d->tncparams[i] = x[i];
    }

    double f = d->pl_->calc_function_gradient(&d->tncparams,&d->tncgrads);
//       cout << "adf: " << f << endl;
    if (isnan(f) || f == std::numeric_limits<double>::infinity()){
	return LARGE;
//    cout << "adf: " << f << endl;
    }

    int tv = d->tncgrads.size()-1;
    for(unsigned int i=0;i<d->tncgrads.size();i++){
	g[i] = d->tncgrads[i];
    }
    return f;
}

double function_plcp_nlopt(unsigned n, const double *x, double *g, void *state){
    data_obj_nlopt *d = (data_obj_nlopt *) state;
    for(unsigned int i=0;i<d->tncparams.size();i++){
	d->tncparams[i] = x[i];
//	cout <<"x " << x[i] << endl;
    }
//    cout << "here" <<endl;
    double f = d->pl_->calc_pl(d->tncparams);
    //   cout << "f: " << f << endl;
    if (isnan(f) || f == std::numeric_limits<double>::infinity()){
	cout << "throwing" << endl;
	return LARGE;
    }else{
	d->pl_->calc_gradient(d->tncparams,&d->tncgrads);
//	cout << "f2 " << f2 << endl;
	for(unsigned int i=0;i<d->tncgrads.size();i++){
	    g[i] = d->tncgrads[i];
//	    cout << "g: " << g[i] << " " << tgrads[i] << endl;
	}
    }
    return f;
}

double function_plcp_nlopt_nograd(unsigned n, const double *x, double *g, void *state){
    data_obj_nlopt *d = (data_obj_nlopt *) state;
    for(unsigned int i=0;i<d->tncparams.size();i++){
	d->tncparams[i] = x[i];
    }
    double f = d->pl_->calc_pl(d->tncparams);
//    cout << "nog: " << f << endl;
    if (isnan(f) || f == std::numeric_limits<double>::infinity()){
	cout << "throwing" << endl;
	return LARGE;
    }
    return f;
}

double function_plcp_nlopt_ad_parallel(unsigned n, const double *x, double *g, void *state){
    data_obj_nlopt *d = (data_obj_nlopt *) state;
    for(unsigned int i=0;i<d->tncparams.size();i++){
	d->tncparams[i] = x[i];
    }
    double f = d->pl_->calc_pl_function_gradient_adolc(&d->tncparams,&d->tncgrads);
    //   cout << "par: " << f << endl;
    if(isnan(f) || f == std::numeric_limits<double>::infinity())
	return LARGE;
    for(unsigned int i=0;i<d->tncgrads.size();i++){
	g[i] = d->tncgrads[i];
    }
    return f;
}

/*
 * which one for all of these are 
 * 1 lbfgs
 * 2 tnewton
 * 3 mma 
 * (0 refers to tnc in the optimize_tnc)
 * more detail does a different scaling of the ftol
 */

int optimize_plcp_nlopt(double *init_x,pl_calc_parallel *pl_,int numiter, int whichone, bool ad,bool moredetail, double ftol,double xtol){
    //tncplcp = pl_;
    double f;
    data_obj_nlopt d;
    d.tncparams.assign(pl_->numparams,0.0);
    d.tncgrads.assign(pl_->numparams,0.0);
    double low[pl_->numparams];
    double up[pl_->numparams];
//    double init_v[pl_->numparams];

    for(int i=0;i<pl_->numparams;i++){
//	cout << init_x[i] <<  endl;
	low[i] = 0;
	up[i] = LARGE;
//	init_v.push_back(init_x[i]);
    }
//    exit(0);
    //if((whichone == 3 || whichone == 1) && ad == true){//sometimes works, sometimes doesn't
	int fcount = 0;
	for(int i=0;i<pl_->rates.size();i++){
	    if(pl_->freeparams[fcount] != -1){
//		if (pl_->lf == false)
//		    low[pl_->freeparams[fcount]] = pl_->minrate;
		low[pl_->freeparams[fcount]] = 1e-40;
		up[pl_->freeparams[fcount]] = LARGE;
	    }
	    fcount += 1;
	}
	
	for(int i=0;i<pl_->numnodes;i++){
	    if(pl_->freeparams[fcount] != -1){
		int curponi = i;
		if(pl_->min->at(curponi) != -1)
		    low[pl_->freeparams[fcount]] = pl_->min->at(curponi);
		if(pl_->max->at(curponi) != -1)
		    up[pl_->freeparams[fcount]] = pl_->max->at(curponi);
	    }
	    fcount += 1;
	}

	//}

    
    double g[pl_->numparams];
    d.pl_ = pl_;
    int rc = 0;
    
    //LN_SBPLX is the fastest of the non deriv

    nlopt_opt opt;
    if(whichone == 1){
	opt = nlopt_create(NLOPT_LD_LBFGS, pl_->numparams);
	cout << "setting NLOPT: LD_LBFGS " << endl;
    }else if(whichone == 2){ 
	opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART, pl_->numparams);
	cout << "setting NLOPT: LD_TNEWTON_PRECOND_RESTART" << endl;
    }else if(whichone == 3){
	opt = nlopt_create(NLOPT_LD_MMA, pl_->numparams);
	cout << "setting NLOPT: LD_MMA" << endl;
    }else if(whichone == 4){
	opt = nlopt_create(NLOPT_LD_VAR2, pl_->numparams);
	cout << "setting NLOPT: LD_VAR2 " << endl;
    }else if(whichone == 5){
	if(pl_->numparams < 10000){
	    opt = nlopt_create(NLOPT_LN_SBPLX, pl_->numparams);
	    cout << "setting NLOPT : LN_SBPLX " << endl;
	}else{
	    whichone = 3;
	    opt = nlopt_create(NLOPT_LD_MMA, pl_->numparams);
	    cout << "setting NLOPT : NLOPT_LD_MMA " << endl;
	    cout << "you may want to set the plsimaniter > 100000 " << endl;
	}
    }
//    opt.set_vector_storage(10);
    
    nlopt_set_lower_bounds(opt,low);
    nlopt_set_upper_bounds(opt,up);
    nlopt_set_maxeval(opt,numiter);
    if (moredetail){
	nlopt_set_ftol_rel(opt,ftol*1e-5);
    }else{
	nlopt_set_ftol_rel(opt,ftol);
    }
    nlopt_set_xtol_rel(opt,xtol);
    if(whichone == 5)
	nlopt_set_min_objective(opt,function_plcp_nlopt_nograd, (&d));
    else if(ad==false)
	nlopt_set_min_objective(opt,function_plcp_nlopt, (&d));
    else
	nlopt_set_min_objective(opt,function_plcp_nlopt_ad, (&d));
    try{
	nlopt_result result = nlopt_optimize(opt,init_x, &f); //was init_v
//	for(int i=0;i<pl_->numparams;i++){
//	    init_x[i] = init_v[i];
//	}
	cout << "result: " << result << endl;
	nlopt_destroy(opt);
	return result;
    }catch(...){
	cout << "failed line search / last value: " << f << endl;
//	for(int i=0;i<pl_->numparams;i++){
//	    init_x[i] = init_v[i];
//	}
	nlopt_destroy(opt);
	return -1;
    }
    return rc;
}


int optimize_plcp_nlopt_ad_parallel(double *init_x,pl_calc_parallel *pl_,int numiter, int whichone, bool moredetail, double ftol, double xtol){
    //tncplcp = pl_;
    double f;
    data_obj_nlopt d;
    d.tncparams.assign(pl_->numparams,0.0);
    d.tncgrads.assign(pl_->numparams,0.0);
    double low[pl_->numparams];
    double up[pl_->numparams];
//    vector<double> init_v;

    for(int i=0;i<pl_->numparams;i++){
	low[i] = 0;
	up[i] = LARGE;
//	init_v.push_back(init_x[i]);
    }

    int fcount = 0;
    for(int i=0;i<pl_->rates.size();i++){
	if (pl_->freeparams[fcount] != -1){
//	    if (pl_->lf == false)
//		low[pl_->freeparams[fcount]] = pl_->minrate;
	    up[pl_->freeparams[fcount]] = LARGE;
	    low[pl_->freeparams[fcount]] = 1e-40;
	}
	fcount += 1;
    }

    for(int i=0;i<pl_->numnodes;i++){
	if(pl_->freeparams[fcount] != -1){
	    int curponi = i;
	    if(pl_->min->at(curponi) != -1)
		low[pl_->freeparams[fcount]] = pl_->min->at(curponi);
	    if(pl_->max->at(curponi) != -1)
		up[pl_->freeparams[fcount]] = pl_->max->at(curponi);
	}
	fcount += 1;
    }

    for(int i=0;i<pl_->numparams;i++){
	if(low[i] > init_x[i])
	    init_x[i] = low[i]+.1;
	if(up[i] < init_x[i])
	    init_x[i] = up[i];
//	cout << "- " << init_x[i] << " " << low[i] << " " << up[i] << endl;
    }

    double g[pl_->numparams];
    d.pl_ = pl_;
    int rc = 0;
    nlopt_opt opt;
    if(whichone == 1){
	opt = nlopt_create(NLOPT_LD_LBFGS, pl_->numparams); 
	cout << "setting NLOPT parallel : LD_LBFGS " << endl;
    }else if(whichone == 2){
	opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART, pl_->numparams);
	cout << "setting NLOPT parallel : LD_TNEWTON_PRECOND_RESTART" << endl;
    }else if(whichone == 3){
	opt = nlopt_create(NLOPT_LD_MMA, pl_->numparams);
	cout << "setting NLOPT parallel : LD_MMA" << endl;
    }else if(whichone == 4){
	opt = nlopt_create(NLOPT_LD_VAR2, pl_->numparams);
	cout << "setting NLOPT parallel : LD_VAR2 " << endl;
    }else if(whichone == 5){
	if(pl_->numparams < 10000){
	    opt = nlopt_create(NLOPT_LN_SBPLX, pl_->numparams);
	    cout << "setting NLOPT parallel : LN_SBPLX " << endl;
	}else{
	    whichone = 3;
	    opt = nlopt_create(NLOPT_LD_MMA, pl_->numparams);
	    cout << "setting NLOPT parallel : LD_MMA " << endl;
	    cout << "you may want to set the plsimaniter > 100000 " << endl;
	}
	
    }
    nlopt_set_lower_bounds(opt,low);
    nlopt_set_upper_bounds(opt,up);
    nlopt_set_maxeval(opt,numiter);
    if(moredetail)
	nlopt_set_ftol_rel(opt,ftol*1e-5);
    else
	nlopt_set_ftol_rel(opt,ftol);
    if(whichone == 5){
	nlopt_set_min_objective(opt,function_plcp_nlopt_nograd, (&d));
    }else{
#if HAVE_LIBADOLC
	nlopt_set_min_objective(opt,function_plcp_nlopt_ad_parallel, (&d));
#else
	nlopt_set_min_objective(opt,function_plcp_nlopt, (&d));
#endif
    }
    nlopt_set_xtol_rel(opt,xtol);
    try{
	nlopt_result result = nlopt_optimize(opt,init_x, &f);
//	for(int i=0;i<pl_->numparams;i++){
//	    init_x[i] = init_v[i];
//	}
	cout << "result: " << result << endl;
	nlopt_destroy(opt);
	return result;
    }catch(...){
	cout << "failed line search / last value: " << f << endl;
//	for(int i=0;i<pl_->numparams;i++){
//	    init_x[i] = init_v[i];
//	}
	nlopt_destroy(opt);
	return -1;
    }
    return rc;
}

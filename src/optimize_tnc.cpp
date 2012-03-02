/*
 * optimize_tnc.cpp
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "optimize_tnc.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include "tnc.h"
#include "pl_calc_parallel.h"
#include <fstream>

using namespace std;

typedef struct{
    vector <double> tncparams;
    vector <double> tncgrads;
    pl_calc_parallel * pl_;
} data_obj;

#define VERBOSE false

#define LARGE 1e+15
#define SMALL 1e-20

static int function_plcp_ad(double x[], double *f, double g[], void *state){
    data_obj *d = (data_obj *) state;
    if (isnan(x[0])){
	return 1;
    }
    for(unsigned int i=0;i<d->tncparams.size();i++){
	d->tncparams[i] = x[i];
    }
    *f = d->pl_->calc_function_gradient(&d->tncparams,&d->tncgrads);
//    cout << "*f: " << (*f) << endl;
    if (isnan(*f))
	return 1;

    int tv = d->tncgrads.size()-1;
    for(unsigned int i=0;i<d->tncgrads.size();i++){
	g[i] = d->tncgrads[i];
	if (g[i] > 1e+40)
	    g[i] = 1e+40;
	if(g[i] < 1e-40)
	    g[i] = 1e-40;
    }
    return 0;
}

static int function_plcp_ad_parallel(double x[], double *f, double g[], void *state){
    data_obj *d = (data_obj *) state;
    for(unsigned int i=0;i<d->tncparams.size();i++){
	d->tncparams[i] = x[i];
    }
    *f = d->pl_->calc_pl_function_gradient_adolc(&d->tncparams,&d->tncgrads);
//    cout << "(*f)  " << (*f) << endl;
    int tv = d->tncgrads.size()-1;
    for(unsigned int i=0;i<d->tncgrads.size();i++){
	g[i] = d->tncgrads[i];
	if (g[i] > 1e+40)
	    g[i] = 1e+40;
	if(g[i] < 1e-40)
	    g[i] = 1e-40;
    }
    return 0;
}

static int function_plcp(double x[], double *f, double g[], void *state){
    data_obj *d = (data_obj *) state;

    for(unsigned int i=0;i<d->tncparams.size();i++){
//	cout << x[i] << " ";
	d->tncparams[i] = x[i];
    }
//    cout << endl;
    *f = d->pl_->calc_pl(d->tncparams);
    //   if (isnan(*f)){
//	cout << "NAN result" << endl;
//	return 1;
    //  }
    if (isinf(*f)){
	cout << "INF result" << endl;
	return 1;
      }
//    cout << (*f) << endl;
    d->pl_->calc_gradient(d->tncparams,&d->tncgrads);
//    cout <<"-";
    for(unsigned int i=0;i<d->tncgrads.size();i++){
	g[i] = d->tncgrads[i];
	if (g[i] > 1e+40)
	    g[i] = 1e+40;
	if(g[i] < 1e-40)
	    g[i] = 1e-40;
//	cout << g[i] << " ";
    }
//    cout << endl;
    return 0;
}

int optimize_plcp_tnc(double *init_x,pl_calc_parallel *pl_,int numiter, bool ad,double inftol){
    //tncplcp = pl_;
    double f;
    data_obj d;
    d.tncparams.assign(pl_->numparams,0.0);
    d.tncgrads.assign(pl_->numparams,0.0);
    double low[pl_->numparams];
    double up[pl_->numparams];
    
    for(int i=0;i<pl_->numparams;i++){
	low[i] = 0;
	up[i] = LARGE;
    }
    int fcount = 0;
    for(int i=0;i<pl_->rates.size();i++){
	if(pl_->freeparams[fcount] !=-1){
//	    if(pl_->lf == false)
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
	if (low[i]==0)
	    low[i] = 1e-40;
//	cout << low[i] << " " << up[i] << " " << init_x[i] << endl;
    }

    for(int i=0;i<pl_->numparams;i++){
	if(low[i] > init_x[i])
	    init_x[i] = low[i];
	if(up[i] < init_x[i])
	    init_x[i] = up[i];
//	cout << "- " << init_x[i] << " " << low[i] << " " << up[i] << endl;
    }

    double g[pl_->numparams];
    d.pl_=pl_;
//maxCGit = -1
    int rc, maxCGit = -1, maxnfeval = 10000, nfeval;
    double eta = -1, stepmx = -1, accuracy = -1, fmin = -1, ftol = inftol, rescale = -1;
    if (ad == false){
	rc = tnc(pl_->numparams, init_x, &f, g, function_plcp, (&d), low, up, NULL, TNC_MSG_NONE,
		 maxCGit, maxnfeval, eta, stepmx, accuracy, fmin, ftol, rescale, &nfeval);
    }else{
	rc = tnc(pl_->numparams, init_x, &f, g, function_plcp_ad,(&d), low,up, NULL, TNC_MSG_NONE,
	     maxCGit, maxnfeval, eta, stepmx, accuracy, fmin, ftol, rescale, &nfeval);
    }

    cout << "nfeval: " << nfeval << " rc: " << tnc_rc_string[rc-TNC_MINRC] << endl;
    return rc;
}

int optimize_plcp_tnc_ad_parallel(double *init_x,pl_calc_parallel *pl_,double inftol){
    //tncplcp = pl_;
    double f;
    data_obj d;
    d.tncparams.assign(pl_->numparams,0.0);
    d.tncgrads.assign(pl_->numparams,0.0);
    double low[pl_->numparams];
    double up[pl_->numparams];
    
    for(int i=0;i<pl_->numparams;i++){
	low[i] = 0;
	up[i] = LARGE;
    }

    int fcount = 0;
    for(int i=0;i<pl_->rates.size();i++){
	if(pl_->freeparams[fcount] != -1){
//	    if(pl_->lf == false)
//		low[pl_->freeparams[fcount]] = pl_->minrate;
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
    for(int i=0;i<pl_->numparams;i++){
	if(low[i] == 0)
	    low[i] = 1e-40;
    }

    for(int i=0;i<pl_->numparams;i++){
	if(low[i] > init_x[i])
	    init_x[i] = low[i];
	if(up[i] < init_x[i])
	    init_x[i] = up[i];
//	cout << "- " << init_x[i] << " " << low[i] << " " << up[i] << endl;
    }

    double g[pl_->numparams];
    d.pl_=pl_;
//maxCGit = -1
    int rc, maxCGit = -1, maxnfeval = 10000, nfeval;
    double eta = -1, stepmx = -1, accuracy = -1, fmin = -1, ftol = inftol, rescale = -1;
#if HAVE_LIBADOLC
	rc = tnc(pl_->numparams, init_x, &f, g, function_plcp_ad_parallel, (&d), low, up, NULL, TNC_MSG_NONE,maxCGit, maxnfeval, eta, stepmx, accuracy, fmin, ftol, rescale, &nfeval);
#else
	rc = tnc(pl_->numparams, init_x, &f, g, function_plcp, (&d), low, up, NULL, TNC_MSG_NONE,maxCGit, maxnfeval, eta, stepmx, accuracy, fmin, ftol, rescale, &nfeval);
#endif

	     
    cout << "nfeval: " << nfeval << " rc: " << tnc_rc_string[rc-TNC_MINRC] << endl;
    return rc;
}


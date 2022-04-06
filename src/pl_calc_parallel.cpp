/*
 * pl_calc.cpp
 *
 *  Created on: Feb 2, 2010
 *      Author: smitty
 */

#include "pl_calc_parallel.h"
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <limits>
#include "myrad.h" //for auto diff
//#include <adolc/adolc.h>
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>

using namespace std;

#define SMALL 1e-20
#define LARGE 1e+15

pl_calc_parallel::pl_calc_parallel(){}

 
void pl_calc_parallel::setup_starting_bits(const vector<int> * parents_nds, const vector<int> * child_cnts, 
					    const vector<int> * vfree, const vector<double> * char_dur, 
					    const vector<double> * log_fact_ch, const vector<double> * vmn, 
					    const vector<double> * vmx, const vector<double> * start_dates, 
					    const vector<double> * start_rates, const vector<double> * start_durations){
    log_pen = false;
    parents_nds_ints = parents_nds;
    free = vfree;
    char_durations = char_dur;
    log_fact_char_durations = log_fact_ch;
    min = vmn;
    max = vmx;
    child_counts = child_cnts;
    numnodes = parents_nds_ints->size();
    dates = vector<double>();
    rates = vector<double>();
    durations = vector<double>();
    cvnodes = vector<int>();
    smoothing = 1;
    for(int i=0; i<numnodes; i++){
        dates.push_back(start_dates->at(i));
        rates.push_back(start_rates->at(i));
        durations.push_back(start_durations->at(i));
        cvnodes.push_back(0);
    }
    //FOR AUTODIFF
    advdurations_adc = new adouble[durations.size()];
    advdates_adc = new adouble[dates.size()];
    advrates_adc= new adouble[rates.size()];
    adc_size = 1;
    adx = new adouble[adc_size];
    paramverbose = false;
    minrate = 1e-10;
    ftol = 1e-7;
    isConstrained = true;
    penalty_boundary = 0.0025;
    //cout << "set everything up" << endl;
}

double pl_calc_parallel::calc_pl(vector<double> & params){
    for (unsigned int i=0;i<params.size();i++){
	if (params[i] < 0){
	    return LARGE;
	}
    }
    int pcount = 0;
    int fcount = 0;
    //rates
    for(int i=0;i<rates.size();i++){
	if(freeparams[fcount] != -1){//TODO: check this cvnodes
	    rates[i] = params[freeparams[fcount]];//params[pcount];
//	    cout << "r: " << rates[i]/1000. << " " << freeparams[fcount]<<  endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //dates

    for(int i=0;i<dates.size();i++){
	if(freeparams[fcount] != -1){//is free, negative is not free
	    dates[i] = params[freeparams[fcount]];//params[pcount];
//	    cout << "d: " << dates[i]/1000. << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    double ret = set_durations();
    if (ret == LARGE){
	return LARGE;
    }
    double ll = calc_log_like();
    //cout << ll << endl;
    if(paramverbose){
	ofstream paramverboseFile;
	if(lf){
	    paramverboseFile.open("paramverboself",ios::app);
	}else{
	    paramverboseFile.open("paramverbosepl",ios::app);
	}
	paramverboseFile << ll;
	for (unsigned int i=0;i<params.size();i++){
	    paramverboseFile << "\t" << params[i];
	}
	paramverboseFile << "\n";
	paramverboseFile.close();
    }
    if (!lf){
	double tp = calc_penalty();
	double rp = calc_roughness_penalty();
	double pl = (ll+tp)+smoothing*rp;//(smoothing*rp);
//	cout << pl << " " << (ll+tp) << " " << tp << " " << rp  << " "<<  smoothing <<  " " << smoothing*rp << endl;
	return pl;
    }else{
	return ll; // may need to seperate this for autodiff lf
    }
    return 0;
}

double pl_calc_parallel::calc_function_gradient(vector<double> * params, vector<double> * g){
    if (lf == true)
	return calc_lf_function_gradient(params,g);
    else
	return calc_pl_function_gradient(params,g);
    return 0;
}

//with rad autodiff
//rad doesn't work with openmp
double pl_calc_parallel::calc_lf_function_gradient(vector<double> * params, vector<double> * g){
    ADvar rv, v[params->size()];
    //initialize independent active variables
    for(int i=0;i<params->size();i++)
	v[i] = params->at(i);
    //calculate the pl
    rv = 0;
    int pcount = 0;
    int fcount = 0;
    //rates
    ADvar advrates[rates.size()];
    for(int i=0;i<rates.size();i++){
	if(freeparams[fcount] != -1){
	    // if(v[freeparams[fcount]] <= minrate)
	    //	advrates[i] = minrate;
	    //else
		advrates[i] = v[freeparams[fcount]];//params[pcount];
	    // cout << "rates " << i << " " << rates[i] << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //dates
    ADvar advdates[dates.size()];
    for(int i=0;i<dates.size();i++){
	//cout << dates[i] << endl;
	advdates[i] = dates[i];
	if(freeparams[fcount] != -1){//is free, negative is not free
	    advdates[i] = v[freeparams[fcount]];//params[pcount];
	    // cout << "dates " << i << " " << dates[i] << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //set durations
    ADvar advdurations[durations.size()];
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if (min->at(curponi) != -1){
	    if(advdates[curponi].val() < min->at(curponi)){
		//cout << "min " << curponi << " " << min->at(curponi) << " " << advdates[curponi].val() << endl;
		return LARGE;
	    }
	}
	if (max->at(curponi) != -1){
	    if(advdates[curponi].val() > max->at(curponi)){
		//	cout << "max " << curponi << " " << max->at(curponi) << " " << advdates[curponi].val() << endl;
		return LARGE;
	    }
	}
	if (i != 0){//not equal to the root
	    advdurations[curponi] = (advdates[parents_nds_ints->at(curponi)]) - (advdates[curponi]);//add fabs
	    if(advdurations[curponi].val() < 0)
		return LARGE;
	}
    }
    
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if(cvnodes[curponi] == 0){
	    if (curponi != 0){//not the root
		double c = char_durations->at(curponi);
		double lf = log_fact_char_durations->at(curponi);
		//cout << rt << " " << d  << " " << c  << " " << x << " " << lf << " " << l << endl;
		//rv += (-(advrates[curponi]*advdurations[curponi])+(c*log(advrates[curponi]*advdurations[curponi])-lf));//was int(c)
		rv += (-(c*log(advrates[curponi]*advdurations[curponi])-(advrates[curponi]*advdurations[curponi])-lf));
	    }
	}
    }

    //compute the gradient
    ADcontext::Gradcomp();

    for(int i=0;i<params->size();i++)
	g->at(i) = v[i].adj(); // need to see the direction of these

    //return the results
    return rv.val();
}

//with autodiff_adolc
//adolc requires a library but allows openmp
int pl_calc_parallel::calc_lf_function_gradient_adolc_void(int size, double * xp, double * yp, int tape){
    adouble y = 0;
    trace_on(tape,1);
    //initialize independent active variables
    for(int i=0;i<size;i++){
	adx[i] <<= xp[i];
    }
    //calculate the pl
    int ret = 1;
    int pcount = 0;
    int fcount = 0;
    //rates
    for(int i=0;i<rates.size();i++){
	if(freeparams[fcount] != -1){
	    // if(adx[freeparams[fcount]] <= minrate)
	    //	advrates_adc[i] = minrate;
	    //else
		advrates_adc[i] = adx[freeparams[fcount]];//params[pcount];
	    // cout << "rates " << i << " " << rates[i] << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //dates
    for(int i=0;i<dates.size();i++){
	advdates_adc[i] = dates[i];
	if(freeparams[fcount] != -1){//is free, negative is not free
	    advdates_adc[i] = adx[freeparams[fcount]];//params[pcount];
	    // cout << "dates " << i << " " << dates[i] << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //set durations
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if (min->at(curponi) != -1){
	    if(advdates_adc[curponi].value() < min->at(curponi)){
//		cout << "min " << curponi << " " << min->at(curponi) << " " << advdates_adc[curponi].value() << endl;
		ret = 0;//NEED TO FIX THIS
	    }
	}
	if (max->at(curponi) != -1){
	    if(advdates_adc[curponi].value() > max->at(curponi)){
//		cout << "max " << curponi << " " << max->at(curponi) << " " << advdates_adc[curponi].value() << endl;
		ret = 0;//NEED TO FIX THIS
	    }
	}
	if (i != 0){//not equal to the root
	    advdurations_adc[curponi] = (advdates_adc[parents_nds_ints->at(curponi)]) - (advdates_adc[curponi]);//add fabs
	    if(advdurations_adc[curponi].value() < 0){
//		cout << "duration " << curponi << " " << advdurations_adc[curponi].value() << " " << advdates_adc[parents_nds_ints->at(curponi)].value() << " " << advdates_adc[curponi].value() << endl;
		ret = 0;//NEED TO FIX THIS
	    }
	}
    }
    if (ret == 0){//TODO: MAKE SURE THIS WORKS, memory 
	trace_off();
	return ret;
    }
    
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if(cvnodes[curponi] == 0){
	    if (curponi != 0){//not the root
		double c = char_durations->at(curponi);
		double lf = log_fact_char_durations->at(curponi);
		//cout << rt << " " << d  << " " << c  << " " << x << " " << lf << " " << l << endl;
		//y += (-(advrates_adc[curponi]*advdurations_adc[curponi])+(c*log(advrates_adc[curponi]*advdurations_adc[curponi])-lf));//was int(c)
		y += (-(c*log(advrates_adc[curponi]*advdurations_adc[curponi]) - (advrates_adc[curponi]*advdurations_adc[curponi]) - lf));
	    }
	}
    }
    y >>= *yp;
    trace_off();
    return ret;
}

//NEED TO DELETE THE DURATION, DATES, RATES
double pl_calc_parallel::calc_lf_function_gradient_adolc(vector<double> * params, vector<double> * g){
    double * xp = new double[params->size()];double yp = 0.0;
    short int tag = 0;
    for(int i=0;i<params->size();i++){
	xp[i] = params->at(i);
    }
    int ps = params->size();
    
    if (ps != adc_size){
	adc_size = ps;
	delete []adx;
	adx = new adouble[ps];
    }

    double ret = calc_lf_function_gradient_adolc_void(ps,xp,&yp,tag);
    if (ret == 0){
	delete []xp;
	return LARGE;
    }
    double * gd = new double[params->size()];
    //cout << yp << " " << ret << endl;
    int b = gradient(tag,ps,xp,gd);
    //cout << b << endl;
    for(int i=0;i<params->size();i++)
	g->at(i) = gd[i]; // need to double check the direction
    
    delete []gd;
    delete []xp;
    //return the results
    return yp;
}

//with autodiff

double pl_calc_parallel::calc_pl_function_gradient(vector<double> * params, vector<double> * g){
    ADvar rv, v[params->size()];

    //initialize independent active variables
    for(int i=0;i<params->size();i++)
	v[i] = params->at(i);
    //calculate the pl
    rv = 0;
    int pcount = 0;
    int fcount = 0;
    //rates
    ADvar advrates[rates.size()];
    for(int i=0;i<rates.size();i++){
	if(freeparams[fcount] != -1){
	    // if (v[freeparams[fcount]] <= minrate)
	    //	advrates[i] = minrate;
	    //else
		advrates[i] = v[freeparams[fcount]];//params[pcount];
	    // cout << "rates " << i << " " << rates[i] << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //dates
    ADvar advdates[dates.size()];
    for(int i=0;i<dates.size();i++){
	advdates[i] = dates[i];
	if(freeparams[fcount] != -1){//is free, negative is not free
	    advdates[i] = v[freeparams[fcount]];//params[pcount];
	    //cout << "dates " << i << " " << dates[i] << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //set durations
    ADvar advdurations[durations.size()];
    ADvar tpen = 0.;
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if (min->at(curponi) != -1){
	    if(advdates[curponi].val() < min->at(curponi)){
		//cout << "min " << curponi << " " << min->at(curponi) << " " << advdates[curponi].val() << endl;
		return LARGE;
	    }
	    if((*pen_min)[curponi] != -1){
		if ((advdates[curponi]-(*pen_min)[curponi]) > 0)
		    tpen += 1./(advdates[curponi]-(*pen_min)[curponi]);
	    } 
	}
	if (max->at(curponi) != -1){
	    if(advdates[curponi].val() > max->at(curponi)){
		//cout << "max" << endl;
		return LARGE;
	    }
	    if((*pen_max)[curponi] != -1){
		if(((*pen_max)[curponi]-advdates[curponi]) > 0)
		    tpen += 1./((*pen_max)[curponi]-advdates[curponi]);
	    } 
	}
	if (i != 0){//not equal to the root
	    advdurations[curponi] = (advdates[parents_nds_ints->at(curponi)]) - (advdates[curponi]);//add fabs
	    if(advdurations[curponi].val() < 0){
		return LARGE;
	    }else if(advdurations[curponi].val() == 0){
		advdurations[curponi] = numeric_limits<double>::min();
	    }
	}
    }
    if(isinf(tpen.val()) || isnan(tpen.val()))
	tpen = 0;
    rv += (penalty_boundary)*tpen;//need to check on this penalty boundary
    //calc_ll
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if(cvnodes[curponi] == 0){
	    if (curponi != 0){//not the root
		double c = char_durations->at(curponi);
		double lf = log_fact_char_durations->at(curponi);
		double x = advrates[curponi].val()*advdurations[curponi].val();
		if(x > 0){
		    rv += (-(c*log(advrates[curponi]*advdurations[curponi])-(advrates[curponi]*advdurations[curponi])-lf));
		}else{
		    cout << c << " " << lf << " " <<advrates[curponi].val() << " " << advdurations[curponi].val() << endl;
		    if (c > 0.0)
			rv += LARGE;
		    else if (c == 0)
			rv += 0;
		}
		//rv += (-(advrates[curponi]*advdurations[curponi])+(c*log(advrates[curponi]*advdurations[curponi])-lf));//was int(c)

	    }
	}
    }
    //calc_roughness
    ADvar rp;
    ADvar ss;
    ADvar s;
    ADvar tomy;
    ADvar r;
    rp = 0.0;
    ss = 0.0;
    s = 0.0;
    tomy = 0;
    for(int i=0;i<numnodes;i++){
	int curponi = i;//start here
	if(cvnodes[curponi]==0){
	    if(curponi != 0){// isn't the root
		if(parents_nds_ints->at(curponi) != 0){
		    rp += (advrates[curponi] - advrates[parents_nds_ints->at(curponi)])*(advrates[curponi] - advrates[parents_nds_ints->at(curponi)]);
		}
	    }else{
		for(int j=0;j<children_vec->at(curponi).size();j++){
		    int curch = children_vec->at(curponi)[j];
		    if(log_pen){//make sure that this is efficient
			r = log(advrates[curch]);
			s += r;
			ss += r*r;
		    }else{
			s += advrates[curch];
			ss += advrates[curch]*advrates[curch]; 
		    }
		    tomy+=1;
		}
	    }
	}
    }
    rp += (ss-s*s/tomy)/tomy;
//    rp += 2*( ss-s*s/tomy);
    rv += (smoothing*rp);
    //rv -= (smoothing*rp);
    //compute the gradient
    ADcontext::Gradcomp();

    for(int i=0;i<params->size();i++)
	g->at(i) = v[i].adj(); // need to double check the direction

    //return the results
    return rv.val();
}

//with autodiff_adolc
int pl_calc_parallel::calc_pl_function_gradient_adolc_void(int size, double * xp, double * yp, int tape){
    adouble y = 0;
    trace_on(tape,1);
    //initialize independent active variables
    for(int i=0;i<size;i++){
	adx[i] <<= xp[i];
    }
    //calculate the pl
    int ret = 1;
    int pcount = 0;
    int fcount = 0;
    //rates
    int rs = rates.size();
    for(int i=0;i<rs;i++){
	if(freeparams[fcount] != -1){
	    //  if(adx[freeparams[fcount]] <= minrate){
	    //	advrates_adc[i] = minrate;
	    //}else
		advrates_adc[i] = adx[freeparams[fcount]];//params[pcount];
	    // cout << "rates " << i << " " << rates[i] << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //dates
    int ds = dates.size();
    for(int i=0;i<ds;i++){
	advdates_adc[i] = dates[i];
	if(freeparams[fcount] != -1){//is free, negative is not free
	    advdates_adc[i] = adx[freeparams[fcount]];//params[pcount];
	    // cout << "dates " << i << " " << dates[i] << endl;
	    pcount += 1;
	}
	fcount += 1;
    }
    //set durations
    adouble tpen = 0.;
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if (min->at(curponi) != -1){
	    if(advdates_adc[curponi].value() < min->at(curponi)){
//		cout << "min " << curponi << " " << min->at(curponi) << " " << advdates_adc[curponi].val() << endl;
		ret = 0;//NEED TO FIX THIS
	    }
	    if((*pen_min)[curponi] != -1){
		adouble tv = (*pen_min)[curponi];
		if((advdates_adc[curponi]-tv) > 0){
		    tpen += 1./(advdates_adc[curponi]-tv);
		}
	    } 
	}
	if (max->at(curponi) != -1){
	    if(advdates_adc[curponi].value() > max->at(curponi)){
//		cout << "max" << endl;
		ret = 0;//NEED TO FIX THIS
	    }
	    if((*pen_max)[curponi] != -1){
		adouble tv = (*pen_max)[curponi];
		if((tv-advdates_adc[curponi])>0){
		    tpen += 1./(tv-advdates_adc[curponi]);
		}
	    } 
	}
	if (i != 0){//not equal to the root
	    advdurations_adc[curponi] = (advdates_adc[parents_nds_ints->at(curponi)]) - (advdates_adc[curponi]);//add fabs
	    if(advdurations_adc[curponi].value() < 0){
//		cout << "duration" << endl;
		ret = 0;//NEED TO FIX THIS
	    }else if(advdurations_adc[curponi].value() == 0){
		advdurations_adc[curponi] = numeric_limits<double>::min();
	    }
	}
    }
    if(isinf(tpen.value()) || isnan(tpen.value()))
	tpen = 0;
    y += (penalty_boundary)*tpen;
    if (ret == 0){//TODO: MAKE SURE THIS WORKS, memory 
	trace_off();
	return ret;
    }
    
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if(cvnodes[curponi] == 0){
	    if (curponi != 0){//not the root
		double c = char_durations->at(curponi);
		double lf = log_fact_char_durations->at(curponi);
		double x = advrates_adc[curponi].value()*advdurations_adc[curponi].value();
		if(x > 0){
		    y += (-(c*log(advrates_adc[curponi]*advdurations_adc[curponi]) - (advrates_adc[curponi]*advdurations_adc[curponi]) - lf));
		}else{
		    cout << c << " " << lf << " " <<advrates_adc[curponi].value() << " " << advdurations_adc[curponi].value() << endl;
		    if (c > 0.0)
			y += LARGE;
		    else if (c == 0)
			y += 0;
		}

		//cout << rt << " " << d  << " " << c  << " " << x << " " << lf << " " << l << endl;
		//y += (-(advrates_adc[curponi]*advdurations_adc[curponi])+(c*log(advrates_adc[curponi]*advdurations_adc[curponi])-lf));//was int(c)
	    }
	}
    }
    //calc_roughness
    adouble rp;
    adouble ss;
    adouble s;
    adouble tomy;
    adouble r;
    rp = 0.0;
    ss = 0.0;
    s = 0.0;
    tomy = 0;
    for(int i=0;i<numnodes;i++){
	int curponi = i;//start here
	if(cvnodes[curponi]==0){//TODO: cv
	    if(curponi != 0){// isn't the root
		if(parents_nds_ints->at(curponi) != 0){
		    rp += (advrates_adc[curponi] - advrates_adc[parents_nds_ints->at(curponi)])*(advrates_adc[curponi] - advrates_adc[parents_nds_ints->at(curponi)]);
		}
	    }else{
		for(int j=0;j<children_vec->at(curponi).size();j++){
		    int curch = children_vec->at(curponi)[j];
		    if(log_pen){
			r = log(advrates_adc[curch]);
			ss += r*r;
			s += r;
		    }else{
			ss += advrates_adc[curch] * advrates_adc[curch];
			s += advrates_adc[curch];
		    }
		    tomy++;
		}
	    }
	}
    }
    rp += (ss-s*s/tomy)/tomy;
    y += (smoothing*rp);

    y >>= *yp;
    trace_off();
    return ret;
}

double pl_calc_parallel::calc_pl_function_gradient_adolc(vector<double> * params, vector<double> * g){
    double * xp = new double[params->size()];double yp = 0.0;
    //TODO: add for openmp
    short int tag = omp_get_thread_num(); 

    //cout << "TAG: "<< tag<< endl; 
    int ps = params->size();
    for(int i=0;i<ps;i++){
	xp[i] = params->at(i);
    }
    
    if (ps != adc_size){
	adc_size = ps;
	delete []adx;
	adx = new adouble[ps];
    }

    double ret = calc_pl_function_gradient_adolc_void(ps,xp,&yp,tag);
    if (ret == 0){
	delete []xp;
	return LARGE;
    }
    double * gd = new double[params->size()];
    //cout << yp << " " << ret << endl;
    int b = gradient(tag,ps,xp,gd);
    //cout << b << endl;
    for(int i=0;i<params->size();i++)
	g->at(i) = gd[i]; // need to double check the direction
    
    delete []gd;
    delete []xp;
    //return the results
    return yp;
}

void pl_calc_parallel::calc_gradient(vector<double> & params, vector<double> * g){
    if (lf == true)
	calc_lf_gradient(params,g);
    else
	calc_pl_gradient(params, g);
}

/*
 * No CV in the calc_lf_gradient because this is typically 
 * going to be called before any of that is done
 */
void pl_calc_parallel::calc_lf_gradient(vector<double> & params, vector<double> * g){
    //rate derivative
    int count = 0;
    double rate = params[0];
    double sumbl = 0; double sumtimes = 0;
    for(unsigned int i=0;i<numnodes;i++){
	int curponi = i;
	//sumbl += int(char_durations->at(curponi));
	sumbl += char_durations->at(curponi);
	sumtimes += durations[curponi];
    }
    g->at(count) =- (sumbl/rate-sumtimes);

    if (isnan(g->at(count))){
	//cout << sumbl << " " << rate << " " << sumtimes << endl;
	//cout << "isnan" << endl;
	//exit(0);
	for (int i=0;i < params.size(); i++)
	    g->at(i) = 1000000;
	return;
    } 
    count += 1;
    
    //time derivative
    for(unsigned int i=0;i<numnodes;i++){
	int curponi = i;
	if (free->at(curponi) == 1){
	    g->at(count) = 0;
	    if(curponi != 0){//not the root
		if(char_durations->at(curponi) == 0.0){  
		    g->at(count) = rate;
		}else{
		    g->at(count) = -char_durations->at(curponi)/durations[curponi]+rate;
		}
	    }
	    for(int j=0;j<children_vec->at(curponi).size();j++){//WAY FASTER THAN ABOVE
		int curch = children_vec->at(curponi)[j];
		if(char_durations->at(curch) == 0.0){
		    g->at(count) -= rate;
		}else{
		    g->at(count) += char_durations->at(curch)/durations[curch]-rate;
		}
	    }
	    
	    g->at(count) *= -1;
	    count += 1;
	}
    }
}

void pl_calc_parallel::calc_pl_gradient(vector<double> & params,vector<double> * g){
    int icount = numnodes-1;
    int subcv = calc_isum(cvnodes);
    icount -= subcv;
    int pos = numnodes;
//time derivative
//#pragma omp parallel for default(none) schedule(static,100) shared(icount,g) num_threads(4)
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if (free->at(curponi) == 1){// greater than 0 is free, negative is not free
	    g->at(icount) = 0.0;
	    if(curponi != 0){//not the root
		if(char_durations->at(curponi) == 0.0){
		    g->at(icount) = rates[curponi];
		}else{
		    //g->at(icount) = -nd->char_duration/nd->duration+nd->rate;
		    g->at(icount) = -char_durations->at(curponi)/durations[curponi]+rates[curponi];
		}
	    }
	    for(int j=0;j<children_vec->at(curponi).size();j++){
		int curch = children_vec->at(curponi)[j];
		if(cvnodes[curch] == 0){
		    if(char_durations->at(curch) == 0.0){
			g->at(icount) -= rates[curch];
		    }else{
			g->at(icount) += char_durations->at(curch)/durations[curch]-rates[curch];
		    }
		}
	    }

	    g->at(icount) *= -1;//-1;
	    icount += 1;
	}
    }
//rate derivative//this one is additive
    icount = 0;
//#pragma omp parallel for default(none) schedule(static,100) shared(icount,g) num_threads(4)
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if(cvnodes[curponi] == 0){
	    int tomy = 0;
	    double meanr = 0.0;
	    double tg = 0.0;
	    double lograte;
	    if(curponi != 0){//not the root
		tg = char_durations->at(curponi)/rates[curponi] - durations[curponi];
		if(log_pen)
		    lograte = log(rates[curponi]);
		if(parents_nds_ints->at(curponi) == 0){//parent is the root
		    int curp = parents_nds_ints->at(curponi);
		    for(int j=0;j<children_vec->at(curp).size();j++){
			int curch = children_vec->at(curp)[j];
			if(cvnodes[curch] == 0){
			    ++tomy; 
			    if(log_pen)
				meanr += log(rates[curch]);
			    else
				meanr += rates[curch];
//			    printf("ch1 %f\t%f\t%i\n",rates[curch],meanr,tomy);
			}
		    }
		    meanr/=tomy;
		    if(log_pen)
			tg += -(2*smoothing/rates[curponi])*(lograte-meanr)/tomy;
		    else
			tg += -(2*smoothing)*(rates[curponi]-meanr)/tomy;
//		    cout << "f "<<  tg << endl;
		    for(int j=0;j<children_vec->at(curponi).size();j++){
			int curch = children_vec->at(curponi)[j];
			if (cvnodes[curch] == 0){
			    if(log_pen)
				tg += 2*smoothing*(log(rates[curch]) - lograte)/rates[curponi];
			    else
				tg += 2*smoothing*(rates[curch] - rates[curponi]);
			}
		    }
//		    cout << tg << endl;
		}else{
		    int curp = parents_nds_ints->at(curponi);
		    if(log_pen)
			tg += (-2*smoothing)*(lograte - log(rates[curp]))/rates[curponi];		  
		    else
			tg += (-2*smoothing)*(rates[curponi] - rates[curp]);
		    for(int j=0;j<children_vec->at(curponi).size();j++){
			int curch = children_vec->at(curponi)[j];
			if (cvnodes[curch] == 0){
			    if(log_pen)
				tg += 2*smoothing*(log(rates[curch])-lograte)/rates[curponi];
			    else
				tg += 2*smoothing*(rates[curch]-rates[curponi]);
			}
		    }
//		    cout << "notroot " << -tg << endl;
		}
		g->at(icount) = -tg;//-1;
		icount += 1;
	    }
	}
    }
//    exit(0);
}


double pl_calc_parallel::calc_log_like(){
    double ll = 0;
//marginal speed up
//#pragma omp parallel for default(none) schedule(dynamic,1000) reduction(+:ll) num_threads(2)
    for(int i=0;i<numnodes;i++){
	int curponi = i;
	if(cvnodes[curponi] == 0){
	    if (curponi != 0){//not the root
		double rt = rates[curponi];
		double d = durations[curponi];
		double c = char_durations->at(curponi);
		double x = rt*d;
		double lf = log_fact_char_durations->at(curponi);
		double l;
		if(x > 0.0)
		    l = -(c*log(x)-x-lf);
		else if(x == 0)
		    if(c > 0.0)
			l = LARGE;
		    else if(c == 0)
			l = 0.0;
//		cout << rt << " " << d  << " " << c  << " " << x << " " << lf << " " << l << endl;
		ll += l;
//		cout << ll << endl;
	    }
	}
    }
    return ll;
}

//need to add whether this will be log penalty
double pl_calc_parallel::calc_roughness_penalty(){
    double su = 0;
    double ss = 0;
    double s = 0;
    double tomy = 0;
    for(int i=0;i<numnodes;i++){
	int curponi = i;//start here
	if(cvnodes[curponi]==0){
	    if(curponi != 0){// is not the root
		if(parents_nds_ints->at(curponi) != 0){
		    double sm = rates[curponi] - rates[parents_nds_ints->at(curponi)] ;
		    su += sm*sm;
		}
	    }else{
		for(int j=0;j<children_vec->at(curponi).size();j++){
		    int curch = children_vec->at(curponi)[j];
		    double r = rates[curch];
		    if(log_pen)
			r = log(r);
		    s += r;
		    ss += r*r; 
		    ++tomy;
		}
	    }
        }
    }
    su += (ss-s*s/tomy)/tomy;
//    cout << "ss " << ss << " s " << s << " tomy " << tomy << " su "<< su << endl;
    return su;
}

double pl_calc_parallel::set_durations(){
    for(int i=0;i<numnodes;i++){
	int curponi = i;
//	if ((*min)[curponi] != -1){//just set to 0 now 
	if(dates[curponi] < (*min)[curponi]){
//		cout << "min " << curponi << " " << min->at(curponi) << " " << dates[curponi] << endl;
		return LARGE;
	    }
//	}
//	if ((*max)[curponi] != -1){//just set to 1e+20 now
	    if(dates[curponi] > (*max)[curponi]){
//		cout << "max " << curponi << " " << max->at(curponi) << " " << dates[curponi] << endl;
		return LARGE;
	    } 
//	}
	if (i != 0){//not equal to the root
	    durations[curponi] = dates[(*parents_nds_ints)[curponi]] - dates[curponi];
	    if (durations[curponi] < 0){
//		cout << "dur " << curponi << " " << dates[parents_nds_ints->at(curponi)] << " " << dates[curponi] << endl;
		return LARGE;
	    }
	}
	if(durations[curponi] == 0 || isnan(durations[curponi])){//TODO:set this locally
	        durations[curponi] = numeric_limits<double>::min();
	}
    }
    return 0.0;
}


int pl_calc_parallel::calc_isum(vector<int> & vec){
    int rsum=0; 
    for(int i=0;i<vec.size();i++){rsum+=vec[i];} 
    return rsum;
}


void pl_calc_parallel::set_freeparams(int nump, bool ilf, vector<int> * freep, vector<double> * params){
    freeparams.clear();
    lf = ilf;
    for(int i=0;i<freep->size();i++){
//	cout << freep->at(i) <<  " " ;
	freeparams.push_back(freep->at(i));
    }
//    cout << endl;
    numparams = nump;
    params->clear();
    int pcount = 0;//param count
    int fcount = 0;//free param count
    //rates
    for(int i=0;i<rates.size();i++){
	if(freeparams[fcount] != -1){//TODO: check this cvnodes bit
	    params->push_back(rates[i]);
	    pcount += 1;
	    if(lf == true){
		fcount = rates.size();
		break;
	    }
	}
	fcount += 1;
    }
//    cout << fcount << endl;
    //dates
    for(int i=0;i<dates.size();i++){
	if(freeparams[fcount] != -1){//is free, negative is not free
	    params->push_back(dates[i]);
	    pcount += 1;
	}
	fcount += 1;
    }
//    cout << numparams << " " << params->size() << " " << (numparams == params->size()) << endl;
    //assert numparams == params.size();
}

vector<int> * pl_calc_parallel::get_cv_nodes(){
    return &cvnodes;
}

void pl_calc_parallel::set_cv_nodes(vector<int> & incvnodes){
    cvnodes = incvnodes;
}

void pl_calc_parallel::delete_ad_arrays(){
    delete []adx;
    delete []advrates_adc;
    delete []advdates_adc;
    delete []advdurations_adc;
}

/*
 * should be able to add this to the duration calculator
 */
double pl_calc_parallel::calc_penalty(){
    double rk = penalty_boundary; // TODO: need to double check this
    double tpen=0;
    for(int i=0;i<numnodes;i++){
	int curponi = i;
//add if free, add if not tip	
	if((*pen_min)[curponi] != -1){
	    tpen += 1./(dates[curponi]-(*pen_min)[curponi]);
	    if(isinf(tpen))
		tpen = 0;
//	    cout << "min " << tpen << " " << dates[curponi] << endl;
	}
	
	if((*pen_max)[curponi] != -1){
	    tpen += 1./((*pen_max)[curponi]-dates[curponi]);
	    if(isinf(tpen))
		tpen = 0;
//	    cout << "max " << tpen << " " << dates[curponi] << endl;
	} 
	
//	cout << "pen: " << tpen << endl;
    }
    return rk* tpen;
}

void pl_calc_parallel::set_log_pen(bool lpen){
    log_pen = lpen;
}

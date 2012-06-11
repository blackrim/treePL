/*
 * pl_calc_parallel.h
 *
 *  Created on: Feb 2, 2010
 *      Author: smitty
 */

#ifndef PL_CALC_PARALLEL_H_
#define PL_CALC_PARALLEL_H_

#include <map>
#include <string>
#include <vector>
#include <adolc/adolc.h>
using namespace std;

class pl_calc_parallel{
private:
	double calc_log_like();
	double calc_roughness_penalty();
	double set_durations();
	double calc_penalty();
	vector<int> cvnodes;
	adouble * adx;
	adouble ady;
	adouble * advrates_adc;
	adouble * advdates_adc;
	adouble * advdurations_adc;
	int adc_size;
	bool log_pen;

public:
	pl_calc_parallel();
	double calc_pl(vector<double> & params);
	double calc_function_gradient(vector<double> * params, vector<double> * g);
	double calc_lf_function_gradient(vector<double> * params, vector<double> * g);
	int calc_lf_function_gradient_adolc_void(int size, double * xp, double * yp, int tape);
	double calc_lf_function_gradient_adolc(vector<double> * params, vector<double> * g);
	int calc_pl_function_gradient_adolc_void(int size, double * xp, double * yp, int tape);
	double calc_pl_function_gradient_adolc(vector<double> * params, vector<double> * g);
	double calc_pl_function_gradient(vector<double> * params, vector<double> * g);
	void calc_gradient(vector<double> & params, vector<double> * g);
	void calc_lf_gradient(vector<double> & params, vector<double> * g);
	void calc_pl_gradient(vector<double> & params,vector<double> * g);
	void setup_starting_bits(const vector<int> * parents_nds, const vector<int> * child_cnts, 
					    const vector<int> * vfree, const vector<double> * char_dur, 
					    const vector<double> * log_fact_ch, const vector<double> * vmn, 
					    const vector<double> * vmx, const vector<double> * start_dates, 
		     const vector<double> * start_rates, const vector<double> * start_durations);
	void set_freeparams(int nump, bool lf, vector<int> * freep, vector<double> * params);
	int calc_isum(vector<int> & inc);
	void set_log_pen(bool); //set whether log penalty is set
	//these should all be hidden but oh well for now
	int numparams;
	bool cvnodeset;
	vector<int> * get_cv_nodes();
	void set_cv_nodes(vector<int> & cvnodes);
	vector<double> dates;
	vector<double> rates;
	vector<double> durations;
	vector<int> freeparams;
	const vector<int> * free;
	const vector<double> * min;
	const vector<double> * max;
	const vector<double> * pen_min;//just the mins that are calculated in the penalty
	const vector<double> * pen_max;//just the maxs that are calculated in the penalty
	const vector<int> * parents_nds_ints;
	const vector<int> * child_counts;
	const vector<double> * char_durations;
	const vector<double> * log_fact_char_durations;
	const vector< vector<int> > * children_vec;
	int numnodes;
	int numsites;
	double smoothing;
	bool lf; 
	void delete_ad_arrays();
	bool paramverbose;
	double minrate;
	double ftol;
	bool isConstrained;
	double penalty_boundary;
};

#endif /* PL_CALC_PARALLEL_H_ */

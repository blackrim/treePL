/*
 * siman_calc.h
 *
 *  Created on: Feb 22, 2010
 *      Author: smitty
 */

#ifndef SIMAN_CALC_PAR_H_
#define SIMAN_CALC_PAR_H_

#include "pl_calc_parallel.h"

class siman_calc_par{
private:
	double temperature;
	double cool;
	double cur_like;
	pl_calc_parallel * pl_;
	double step_size_rt;
	double step_size_dt;
	int maxit;
	double stop_temp;
	vector<double> cur_values;
	vector<double> test_values;
	vector<double> last_print_values;
	void make_step_plcp(double, int);
	vector<double> index_probs; //this is to determine which index
	int rate_accepts;
	int rate_trials;
	int date_accepts; 
	int date_trials;
	double rate_accept;double rate_trial;
	double date_accept;double date_trial;
	bool converged;
	double calc_difference_between_vecs(vector<double> &,vector<double> &, int, int);
	bool verbose;
	double calc_sum_part_d(int, int,int);

public:
	siman_calc_par();
	void set_verbose(bool v);
	void set_pl(pl_calc_parallel * plcp,double temp,double cl, double rt_step, double dt_step, int derit, double sttemp);
	void optimize(vector<double> & init_params);
	void mapspace(pl_calc_parallel * plp,vector<double> & init_params);
};

#endif /* SIMAN_CALC_H_ */

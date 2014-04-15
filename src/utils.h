
#ifndef UTILS_H_
#define UTILS_H_


#include <string>
#include <vector>
#include <map>
#include "pl_calc_parallel.h"
#include "optim_options.h"
#include "tree.h"

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);
void TrimSpaces( string& str);
void set_mins_maxs(Tree * tr,map<Node *,double> * mins, map<Node *,double> * maxs);
double get_start_rate(Tree * tree, vector<double> * durations);
void process_initial_branch_lengths(Tree * tree, bool collapse, int numsites);
void apply_node_label_preorder(Tree * tr);
void calc_char_durations(Tree * tr, int numsites);
void setup_date_constraints(Tree * tr, map<Node *, double> * inmins, map<Node *, double> * inmaxs);
void extract_tree_info(Tree * tr, vector<int> * free,
		       vector<int> * parent_nds_ints, vector<int> * child_c, vector<double> * char_durations, 
		       vector<double> * log_fact_char_durations, vector<double> * min,
		       vector<double> * max, vector< vector<int> > * children_vec,vector<double>* pen_min, vector<double>* pen_max);
double logFact(double k);
int generate_param_order_vector(vector<int> * freeparams, bool lf, 
				vector<int> * cvnodes, vector<int> * free);
void get_feasible_start_dates(Tree * tr, vector<double> * dates, vector<double> * rates,
			      vector<double> * durations);
int a_feasible_time(Node * node,double timeAnc);
double round(double x);
double myRand(void);
float rand_float_range(float a, float b);
void set_node_order(Tree * tr);
int optimize_best(int whichone, bool ad, double * init_x, pl_calc_parallel * plp, int numiter, bool moredetail,double ftol,double xtol);
int optimize_best_parallel(int whichone, double * init_x, pl_calc_parallel * plp, int numiter, bool moredetail,double ftol,double xtol);
void prime_optimization(pl_calc_parallel & plp, vector<double> & params);

bool check_possible_cv_nodes(const int cvnode,const  vector<int> & samps, const vector<int> & parents);
double get_median(vector<double> & container);
double get_gmean(vector<double> & container);
double get_sum(vector<double> & container);
void process_ind8s_inr8s(string ind8s, string inr8s,vector<double> * params, Tree * intree);
bool optimize_full(pl_calc_parallel & plp, vector<double> * params, const OptimOptions * oopt, bool cv);
void mapspace(pl_calc_parallel &plp, vector<double> * params);
int count_trees (string & treefile);
#endif

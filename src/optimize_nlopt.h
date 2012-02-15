
#ifndef OPTIMIZE_NLOPT_H_
#define OPTIMIZE_NLOPT_H_

#include "pl_calc_parallel.h"

double function_plcp_nlopt_ad(unsigned n, const double *x, double *g, void *state);
double function_plcp_nlopt(unsigned n, const double *x, double *g, void *state);
double function_plcp_nlopt_nograd(unsigned n, const double *x, double *g, void *state);
double function_plcp_nlopt_ad_parallel(unsigned n, const double *x, double *g, void *state);
int optimize_plcp_nlopt(double * init_x, pl_calc_parallel *pl_, int numiter, int whichone, bool ad, bool moredetail,double ftol,double xtol);
int optimize_plcp_nlopt_ad_parallel(double * init_x, pl_calc_parallel *pl_,int numiter, int whichone, bool moredetail, double ftol,double xtol);
#endif

/*
 * optimize_tnc.h
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */

#ifndef OPTIMIZE_TNC_H_
#define OPTIMIZE_TNC_H_

#include "tnc.h"
#include "pl_calc_parallel.h"

static tnc_function function_plcp_ad;
static tnc_function function_plcp_ad_parallel;
static tnc_function function_plcp;
int optimize_plcp_tnc_ad_parallel(double *init_x, pl_calc_parallel *pl,double ftol);
int optimize_plcp_tnc(double *init_x, pl_calc_parallel *pl_, int numiter,bool ad,double ftol);

#endif /* OPTIMIZE_TNC_H_ */

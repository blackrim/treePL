
#include "optim_options.h"

OptimOptions::OptimOptions(){
    lftemp = 50;
    lfcool = 0.999;
    lfstoptemp = 0.001;
    plstoptemp = 0.0001;
    pltemp = 20;
    plcool = 0.99975;
    lfrtstep = 0.07;
    lfdtstep = 0.25;
    plrtstep = 0.1;
    pldtstep = 0.1;
    LARGE =10e+15;
    thorough = false;
    lfiter = 3;
    pliter = 5; 
    cviter = 3;
    lfsimaniter = 15000;
    plsimaniter = 15000;
    cvsimaniter = 7500;
    calcgrad = false;
    verbose = false;
    
    bestopt = 0; //0 = tnc, 1 = lbfgs, 2 = tnewton, 3 = mma
    bestadopt = 0; //0 = tnc, 1 = lbfgs, 2 = tnewton, 3 = mma
    bestcvopt = 3; //0 = tnc, 1 = lbfgs, 2 = tnewton, 3 = mma, these are all ad
    moredetail = false;
    moredetailad = false;
    moredetailcvad = false;
    ftol = 1e-10;
    xtol = -1;
    maxoptimiters = 10000;
}

OptimOptions::~OptimOptions(){}

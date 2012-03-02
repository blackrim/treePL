#ifndef OPTIM_OPTIONS_H_
#define OPTIM_OPTIONS_H_

/*
 * This is just optimization options storage.
 * Could be a struct, but oh well
 */

class OptimOptions{
public:
    OptimOptions();
    double lftemp;
    double lfcool;
    double lfstoptemp;
    double plstoptemp;
    double pltemp;
    double plcool;
    double lfrtstep;
    double lfdtstep;
    double plrtstep;
    double pldtstep;
    double LARGE;
    bool thorough;
    int lfiter;
    int pliter; 
    int cviter;
    int lfsimaniter;
    int plsimaniter;
    int cvsimaniter;
    bool calcgrad;
    bool verbose;

    int bestopt; //0 = tnc, 1 = lbfgs, 2 = tnewton, 3 = mma
    int bestadopt; //0 = tnc, 1 = lbfgs, 2 = tnewton, 3 = mma
    int bestcvopt; //0 = tnc, 1 = lbfgs, 2 = tnewton, 3 = mma, these are all ad
    bool moredetail;
    bool moredetailad;
    bool moredetailcvad;
    double ftol;
    double xtol;
    int maxoptimiters;
    int nthreads;
    ~OptimOptions();

};

#endif /* OPTIM_OPTIONS_H_ */

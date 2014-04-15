/*
 * siman_calc.cpp
 *
 *  Created on: Feb 22, 2010
 *      Author: smitty
 */

#include "siman_calc_par.h"
#include "pl_calc_parallel.h"
#include "optimize_tnc.h"

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>

#define LARGE 10e+13

siman_calc_par::siman_calc_par():index_probs(vector<double>()), verbose(false){

}

void siman_calc_par::set_verbose(bool v){
    verbose = v;
}

void siman_calc_par::set_pl(pl_calc_parallel * pl, double temp, double cl, double rt_step, double dt_step, int inmaxit, double sttemp){
    temperature=  temp;
    cool = cl;
//    if (pl->smoothing <= 1){
    step_size_rt = rt_step;
    step_size_dt = dt_step;//0.0001 for 1
//    }else{
//    step_size_rt = 0.00001;//0.0025;
//    step_size_dt = 0.075;//0.0025;
//   }
    stop_temp = sttemp;
    maxit = inmaxit;
    //step_size = 0.0000025;
    pl_ = pl;

    date_accepts = 0;
    date_trials = 0;
    rate_accepts = 0;
    rate_trials = 0;
    
    index_probs.clear();    
    /*
     * setting up the proportion of choosing particular indices
     */
    //make sure that this is correct given the number of cv
    int fcount = 0;
    if(pl_->lf == true){
        cout << "LF SIM" << endl;
        //fcount = pl_->rates.size();
        for(int i=0; i < pl_->numparams;i++){
            if (i == 0){
                index_probs.push_back(0.01);
            }else{
                index_probs.push_back(0.99/(pl_->numparams-1));
            }
        }
    }else{
        for(int i=0;i<pl_->rates.size();i++){
            if(pl_->freeparams[fcount] != -1){
                index_probs.push_back(0.5/(pl_->numnodes-1));
            }
            fcount += 1;
        }
        for(int i=0; i < pl_->dates.size();i++){
            if (pl_->freeparams[fcount] != -1){
                index_probs.push_back(0.5/(pl_->numparams-(pl_->numnodes-1)));
            }
            fcount += 1;
        }
    }
}

void siman_calc_par::optimize(vector<double> &init_params){
    if (index_probs.size() != pl_->numparams){
        cout << "something wrong with the number of params" << endl;
        exit(0);
    }
    cur_values = vector<double>(init_params);
    last_print_values = cur_values;
    cur_like = pl_->calc_pl(cur_values);
    //cout << "start sim like: " << cur_like << endl;
    test_values = vector<double>(init_params);
    converged = false;
    //srand( time(NULL) );        // init the random generator. no - this is done in main.
    int ind,iteration = 0;double sa,r,step,test_like;
    int same_value = 0;
    while (temperature > stop_temp){
        ind = rand() % pl_->numparams;// choose one index to change
        //proper way to do this when there are different probs
        /*
         * get a random double
         */

        double id = rand()/(double)RAND_MAX;
        bool kp = true;
        unsigned int tid = 0;
        double cur = 0;
        while(kp){
    //        cout << id << " " << index_probs[tid] << " " << tid << " " << cur << endl;
            if(id >= cur){
                if((tid + 1) == index_probs.size() ){
                    ind = tid;
                    kp = false;
                    break;
                }else if( id < (cur+index_probs[tid]) ){
                    ind = tid;
                    kp = false;
                    break;
                }
            }
            cur += index_probs[tid];
            tid += 1;
        }
    //    cout << "index: " << ind << endl;
        /*
         * end choosing index
         */
        int index_type; //0 = rate, 1 = date
        double step_size;
        if (pl_->lf == true){
            if(ind == 0){
                index_type = 0;
                step_size = step_size_rt;
                rate_trials += 1;
            }else{
                index_type = 1;
                step_size = step_size_dt;
                date_trials += 1;
            }
        }else{
            if(ind < pl_->numnodes-1){
                index_type = 0;
                step_size = step_size_rt;
                rate_trials += 1;
            }else{
                index_type = 1;
                step_size = step_size_dt;
                date_trials += 1;
            }
        }

        double rna = (rand() / (RAND_MAX + 1.0));
        step = (2 * step_size * rna) - step_size;

        make_step_plcp(step,ind);

        test_like = pl_->calc_pl(test_values);

    //    cout << "tl " << iteration << " : " << test_like << endl;
        if(test_like >= LARGE){//if the return is inf then just continue without counting an iteration
            test_values[ind] = cur_values[ind];
            converged = false;
            continue;
        }
        sa = exp((-cur_like-test_like)/temperature);
        if (test_like<cur_like){
            if((cur_like-test_like) > 0.01){
                same_value = 0;
            }else{
                same_value += 1;
            }
            cur_values[ind] = test_values[ind];
            cur_like = test_like;
            if(index_type == 0) {
                rate_accepts += 1;
            }else{
                date_accepts += 1;
            }
        } // simulated annealing
        else if ((r = (rand()/(RAND_MAX + 1.0))) < sa){
            cur_like = test_like;
            cur_values[ind] = test_values[ind];
        }else{
            converged = false;
            test_values[ind] = cur_values[ind];
        }
        iteration++;
        if(iteration % 100 == 0){
            temperature *= cool;          // reduce the temperature
        }
        if(iteration % 100 == 0){
            //   cout << date_accepts << " " << date_accepts/float(date_trials) << " " << float(date_trials) << " " << step_size_dt<< " " << step_size_dt * 0.99<< endl;
            if((rate_accepts/float(rate_trials)) < 0.15)
            step_size_rt *= 0.99;
            else
            step_size_rt *= 1.01;
            if((date_accepts/float(date_trials)) < 0.15)
            step_size_dt *= 0.99;
            else
            step_size_dt *= 1.01;
        }
        if(iteration % 1000 == 0){
            if(verbose){
            cout << iteration<< " cur_like: " << cur_like << " test_like: "<< test_like<<
                " same_value: "<< same_value<< " temp: "<< temperature << " rate_acc: ("<< step_size_rt<< ") " << rate_accepts/float(rate_trials) <<
                " date_acc: ("<< step_size_dt << ") " << date_accepts/float(date_trials)  << " vecdiff: " << calc_difference_between_vecs(cur_values,last_print_values,0,1) <<
                "," << calc_difference_between_vecs(cur_values,last_print_values,1,cur_values.size()) << endl;
            }
            last_print_values = cur_values;
        }
        if(same_value >= 5000){
            if(converged == false){
                if(verbose){
                    cout << "converged" << endl;
                }
                converged = true;
                same_value = 0;
            }
            if((converged == true && same_value >= 5000) || iteration > maxit){
                break;
            }
        }
        if (iteration > maxit){
            break;
        }
    }// while
    for(int i=0;i<init_params.size();i++){
        init_params[i] = cur_values[i];
    }
    //cout << "end siman: " << cur_like << endl;
    //exit(0);
}



//index should be the index of the free params
void siman_calc_par::make_step_plcp(double step_size, int index){
    //should have a step size per variable
//step size should be (1.+step_size) * test_values[index]
    //rate or date
//need to make sure that this works with the cv
    
    //   cout << "dt: " << index  << " " << test_values[index] << " " << (1.+step_size)+test_values[index] << " " << step_size << endl;
test_values[index] = (1.+step_size)*test_values[index];
    /*  
    if((pl_->lf == true && index == 0) || (pl_->lf == false && index < pl_->numnodes-1)){//rate
        test_values[index] = (1.+step_size)*test_values[index];
    }else{//date
    double mi = -1;
    double mx = -1;
    int fcount;
    if (pl_->lf == true)
        fcount = 1;
    else
        fcount = pl_->numnodes-1;//TODO : make sure this is right with cv
    for(int j=0;j<pl_->numnodes;j++){
        if((*pl_->get_cv_nodes())[j] == 0 && pl_->freeparams[pl_->numnodes+j] != -1){
            if(index == fcount+j){
                mi = pl_->min->at(j);//nd->min;//pl_->mins[nd];
                mx = pl_->max->at(j);//nd->max;//pl_->maxs[nd];
            }
        }
    }
    double st = (1.+step_size)*test_values[index];
//        double st = (step_size)+cu;
    if (mx != -1){
        if(st > mx)
        st = mx;
    }
    if(mi != -1 ){
        if(st < mi)
        st = mi;
    }
    test_values[index] = st;
//    cout << index << " " << test_values[index] <<" " << mx << " " << mi  << " " << endl;
    //exit(0);
    }*/

}


double siman_calc_par::calc_difference_between_vecs(vector<double> & vec1,vector<double> & vec2, int start, int end){
    double ret = 0;
    for(int i=start;i<end;i++){
        ret += fabs(vec1[i]-vec2[i]);
    }
    return ret;
}


void siman_calc_par::mapspace(pl_calc_parallel * plp,vector<double> & init_params){
    
    maxit = 10000;
    pl_ = plp;
    step_size_rt = 0.0001;
    step_size_dt = 0.0001;
    date_accepts = 0;
    date_trials = 0;
    rate_accepts = 0;
    rate_trials = 0;
    verbose = true;
    index_probs.clear();
    ofstream mapspacefile;
    mapspacefile.open("mapspace.txt",ios::out);
    mapspacefile.precision(10);
    /*
     * setting up the proportion of choosing particular indices
     */
    int fcount = 0;
    for(int i=0;i<pl_->rates.size();i++){
        if(pl_->freeparams[fcount] != -1){
            index_probs.push_back(0.5/(pl_->numnodes-1));
        }
        fcount += 1;
    }
    for(int i=0; i < pl_->dates.size();i++){
        if (pl_->freeparams[fcount] != -1){
            index_probs.push_back(0.5/(pl_->numparams-(pl_->numnodes-1)));
        }
        fcount += 1;
    }

    if (index_probs.size() != pl_->numparams){
        cout << "something wrong with the number of params" << endl;
        exit(0);
    }
    cur_values = vector<double>(init_params);
    last_print_values = cur_values;
    cur_like = pl_->calc_pl(cur_values);
    double INLIKE = cur_like;
    cout << "start map like: " << cur_like << endl;
    test_values = vector<double>(init_params);
    //srand( time(NULL) );        // init the random generator. no - this is done in main.
    int ind,iteration = 0;double sa,r,step,test_like;
    int same_value = 0;
    while (iteration < maxit){
        ind = rand() % pl_->numparams;// choose one index to change
        //proper way to do this when there are different probs
        /*
         * get a random double
         */

        double id = rand()/(double)RAND_MAX;
        bool kp = true;
        unsigned int tid = 0;
        double cur = 0;
        while(kp){
            if(id >= cur){
                if((tid + 1) == index_probs.size() ){
                    ind = tid;
                    kp = false;
                    break;
                }else if( id < (cur+index_probs[tid]) ){
                    ind = tid;
                    kp = false;
                    break;
                }
            }
            cur += index_probs[tid];
            tid += 1;
        }
        /*
         * end choosing index
         */
        int index_type; //0 = rate, 1 = date
        double step_size;
        if(ind < pl_->numnodes-1){
            index_type = 0;
            step_size = step_size_rt;
            rate_trials += 1;
        }else{
            index_type = 1;
            step_size = step_size_dt;
            date_trials += 1;
        }

        double rna = (rand() / (RAND_MAX + 1.0));
        step = (2 * step_size * rna) - step_size;

        make_step_plcp(step,ind);

        test_like = pl_->calc_pl(test_values);

        if(test_like >= LARGE){//if the return is inf then just continue without counting an iteration
            test_values[ind] = cur_values[ind];
            converged = false;
            continue;
        }
        if ((test_like-2) < cur_like){
            if((cur_like-test_like) > 0.01){
                same_value = 0;
            }else{
                same_value += 1;
            }
    //        cur_values[ind] = test_values[ind];
    //        cur_like = test_like;
            if(index_type == 0){
                rate_accepts += 1;
            }else{
                date_accepts += 1;
            }
            //write to the file
            mapspacefile << test_like;
            for(int n = 0; n < pl_->dates.size();n++){
                if(pl_->dates[n] != 0){
                    mapspacefile << "\t" <<  pl_->dates[n];
                }
            }
            mapspacefile << endl;
        }else{
            converged = false;
            test_values[ind] = cur_values[ind];
        }
        iteration++;
        if(iteration % 100 == 0){
            //   cout << date_accepts << " " << date_accepts/float(date_trials) << " " << float(date_trials) << " " << step_size_dt<< " " << step_size_dt * 0.99<< endl;
            if((rate_accepts/float(rate_trials)) < 0.15){
                step_size_rt *= 0.99;
            }else{
                step_size_rt *= 1.01;
            }
            if((date_accepts/float(date_trials)) < 0.15){
                step_size_dt *= 0.99;
            }else{
                step_size_dt *= 1.01;
            }
        }
        if(iteration % 1000 == 0){
            if(verbose){
            cout << iteration<< " cur_like: " << cur_like << " test_like: "<< test_like<<
                " same_value: "<< same_value << " rate_acc: ("<< step_size_rt<< ") " << rate_accepts/float(rate_trials) <<
                " date_acc: ("<< step_size_dt << ") " << date_accepts/float(date_trials)  << " vecdiff: " << calc_difference_between_vecs(cur_values,last_print_values,0,1) <<
                "," << calc_difference_between_vecs(cur_values,last_print_values,1,cur_values.size()) << endl;
            }
            last_print_values = cur_values;
        }
    }// while
    mapspacefile.close();
    for(int i=0;i<init_params.size();i++){
        init_params[i] = cur_values[i];
    }
}

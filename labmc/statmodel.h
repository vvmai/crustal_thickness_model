/*
 * statmodel.h
 *
 * Jun Korenaga
 * Summer 2008
 */

#ifndef _LABMC_STATMODEL_H_
#define _LABMC_STATMODEL_H_

#include <map>
#include <iostream>
#include "array.h"
#include "parameter.h"
#include "state.h"
#include "flowlaw.h"

class StatModel {
public:
    StatModel(); 
    ~StatModel();

    void addParallel(FlowLaw* f) { para.push_back(f); }
    void addSequential(FlowLaw* f) { seq.push_back(f); }
    int numParallel() const { return para.size(); }
    int numSequential() const { return seq.size(); }
    const FlowLaw* parallel(int i) const { return para(i); }
    const FlowLaw* sequential(int i) const { return seq(i); }

    void readData(char *fn);
    void stateToRead(int i);
    void printRunIds(ostream& os);
    int numData() const { return data.size(); }
    void setupBiasCorrection(double, double); 
    void randomizeData();

    void setUp(); 
    void addConstraints(const Array1d<int>& ifixed); 
    int numUnfixedParams() const { return unfixed_params.size(); }
    Parameter* unfixedParam(int i) { return unfixed_params(i); }
    double calcChiSq(int k, double& chi2_orig);
    void useConjugateGradient(int ntrial=10, int imax=20, 
			      double tol1=1e-4, double tol2=1e-7) 
	{ do_CG = true; cg_ntrial = ntrial; cg_iter_max = imax;
	    cg_tol = tol1; brent_tol = tol2; }
private:
    typedef double (StatModel::*PF1DIM)(double);

    void fixParams(); 
    void beforeBestFitA(int k);
    bool calcBestFitA();
    bool calcBestFitA_CG();
    void afterBestFitA();
    bool removeNonPositive(Array1d<double>& v);
    void setA(const Array1d<double>& logA_);
    double CG_calc_cost_and_grad(bool needGrad, 
				 Array1d<double>& grad);
    double CG_line_min();
    double f1dim(double x);
    void mnbrak(double *ax, double *bx, double *cx,
		double *fa, double *fb, double *fc, PF1DIM);
    double brent(double ax, double bx, double cx, double *xmin, PF1DIM);
    void setMoreState(State& d, int i, double v, double dv);
    double randomizeValue(double val, double dval);
  
    Array1d<FlowLaw*> para; 
    Array1d<FlowLaw*> seq;
    Array1d<int> i_para_scaling;
    Array1d<int> i_para_noscaling;
    Array1d<int> i_seq_scaling;
    Array1d<int> i_seq_noscaling;

    bool is_constrained;
    Array1d< Array1d<int> > isrc;

    Array1d<State> orig_data;
    Array1d<State> data;
    Array1d<int> more_state;
    Array1d<int> runid;
    map<int,int> uid;

    Array2d<double> fmat_para, fmat_seq;
    Array2d<double> edot_para, edot_seq;
    Array1d<double> bvec, tmp_bvec, inv_dedot, alpha1, alpha2;
    Array2d<double> Fmat1, Fmat2;
    Array2d<double> AtA1, AtA2, Atb1, Atb2; 
    Array1d<double> chivec;

    bool do_CG;
    Array1d<double> logA, tmp_logA, grad, new_grad, direc;
    Array1d<double> maxlogA, minlogA, bestlogA;
    int cg_iter_max, cg_ntrial;
    double cg_tol, brent_tol, data_lognorm2;
    long cg_idum, data_idum; 

    bool do_bias;
    double min_bias, max_bias;
    Array1d<Parameter*> bias;

    Array1d<Parameter*> unfixed_params;
    Array1d<int> param_type;
    enum {type_para, type_seq, type_bias};
    Array1d<int> param_id;

    static const double max_chi2;
};

#endif /* _LABMC_STATMODEL_H_ */

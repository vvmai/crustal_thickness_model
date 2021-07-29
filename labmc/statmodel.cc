/*
 * statmodel.cc
 *
 * Jun Korenaga
 * Summer 2008
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctime>
#include "labmc_nr.h"
#include "statmodel.h"
#include "util.h"

const double StatModel::max_chi2 = 1e30;

StatModel::StatModel()
{
    do_bias = false;
    do_CG = false;
    is_constrained = false;
}

StatModel::~StatModel()
{
    for (int i=1; i<=para.size(); i++){
	delete para(i);
    }
    for (int i=1; i<=seq.size(); i++){
	delete seq(i);
    }
}

void StatModel::readData(char *fn)
{
    int ndata = countLines(fn);
    orig_data.resize(ndata);
    data.resize(ndata);
    runid.resize(ndata);

    istream* pin = new ifstream(fn);
    for (int i=1; i<=ndata; i++){
	double T, dT, p, dp, e, de, sig, dsig;
	double id;
	*pin >> T >> dT >> p >> dp >> e >> de >> sig >> dsig;
	orig_data(i).setTemperature(T,dT);
  	orig_data(i).setPressure(p*1e9,dp*1e9);
	orig_data(i).setStrainRate(e,de);
	orig_data(i).setStress(sig,dsig);

	for (int j=1; j<=more_state.size(); j++){
	    double val,dval;
	    *pin >> val >> dval;
	    setMoreState(orig_data(i),more_state(j),val,dval);
	}

	*pin >> id;
	runid(i) = int(id);
    }

    data = orig_data;
}

void StatModel::randomizeData()
{
    for (int i=1; i<=orig_data.size(); i++){
	double T = orig_data(i).temperature();
	double dT = orig_data(i).temperatureError();
	double p = orig_data(i).pressure();
	double dp = orig_data(i).pressureError();
	double e = orig_data(i).strainRate();
	double de = orig_data(i).strainRateError();
	double sig = orig_data(i).stress();
	double dsig = orig_data(i).stressError();

	double newT = randomizeValue(T,dT);
	double newp = randomizeValue(p,dp);
	double newe = randomizeValue(e,de);
	double newsig = randomizeValue(sig,dsig);
	
	data(i).setTemperature(newT,dT);
	data(i).setPressure(newp,dp);
	data(i).setStrainRate(newe,de);
	data(i).setStress(newsig,dsig);

	double val, dval, newval;
	for (int j=1; j<=more_state.size(); j++){
	    switch(more_state(j)){
	    case State::grain_size:
		val = orig_data(i).grainSize();
		dval = orig_data(i).grainSizeError();
		newval = randomizeValue(val,dval);
		data(i).setGrainSize(newval,dval);
		break;
	    case State::water_content:
		val = orig_data(i).waterContent();
		dval = orig_data(i).waterContentError();
		newval = randomizeValue(val,dval);
		data(i).setWaterContent(newval,dval);
		break;
	    case State::oxygen_fugacity:
		val = orig_data(i).oxygenFugacity();
		dval = orig_data(i).oxygenFugacityError();
		newval = randomizeValue(val,dval);
		data(i).setOxygenFugacity(newval,dval);
		break;
	    default:
		error("StatModel::randomizeData - invalid istate");
		break;
	    }
	}
    }

    beforeBestFitA(-1); // recalculate with newly randomized data
}

double StatModel::randomizeValue(double v, double dv)
{
//    double newv = v + dv*gasdev(&data_idum);
    double newv = v + dv*2.0*(ran2(&data_idum)-0.5);
//    cerr << "rG: " << v << " " << dv << " " << newv << '\n';
    return newv;
}

void StatModel::printRunIds(ostream& os)
{
    os << "runids: ";
    for (int i=1; i<=runid.size(); i++){
	os << runid(i) << " ";
    }
    os << '\n';
}

void StatModel::stateToRead(int i)
{
    more_state.push_back(i);
}

void StatModel::setMoreState(State& d, int istate, 
			     double val, double dval)
{
    switch(istate){
    case State::grain_size:
	d.setGrainSize(val,dval);
	break;
    case State::water_content:
	d.setWaterContent(val,dval);
	break;
    case State::oxygen_fugacity:
	d.setOxygenFugacity(val,dval);
	break;
    default:
	error("StatModel::setMoreState - invalid istate");
	break;
    }
}
void StatModel::setupBiasCorrection(double bmin, double bmax)
{
    do_bias = true;
    min_bias = bmin;
    max_bias = bmax;
}

void StatModel::setUp()
{
    //
    // check if all of required states are given
    //
    for (int i=1; i<=para.size(); i++){
	if (para(i)->checkStates(more_state)==false){
	    error("StatModel::setUp - missing state(s)");
	}
    }
    for (int i=1; i<=seq.size(); i++){
	if (seq(i)->checkStates(more_state)==false){
	    error("StatModel::setUp - missing state(s)");
	}
    }

    //
    // classify flow laws
    //
    i_para_scaling.resize(0);
    i_para_noscaling.resize(0);
    i_seq_scaling.resize(0);
    i_seq_noscaling.resize(0);

    for (int i=1; i<=para.size(); i++){
	if (para(i)->needScaling()){
	    i_para_scaling.push_back(i);
	}else{
	    i_para_noscaling.push_back(i);
	}
    }
    for (int i=1; i<=seq.size(); i++){
	if (seq(i)->needScaling()){
	    i_seq_scaling.push_back(i);
	}else{
	    i_seq_noscaling.push_back(i);
	}
    }

    // allocate relevant arrays
    fmat_para.resize(orig_data.size(),para.size());
    fmat_seq.resize(orig_data.size(),seq.size());
    edot_para.resize(orig_data.size(),para.size());
    edot_seq.resize(orig_data.size(),seq.size());
    Fmat1.resize(orig_data.size(),i_para_scaling.size());
    Fmat2.resize(orig_data.size(),i_seq_scaling.size());
    AtA1.resize(i_para_scaling.size(),i_para_scaling.size());
    AtA2.resize(i_seq_scaling.size(),i_seq_scaling.size());
    Atb1.resize(i_para_scaling.size(),1);
    Atb2.resize(i_seq_scaling.size(),1);
    alpha1.resize(i_para_scaling.size());
    alpha2.resize(i_seq_scaling.size());
    chivec.resize(orig_data.size());
    int n_scaling=i_para_scaling.size()+i_seq_scaling.size();
    logA.resize(n_scaling);
    maxlogA.resize(n_scaling);
    minlogA.resize(n_scaling);
    bestlogA.resize(n_scaling);
    tmp_logA.resize(n_scaling);
    grad.resize(n_scaling);
    new_grad.resize(n_scaling);
    direc.resize(n_scaling);

    //
    // count unfixed parameters
    //
    unfixed_params.resize(0);
    param_type.resize(0);
    param_id.resize(0);

    for (int i=1; i<=para.size(); i++){
	for (int j=1; j<=para(i)->numUnfixedParams(); j++){
	    unfixed_params.push_back(para(i)->unfixedParam(j));
	    param_type.push_back(type_para);
	    param_id.push_back(i);
	}
    }
    for (int i=1; i<=seq.size(); i++){
	for (int j=1; j<=seq(i)->numUnfixedParams(); j++){
	    unfixed_params.push_back(seq(i)->unfixedParam(j));
	    param_type.push_back(type_seq);
	    param_id.push_back(i);
	}
    }

    if (do_bias){
	Array1d<int> tmpid(runid.size());
	tmpid = runid;
	tmpid.sort_unique();
	
	for (int i=1; i<=tmpid.size(); i++){
	    uid[tmpid(i)] = i;
	}
	
	bias.resize(tmpid.size());
	for (int i=1; i<=bias.size(); i++){
	    char str[MaxStr];
	    sprintf(str,"Bias[%d]",i);
	    if (i==1){
		bias(i) = new Parameter(0.0,0.0,str); 
	    }else{
		bias(i) = new Parameter(min_bias,max_bias,str);	   
	    }
	    bias(i)->randomize();
	    unfixed_params.push_back(bias(i));
	    param_type.push_back(type_bias);
	    param_id.push_back(i);
	}
    }

    //
    // data-related pre-calculation
    //
    inv_dedot.resize(orig_data.size());
    bvec.resize(orig_data.size());
    tmp_bvec.resize(orig_data.size());
    data_lognorm2=0.0;
    for (int i=1; i<=orig_data.size(); i++){
	double val = 1.0/data(i).strainRateError();
	inv_dedot(i) = val;
	bvec(i) = data(i).strainRate()*val;
	
	double val2 = log(data(i).strainRate());
	data_lognorm2 += val2*val2;
    }

    isrc.resize(unfixed_params.size());
    for (int i=1; i<=isrc.size(); i++){
	isrc(i).resize(0);
	isrc(i).push_back(i);
    }

    //
    // initialize chi2 calculation
    //
    beforeBestFitA(-1);

    cg_idum = long(-(abs(data_lognorm2*1000)+1)); // for ran2();
    data_idum = long(-(abs(data_lognorm2*1000)+2)); // for ran2();
}

void StatModel::addConstraints(const Array1d<int>& ifixed)
{
    for (int i=1; i<=ifixed.size(); i++){
	if (ifixed(i)>0) isrc(ifixed(i)).push_back(i);
    }
    is_constrained=true;
    fixParams();
}

void StatModel::fixParams()
{
    if (is_constrained){
	for (int i=1; i<=isrc.size(); i++){
	    if (isrc(i).size()>1){
		for (int j=2; j<=isrc(i).size(); j++){
		    unfixed_params(isrc(i)(j))->setValue(unfixed_params(isrc(i)(1))->value());
		}
	    }
	}
    }
}

void StatModel::beforeBestFitA(int k)
{
    if (k<0){
	for (int i=1; i<=para.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		fmat_para(j,i) = para(i)->prepPredict(data(j));
	    }      
	}

	for (int i=1; i<=seq.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		fmat_seq(j,i) = seq(i)->prepPredict(data(j));
	    }      
	}
    }else{
	for (int i=1; i<=isrc(k).size(); i++){
	    int kk = isrc(k)(i);
	    switch(param_type(kk)){
	    case type_para:
	    {
		int i = param_id(kk);
		for (int j=1; j<=data.size(); j++){
		    fmat_para(j,i) = para(i)->prepPredict(data(j));
		}      
		break;
	    }
	    case type_seq:
	    {
		int i = param_id(kk);
		for (int j=1; j<=data.size(); j++){
		    fmat_seq(j,i) = seq(i)->prepPredict(data(j));
		}      
		break;
	    }
	    case type_bias:
		// do nothing
		break;
	    }
	}
    }
}

bool StatModel::calcBestFitA()
{
    // 1. para_scaling==0 && seq_scaling==0 
    //    (this is not an impossible case. other parameters can vary)
    // 2. para_scaling>0 && seq_scaling>0 
    //    ---> do_CG
    // 3. para_scaling>0 && seq_scaling==0
    //    ---- but seq_nosclaing can still be non-zero
    // 4. para_scaling==0 && seq_scaling>0 
    //    ---- but para_noscaling can still be non-zero

    if (i_para_scaling.size()==0 && i_seq_scaling.size()==0){
	// need to do nothing
	return true;
    }

    if (do_CG) return calcBestFitA_CG();
    if (i_para_scaling.size()>0 && i_seq_scaling.size()>0){
	if (!do_CG){
	    cerr << "do_CG needs to be set\n"; 
	    exit(1);
	}
    }

    if (i_para_scaling.size()>0){
	// i_seq_scaling.size() == 0 but
	// i_seq_noscaling.size() may not be zero.
	//
	// edot = x1*f1+x2*f2+...+e1+e2+... + 1/(1/g1+1/g2+...)
	// --> edot-(e1+e2+...)-(1/(1/g1+..)) = x1*f1+x2*f2+...
	// --> tmp_bvec = Fmat*xvec

	// prepare tmp_bvec
	tmp_bvec = bvec; 
	if (do_bias){
	    for (int j=1; j<=data.size(); j++){
		tmp_bvec(j) *= exp(bias(uid[runid(j)])->value()*(-1.0));
	    }
	}

	for (int i=1; i<=i_para_noscaling.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		tmp_bvec(j) -= 
		    fmat_para(j,i_para_noscaling(i))*inv_dedot(j);
	    }
	}

	if (i_seq_noscaling.size()>0){
	    for (int j=1; j<=data.size(); j++){
		double val=0.0;
		for (int i=1; i<=i_seq_noscaling.size(); i++){		
		    val+=1.0/fmat_seq(j,i_seq_noscaling(i));
		}
		tmp_bvec(j) -= inv_dedot(j)/val;
	    }
	}
	if (removeNonPositive(tmp_bvec)==false) return false;

	// preprare Fmat
	for (int i=1; i<=i_para_scaling.size(); i++){
	    double val=0.0;
	    for (int j=1; j<=data.size(); j++){
		Fmat1(j,i) = 
		    fmat_para(j,i_para_scaling(i))*inv_dedot(j);
		val += log(Fmat1(j,i));
	    }
	    alpha1(i) = exp(val/data.size());
	}
	for (int i=1; i<=i_para_scaling.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		Fmat1(j,i) /= alpha1(i);
	    }
	}

	// solve linear system
	for (int k=1; k<=i_para_scaling.size(); k++){
	    for (int l=1; l<=k; l++){
		AtA1(k,l) = 0.0;
		for (int i=1; i<=data.size(); i++){
		    AtA1(k,l) += Fmat1(i,k)*Fmat1(i,l);
		}
		AtA1(l,k) = AtA1(k,l);
	    }
	    Atb1(k,1) = 0.0;
	    for (int i=1; i<=data.size(); i++){
		Atb1(k,1) += Fmat1(i,k)*tmp_bvec(i);
	    }
	}

	gaussj(AtA1.toRecipe(), i_para_scaling.size(), Atb1.toRecipe(), 1);
	// check if there's a negative scaling constant...
	for (int i=1; i<=i_para_scaling.size(); i++){
	    if (Atb1(i,1)<0) return false;
	}

	// post-scaling 
	double tmpsum=0.0;
	for (int i=1; i<=data.size(); i++){
	    double pre=0.0;
	    for (int j=1; j<=i_para_scaling.size(); j++){
		pre += Fmat1(i,j)*Atb1(j,1);
	    }
	    tmpsum += log(tmp_bvec(i)/pre);
	}
	double beta = exp(tmpsum/data.size());

	for (int i=1; i<=i_para_scaling.size(); i++){
	    Atb1(i,1) *= beta/alpha1(i);
	    para(i_para_scaling(i))->setScaling(Atb1(i,1));
	}
    }else{
	// i_seq_scaling.size()>0 && i_para_scaling.size()==0 
	// but i_para_noscaling.size() may not be zero.
	//
	// edot = e1+e2+.. 1/(1/g1+1/g2+... +1/x1*h1+1/x2*h2+...)
	// 1/(edot-(e1+e2+...)) - (1/g1+1/g2+...) = 1/x1*h1+1/x2*h2+...
	//                                        = y1/h1+y2/h2+...
	// ---> tmp_bvec = Fmat*yvec
	
	// prepare tmp_bvec
	tmp_bvec = bvec; 
	if (do_bias){
	    for (int j=1; j<=data.size(); j++){
		tmp_bvec(j) *= exp(bias(uid[runid(j)])->value()*(-1.0));
	    }
	}
	for (int i=1; i<=i_para_noscaling.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		tmp_bvec(j) -= 
		    fmat_para(j,i_para_noscaling(i))*inv_dedot(j);
	    }
	}
	if (removeNonPositive(tmp_bvec)==false) return false;
	for (int j=1; j<=data.size(); j++){
	    tmp_bvec(j) = 1.0/tmp_bvec(j);
	}
	for (int i=1; i<=i_seq_noscaling.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		tmp_bvec(j) -= data(j).strainRateError()
		    /fmat_seq(j,i_seq_noscaling(i));
	    }
	}
	if (removeNonPositive(tmp_bvec)==false) return false;

	// prepare Fmat
	for (int i=1; i<=i_seq_scaling.size(); i++){
	    double val=0.0;
	    for (int j=1; j<=data.size(); j++){
		Fmat2(j,i) =
		    data(j).strainRateError()/fmat_seq(j,i_seq_scaling(i));
		val += log(Fmat2(j,i));
	    }
	    alpha2(i) = exp(val/data.size());
	}
	for (int i=1; i<=i_seq_scaling.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		Fmat2(j,i) /= alpha2(i);
	    }
	}

	// solve linear system
	for (int k=1; k<=i_seq_scaling.size(); k++){
	    for (int l=1; l<=k; l++){
		AtA2(k,l) = 0.0;
		for (int i=1; i<=data.size(); i++){
		    AtA2(k,l) += Fmat2(i,k)*Fmat2(i,l);
		}
		AtA2(l,k) = AtA2(k,l);
	    }
	    Atb2(k,1) = 0.0;
	    for (int i=1; i<=data.size(); i++){
		Atb2(k,1) += Fmat2(i,k)*tmp_bvec(i);
	    }
	}
	gaussj(AtA2.toRecipe(), i_seq_scaling.size(), Atb2.toRecipe(), 1);
	// check if there's a negative scaling constant...
	for (int i=1; i<=i_seq_scaling.size(); i++){
	    if (Atb2(i,1)<0) return false;
	}

	// post-scaling 
	double tmpsum=0.0;
	for (int i=1; i<=data.size(); i++){
	    double pre=0.0;
	    for (int j=1; j<=i_seq_scaling.size(); j++){
		pre += Fmat2(i,j)*Atb2(j,1);
	    }
	    tmpsum += log(tmp_bvec(i)/pre);
	}
	double beta = exp(tmpsum/data.size());
	for (int i=1; i<=i_seq_scaling.size(); i++){
	    double val = alpha2(i)/(beta*Atb2(i,1));
	    seq(i_seq_scaling(i))->setScaling(val);
	}
    }

    return true;
}

bool StatModel::removeNonPositive(Array1d<double>& v)
{	
    double vmin=1e30;
    for (int j=1; j<=v.size(); j++){
	if (v(j)>0 && v(j)<vmin){
	    vmin = v(j);
	}
    }

    if (vmin==1e30){ // i.e., all tmp_bvec is negative
	return false;
    }

    for (int j=1; j<=v.size(); j++){
	if (v(j)<0)
	    v(j) = vmin;
    }
    
    return true;
}

bool StatModel::calcBestFitA_CG()
{
    // set bounds on logA
    maxlogA = -100;
    minlogA = 100;
    int kk=1;
    for (int i=1; i<=i_para_scaling.size(); i++){
	for (int j=1; j<=data.size(); j++){
	    double diff = log(data(j).strainRate()) 
		-log(fmat_para(j,i_para_scaling(i)));
	    if (do_bias) diff -= bias(uid[runid(j)])->value();
	    if (diff>maxlogA(kk)) maxlogA(kk) = diff;
	    if (diff<minlogA(kk)) minlogA(kk) = diff;
	}
	maxlogA(kk) += 1.0; // extend the range a bit
	minlogA(kk) -= 1.0;
	kk++;
    }
    for (int i=1; i<=i_seq_scaling.size(); i++){
	for (int j=1; j<=data.size(); j++){
	    double diff = log(data(j).strainRate()) 
		-log(fmat_seq(j,i_seq_scaling(i)));
	    if (do_bias) diff -= bias(uid[runid(j)])->value();
	    if (diff>maxlogA(kk)) maxlogA(kk) = diff;
	    if (diff<minlogA(kk)) minlogA(kk) = diff;
	}
	maxlogA(kk) += 1.0; // extend the range a bit
	minlogA(kk) -= 1.0;
	kk++;
    }
    
    double min_prev_cost = 1e99;
    for (int itrial=1; itrial<=cg_ntrial; itrial++){
	// randomly initialize logA
	for (kk=1; kk<=logA.size(); kk++){
	    logA(kk) = minlogA(kk) + (maxlogA(kk)-minlogA(kk))*ran2(&cg_idum);
	}
	setA(logA);

	double prev_cost = CG_calc_cost_and_grad(true,grad);
	direc = -grad;
    
	// conjugate gradient search
	for (int iter=1; iter<=cg_iter_max; iter++){
	    // minimize along given gradient
	    double new_cost = CG_line_min();
	    if ((prev_cost-new_cost)<=cg_tol*data_lognorm2){
		break;
	    }
	
	    // calc conjugate gradient
	    CG_calc_cost_and_grad(true,new_grad);
	    double gg=0.0, dgg=0.0;
	    for (int i=1; i<=grad.size(); i++){
		gg += grad(i)*grad(i);
		dgg += new_grad(i)*new_grad(i);
	    }
	    if (abs(gg)<1e-10){
		break;
	    }else{
		dgg /= gg;
	    }
	    for (int i=1; i<=grad.size(); i++){
		direc(i) = -new_grad(i)+dgg*direc(i);
	    }
	    grad = new_grad;
	    prev_cost = new_cost;
	}

	for (int i=1; i<=logA.size(); i++){
	    // check for bounds
	    if ( (logA(i) > maxlogA(i)) || (logA(i) < minlogA(i)))
		prev_cost *= 10;
	}
	if (prev_cost<min_prev_cost){
	    min_prev_cost = prev_cost;
	    bestlogA = logA;
	}
    }

    setA(bestlogA);
    return true;
}

void StatModel::setA(const Array1d<double>& logA_)
{
    int kk=1;
    for (int i=1; i<=i_para_scaling.size(); i++){
	para(i_para_scaling(i))->setScaling(exp(logA_(kk)));
	kk++;
    }
    for (int i=1; i<=i_seq_scaling.size(); i++){
	seq(i_seq_scaling(i))->setScaling(exp(logA_(kk)));
	kk++;
    }
}

double StatModel::CG_calc_cost_and_grad(bool needGrad, 
					Array1d<double>& grad)
{
    double cost=0.0;
    if (needGrad) grad=0.0;

    for (int j=1; j<=data.size(); j++){
	double edot_total=0.0;
	double edot_total_seq=0.0;
	double edot_total_seq2=0.0;
	
	for (int i=1; i<=i_para_noscaling.size(); i++){
	    edot_total += fmat_para(j,i_para_noscaling(i));
	}
	for (int i=1; i<=i_para_scaling.size(); i++){
	    int ii = i_para_scaling(i); 
	    edot_total += para(ii)->predict(fmat_para(j,ii));
	}
	if (seq.size()>0){
	    double tmp=0.0;
	    for (int i=1; i<=i_seq_noscaling.size(); i++){
		tmp += 1.0/fmat_seq(j,i_seq_noscaling(i));
	    }
	    for (int i=1; i<=i_seq_scaling.size(); i++){
		int ii = i_seq_scaling(i); 
		tmp += 1.0/seq(ii)->predict(fmat_seq(j,ii));
	    }
	    edot_total_seq = 1.0/tmp;
	    edot_total_seq2 = edot_total_seq*edot_total_seq;
	    edot_total += edot_total_seq;
	}
	
	double misfit = log(data(j).strainRate())-log(edot_total);
	if (do_bias) misfit -= bias(uid[runid(j)])->value();
	cost += misfit*misfit;

	// calculate gradients
	if (needGrad){
	    int kk=1;
	    for (int i=1; i<=i_para_scaling.size(); i++){
		grad(kk) += misfit*fmat_para(j,i_para_scaling(i))/edot_total;
		kk++;
	    }
	    for (int i=1; i<=i_seq_scaling.size(); i++){
		int ii = i_seq_scaling(i); 
		double tmp = 1.0/seq(ii)->predict(fmat_seq(j,ii));
		grad(kk) += misfit*(fmat_seq(j,ii)/edot_total)
		    *(edot_total_seq2*tmp*tmp);
		kk++;
	    }
	}
    }
    
    // finish up gradients
    if (needGrad){
	int kk=1;
	for (int i=1; i<=i_para_scaling.size(); i++){
	    grad(kk) *= (-2.0)*para(i_para_scaling(i))->scaling();
	    kk++;
	}
	for (int i=1; i<=i_seq_scaling.size(); i++){
	    grad(kk) *= (-2.0)*seq(i_seq_scaling(i))->scaling();
	    kk++;
	}
    }

    return cost;
}

double StatModel::CG_line_min()
{
    double ax = 0.0;
    double xx = 1.0;
    double bx, fa, fx, fb, xmin;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, &StatModel::f1dim);
    double new_val = brent(ax,xx,bx,&xmin, &StatModel::f1dim);
    for (int i=1; i<=logA.size(); i++){
	logA(i) += xmin*direc(i);
    }
    setA(logA);

    return new_val;
}

double StatModel::f1dim(double x)
{
    for (int i=1; i<=logA.size(); i++){
	tmp_logA(i) = logA(i)+x*direc(i);
    }
    setA(tmp_logA);
    double tmp = CG_calc_cost_and_grad(false,grad);
    return tmp;
}

void StatModel::afterBestFitA()
{
    for (int i=1; i<=i_para_noscaling.size(); i++){
	int ii = i_para_noscaling(i);
	for (int j=1; j<=data.size(); j++){
	    edot_para(j,ii) = fmat_para(j,ii);
	}
    }
    for (int i=1; i<=i_para_scaling.size(); i++){
	int ii = i_para_scaling(i);
	for (int j=1; j<=data.size(); j++){
	    edot_para(j,ii) 
		= para(ii)->predict(fmat_para(j,ii));
	}
    }
    for (int i=1; i<=i_seq_noscaling.size(); i++){
	int ii = i_seq_noscaling(i);
	for (int j=1; j<=data.size(); j++){
	    edot_seq(j,ii) = fmat_seq(j,ii);
	}
    }
    for (int i=1; i<=i_seq_scaling.size(); i++){
	int ii = i_seq_scaling(i);
	for (int j=1; j<=data.size(); j++){
	    edot_seq(j,ii) 
		= seq(ii)->predict(fmat_seq(j,ii));
	}
    }
}

double StatModel::calcChiSq(int k, double& chi2_orig)
{
    fixParams();

    if (k<0 || !unfixed_params(k)->isScaling()){
	beforeBestFitA(k); 
	if (calcBestFitA() == false){
	    return max_chi2;
	}
    }
    afterBestFitA(); 
  
    double chi2=0.0;
    chi2_orig=0.0;
    for (int j=1; j<=data.size(); j++){
	double edot_total=0.0;
	double edot_total_seq=0.0;
	double edot_total_seq2=0.0;
	double var_sum=0.0;

	// calc edot_total (predicted) 
	for (int i=1; i<=para.size(); i++){
	    edot_total += edot_para(j,i);
	}
	if (seq.size()>0){
	    double tmp=0.0;
	    for (int i=1; i<=seq.size(); i++){
		tmp += 1.0/edot_seq(j,i);
	    }
	    edot_total_seq = 1.0/tmp;
	    edot_total_seq2 = edot_total_seq*edot_total_seq;
	    edot_total += edot_total_seq;
	}

	// calc total variance
	//	double val0 = data(j).strainRateError(); 
//	double val0 = data(j).strainRateError(); 
	//	double val0 = data(j).strainRateError(); 
	double val0 = data(j).strainRateRelError(); 
	double val0sq = val0*val0;
	var_sum += val0sq;
	for (int i=1; i<=para.size(); i++){
	    double val1 = edot_para(j,i)/edot_total;
	    var_sum += val1*val1*para(i)->variance(data(j));
	}
	for (int i=1; i<=seq.size(); i++){
	    double val2 = edot_total_seq2/(edot_total*edot_seq(j,i));
	    var_sum += val2*val2*seq(i)->variance(data(j));
	}

	// calc data misfit
//	if (do_bias) edot_total *= exp(bias(uid[runid(j)])->value());
//	double misfit = data(j).strainRate()-edot_total;
//	double misfit_sq = misfit*misfit;
	double misfit = log(data(j).strainRate())-log(edot_total);
	if (do_bias) misfit -= bias(uid[runid(j)])->value();
	double misfit_sq = misfit*misfit;

	chivec(j) = misfit_sq/val0sq;
	chi2 += chivec(j);
	chi2_orig += misfit_sq/var_sum;
    }

    return chi2;
}


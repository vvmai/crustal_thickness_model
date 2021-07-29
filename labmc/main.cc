/*
 * main.cc
 *
 * usage: labmc -D<data> -P<flowlaw_para> -S<flowlaw_seq>
 *              -M<niter>/<dn>/<m>/<dr> [ -F<i1/i2> -B<min>/<max> ]
 *              [ -R<seed> -C<n> -V ]
 *
 *        -D - sets input data file
 *             First 8 items should be (T dT p dp e de sig dsig)
 *             Last item should be run id (must be integer)
 *
 *        -d1  grain size (d dd)
 *        -d2  water content (COH dCOH)
 *        -d3  oxygen fugacity (fO2 dfO2)
 *        -d4  melt fraction (phi dphi)
 *
 *        -M - sets MCMC configuration
 *             <ninter> = maximum number of MC iterations
 *             <dn>     = output interval
 *             <m>      = number of trial calculations in rejection loop
 *             <dr>     = data randomization interval
 *
 *        -P - sets a parallel flow law
 *
 *             [dry diffusion creep]
 *             -Pa<facA_min>/<facA_max>/<mmin>/<mmax>
 *                /<Emin>/<Emax>/<Vmin>/<Vmax>
 *
 *             [dry dislocation creep]
 *             -Pb<facA_min>/<facA_max>/<nmin>/<nmax>
 *                /<Emin>/<Emax>/<Vmin>/<Vmax>
 *
 *             [wet diffusion creep]
 *             -Pc<facA_min>/<facA_max>/<mmin>/<mmax>
 *                 /<rmin>/<rmax>/<Emin>/<Emax>/<Vmin>/<Vmax>
 *
 *             [wet dislocation creep]
 *             -Pd<facA_min>/<facA_max>/<nmin>/<nmax>
 *                 /<rmin>/<rmax>/<Emin>/<Emax>/<Vmin>/<Vmax>
 *
 *             [gen flow creep]
 *             -Pe<facA_min>/<facA_max>/<nmin>/<nmax>
 *                 /<smin>/<smax>/<Qmin>/<Qmax>
 *
 *        -S - sets a sequential flow law
 *
 *             [gen flow creep]
 *             -Se<facA_min>/<facA_max>/<nmin>/<nmax>
 *                 /<smin>/<smax>/<Qmin>/<Qmax>
 *
 *        -F - fix param(i1) to param(i2)
 *             (i.e., two parameters always share the same value
 *              e.g., -F3/2 & -F4/2 --- param(3) and param(4) will be
 *                                      fixed to param(2)
 *              note: -F3/2 & -F4/3 probably won't give the correct result.)
 *
 *        -B - sets a priori bound for bias estimate (note: in term of log value)
 *        -R - sets seed for rand()
 *        -C - forces to use conjugate gradient search for `best-fit'
 *             scaling factors
 *             <n> = the number of different initial values to be tested
 *        -V - sets verbose mode
 *
 * Jun Korenaga
 * Summer 2008
 * Revised: Fall 2012 (for lambda-approach)
 * Revised: Fall 2014 (for data-randomization approach)
 */

#include <iostream>
#include <fstream>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "labmc_nr.h"
#include "constants.h"
#include "array.h"
#include "util.h"
#include "parameter.h"
#include "state.h"
#include "flowlaw.h"
#include "statmodel.h"

int main(int argc, char **argv)
{
    int nerror=0;
    char *dfn;
    bool getD=false, getM=false, useCG=false, verbose=false;
    int max_iter, dn_out, max_m, dn_ran;
    Array1d<int> ifix, isrc;
    int ran_seed=1;
    StatModel model;

    for (int i=1; i<argc; i++){
	if (argv[i][0] == '-'){
	    switch(argv[i][1]){
	    case 'D':
		dfn = &argv[i][2];
		getD = true;
		break;
	    case 'd':
	    {
		int istate = atoi(&argv[i][2]);
		switch(istate){
		case 1:
		    model.stateToRead(State::grain_size);
		    break;
		case 2:
		    model.stateToRead(State::water_content);
		    break;
		case 3:
		    model.stateToRead(State::oxygen_fugacity);
		    break;
		case 4:
		    model.stateToRead(State::melt_fraction);
		    break;
		}
		break;
	    }
	    case 'M':
		if (sscanf(&argv[i][2],
			   "%d/%d/%d/%d", &max_iter, &dn_out, 
			   &max_m, &dn_ran) != 4){
		    cerr << "invalid -M option\n";
		    nerror++;
		}
		getM = true;
		break;
	    case 'P':
		switch(argv[i][2]){
		case 'a':
		{
		    double a1,a2,m1,m2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&m1,&m2,&e1,&e2,&v1,&v2) != 8){
			cerr << "invalid -Pa option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawDiffDry(a1,a2,m1,m2,e1,e2,v1,v2));
		    break;
		}

		case 'b':
		{
		    double a1,a2,n1,n2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&e1,&e2,&v1,&v2) != 8){
			cerr << "invalid -Pb option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawDisDry(a1,a2,n1,n2,e1,e2,v1,v2));
		    break;
		}
		case 'c':
		{
		    double a1,a2,m1,m2,r1,r2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&m1,&m2,&r1,&r2,&e1,&e2,&v1,&v2) != 10){
			cerr << "invalid -Pc option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawDiffWet(a1,a2,m1,m2,r1,r2,
						     e1,e2,v1,v2));
		    break;
		}
		case 'd':
		{
		    double a1,a2,r1,r2,n1,n2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&r1,&r2,&e1,&e2,&v1,&v2) != 10){
			cerr << "invalid -Pd option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawDisWet(a1,a2,n1,n2,r1,r2,
						     e1,e2,v1,v2));
		    break;
		}
		case 'e':
		{
		  double a1,a2,s1,s2,n1,n2,e1,e2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&s1,&s2,&e1,&e2) != 8){
			cerr << "invalid -Pe option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawGen(a1,a2,n1,n2,s1,s2,
						     e1,e2));
		    break;
		}
		// Vuong add case f (disGBS)
		case 'f':
		{
		    double a1,a2,m1,m2,n1,n2,r1,r2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&m1,&m2,&n1,&n2,&r1,&r2,&e1,&e2,&v1,&v2) != 12){
			cerr << "invalid -Pf option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawdisGBS(a1,a2,m1,m2,n1,n2,r1,r2,
						     e1,e2,v1,v2));
		    break;
		}
		// end vuong

		default:
		    cerr << "unknown -P option\n";
		    nerror++;
		    break;
		}
		break;
	    case 'S':
	    switch(argv[i][2]){
	    	case 'd':
		{
		    double a1,a2,r1,r2,n1,n2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&r1,&r2,&e1,&e2,&v1,&v2) != 10){
			cerr << "invalid -Sd option\n";
			nerror++;
		    }
		    model.addSequential(new
					FlowLawDisWet(a1,a2,n1,n2,r1,r2,
						      e1,e2,v1,v2));
		    break;
		}
	    	case 'e':
	    	 {
	    	 double a1,a2,s1,s2,n1,n2,e1,e2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&s1,&s2,&e1,&e2) != 8){
			cerr << "invalid -Se option\n";
			nerror++;
		    }
		    model.addSequential(new
					FlowLawGen(a1,a2,n1,n2,s1,s2,
						   e1,e2));
		    break;
		}


		default:
		    cerr << "unknown -S option\n";
		    nerror++;
		    break;
		}
		break;

	    case 'F':
	    {
		int i1, i2;
		if (sscanf(&argv[i][2],
			   "%d/%d", &i1, &i2) != 2){
		    cerr << "invalid -F option\n";
		    nerror++;
		}
		ifix.push_back(i1);
		isrc.push_back(i2);
		break;
	    }
	    case 'B':
	    {
		double min_bias, max_bias;
		if (sscanf(&argv[i][2],
			   "%lf/%lf", &min_bias, &max_bias) != 2){
		    cerr << "invalid -B option\n";
		    nerror++;
		}
		model.setupBiasCorrection(min_bias, max_bias);
		break;
	    }
	    case 'R':
		ran_seed = atoi(&argv[i][2]);
		Parameter::setRanSeed(ran_seed);
		break;
	    case 'C':
	    {
		useCG = true;
		int ntrial = atoi(&argv[i][2]);
		model.useConjugateGradient(ntrial);
		break;
	    }
	    case 'V':
		verbose = true;
 		break;
	    default:
		cerr << "unknown options\n";
		nerror++;
		break;
	    }
	}else{
	    cerr << "invalid syntax\n";
	    nerror++;
	}
    }

    if (!getD || !getM) nerror++;
    if (model.numParallel()+model.numSequential()==0){

	cerr << "no flow law specified\n";
	nerror++;
    }
    if (nerror){
	cerr << "invalid command option[s] - abort\n";
	exit(1);
    }
    if (verbose){
	cerr << "the number of parallel flow laws: "
	     << model.numParallel() << '\n';
	cerr << "the number of sequential flow laws: "
	     << model.numSequential() << '\n';
    }

    //
    // read experimental data
    //
    model.readData(dfn);
    int ndata = model.numData();
    if (verbose){
	cerr << "the number of data: " << ndata << '\n';
	model.printRunIds(cerr);
    }

    //
    // run MCMC with Gibbs sampling
    //

    long idum = long(-(abs(ran_seed)+1)); // for ran2()
    model.setUp();

    int np = model.numUnfixedParams();
    Array1d<int> ifixed(np);
    ifixed=0;

    for (int i=1; i<=ifix.size(); i++){
	if (ifix(i)>0 && ifix(i)<=np
	    && isrc(i)>0 && isrc(i)<=np){
	    ifixed(ifix(i)) = isrc(i);
	    if (verbose){
		cerr << model.unfixedParam(ifix(i))->name()
		     << " will be fixed to "
		     << model.unfixedParam(isrc(i))->name()
		     << "\n";
	    }
	}else{
	    error("invalid -F option");
	}
    }

    model.addConstraints(ifixed);
    
    Array1d<double> sampled_mval(max_m);

    for (int ii=1; ii<=max_iter;ii++){
	double chi2, chi2_orig;
	int m_chi2min;

	//
	// randomize data
	//
	if (ii%dn_ran==0) model.randomizeData();

	//
	// single scan for ordinary model parameters
	//
	int k = int(round(ran2(&idum)*(np-1)+1));
	if (ifixed(k)>0) k = ifixed(k);

	// estimate conditional probability distribution
	double chi2min = 1e30;
	for (int m=1; m<=max_m; m++){
	    model.unfixedParam(k)->randomize();

	    chi2 = model.calcChiSq(k,chi2_orig);
	    if (chi2<chi2min){
		chi2min = chi2;
		m_chi2min = m;
	    }

	    sampled_mval(m) = model.unfixedParam(k)->value();
	}

	// now with the rejection method
	int ir=0;
	while (true){
	    model.unfixedParam(k)->randomize();
	    chi2 = model.calcChiSq(k,chi2_orig);
	    ir++;
	    
	    double prob = exp(-0.5*(chi2-chi2min));
	    if (ran2(&idum)<prob){
		break;
	    }
	    if (ir>max_m){ // this extra if-clause is for efficiency
		model.unfixedParam(k)->setValue(sampled_mval(m_chi2min));
		chi2 = model.calcChiSq(k,chi2_orig);
		// note: calling calcChiSq() is important because scaling factors
		// are calculated here. 
		break;
	    }
	}

	// output the current solution
	if (ii%dn_out == 0){
	    char line[MaxStr];
 	    sprintf(line, "%6d %5.3f %5.3f ", ii, chi2/ndata, chi2_orig/ndata);
 	    cout << line;

	    for (int i=1; i<=model.numUnfixedParams(); i++){
		sprintf(line, "%6.5e ", model.unfixedParam(i)->value());
		cout << line;
	    }

	    for (int i=1; i<=model.numParallel(); i++){
		if (model.parallel(i)->needScaling()){
		    sprintf(line, "%6.5e ", model.parallel(i)->scaling());
		    cout << line;
		}

	    }

	    for (int i=1; i<=model.numSequential(); i++){
		if (model.sequential(i)->needScaling()){
		    sprintf(line, "%6.5e ", model.sequential(i)->scaling());
		    cout << line;
		}
	    }
	    cout << '\n';
	    cout.flush();
	}
    }
}


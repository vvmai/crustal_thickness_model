/*
 * flowlaw.h
 *
 * various flow-law-related classes
 *
 * Jun Korenaga, Summer 2008
 */

#ifndef _LABMC_FLOWLAW_H_
#define _LABMC_FLOWLAW_H_

#include "array.h"
#include "parameter.h"
#include "state.h"
#include "labmc_nr.h"

// individual flow law base class
class FlowLaw {

public:
    virtual double prepPredict(const State& s) const = 0;
    virtual bool needScaling() const = 0;
    virtual void setScaling(double) = 0;
    virtual double scaling() const = 0;
    virtual double predict(const State& s) const = 0;
    virtual double predict(double prep) const = 0;
    virtual double variance(const State& s) const = 0;

    int numUnfixedParams() const { return unfixed_params.size(); }
    Parameter* unfixedParam(int i) { return unfixed_params(i); }
    bool checkStates(const Array1d<int>& states) const;

protected:
    Array1d<Parameter*> unfixed_params;
    Array1d<int> required_states;

};

class FlowLawGen : public FlowLaw {
public:
	//edot = A*sigma^n*fO2^s*exp(-E/(RT))
	//
	//state variables: sigma, fO2, T
	//parameters: A, n, s, E

	FlowLawGen(double facA_min, double facA_max,
		   double n_min, double n_max,
		   double s_min, double s_max,
		   double E_min, double E_max);
	~FlowLawGen();
	double prepPredict(const State& s) const;
    	bool needScaling() const { return true; }
    	void setScaling(double val) { A = val; }
    	double scaling() const { return A; }
    	double predict(const State& s) const;
    	double predict(double prep) const;
    	double variance(const State& s) const;

private:

    	double A;
    	Parameter* facA__;
    	Parameter* E__; // activation energy (J/mol)
    	Parameter* n__; // stress exponent
    	Parameter* s__; // fO2 exponent
};

class FlowLawDiffDry : public FlowLaw {
public:
    // edot = A*d^(-m)*sigma*exp(-(E+pV)/(RT))
    //
    // state variables: d, sigma, p, T
    // parameters: A, m, E, V
    //
    // note: Units for E and V in the constructor are kJ/mol and cm^3/mol,
    // respectively.
    FlowLawDiffDry(double facA_min, double facA_max,
		   double m_min, double m_max,
		   double E_min, double E_max,
		   double V_min, double V_max);
    ~FlowLawDiffDry();

    double prepPredict(const State& s) const;
    bool needScaling() const { return true; }
    void setScaling(double val) { A = val; }
    double scaling() const { return A; }
    double predict(const State& s) const;
    double predict(double prep) const;
    double variance(const State& s) const;

private:

    double A;

    Parameter* facA__;

    Parameter* E__; // activation energy (J/mol)
    Parameter* V__; // activation volume (m^3/mol)

    Parameter* m__; // grain-size exponent
};



class FlowLawDisDry : public FlowLaw {
public:
    // edot = A*sigma^n*exp(-(E+pV)/(RT))
    //
    // state variables: sigma, p, T
    // parameters: A, n, E, V
    //
    // note: Units for E and V in the constructor are kJ/mol and cm^3/mol,
    // respectively.
    FlowLawDisDry(double facA_min, double facA_max,
		  double n_min, double n_max,
		  double E_min, double E_max,
		  double V_min, double V_max);
    ~FlowLawDisDry();

    double prepPredict(const State& s) const;
    bool needScaling() const { return true; }
    void setScaling(double val) { A = val; }
    double scaling() const { return A; }
    double predict(const State& s) const;
    double predict(double prep) const;
    double variance(const State& s) const;

private:
    double A;

    Parameter* facA__;
    Parameter* n__; // stress exponent
    Parameter* E__; // activation energy (J/mol)
    Parameter* V__; // activation volume (m^3/mol)
};

class FlowLawDiffWet : public FlowLaw {
public:
    // edot = A*d^(-m)*COH^r*sigma*exp(-(E+pV)/(RT))
    //
    // state variables: d, COH, sigma, p, T
    // parameters: A, m, r, E, V
    //
    // note: Units for E and V in the constructor are kJ/mol and cm^3/mol,
    // respectively.
    FlowLawDiffWet(double facA_min, double facA_max,
		   double m_min, double m_max,
		   double r_min, double r_max,
		   double E_min, double E_max,
		   double V_min, double V_max);
    ~FlowLawDiffWet();

    double prepPredict(const State& s) const;
    bool needScaling() const { return true; }
    void setScaling(double val) { A = val; }
    double scaling() const { return A; }
    double predict(const State& s) const;
    double predict(double prep) const;
    double variance(const State& s) const;

private:
    double A;

    Parameter* facA__;
    Parameter* m__; // grain-size exponent
    Parameter* r__; // water-content exponent
    Parameter* E__; // activation energy (J/mol)
    Parameter* V__; // activation volume (m^3/mol)
};

class FlowLawDisWet : public FlowLaw {
public:
    // edot = A*COH^r*sigma^n*exp(-(E+pV)/(RT))
    //
    // state variables: COH, sigma, p, T
    // parameters: A, r, n, E, V
    //
    // note: Units for E and V in the constructor are kJ/mol and cm^3/mol,
    // respectively.
    FlowLawDisWet(double facA_min, double facA_max,
		  double n_min, double n_max,
		  double r_min, double r_max,
		  double E_min, double E_max,
		  double V_min, double V_max);
    ~FlowLawDisWet();

    double prepPredict(const State& s) const;
    bool needScaling() const { return true; }
    void setScaling(double val) { A = val; }
    double scaling() const { return A; }
    double predict(const State& s) const;
    double predict(double prep) const;
    double variance(const State& s) const;

private:
    double A;

    Parameter* facA__;
    Parameter* r__; // water-content exponent
    Parameter* n__; // stress exponent
    Parameter* E__; // activation energy (J/mol)
    Parameter* V__; // activation volume (m^3/mol)
};

// Vuong flow law disGBS18
class FlowLawdisGBS : public FlowLaw {
public:
    // edot = A*d^(-m)*COH^r*sigma^n*exp(-(E+pV)/(RT))
    //
    // state variables: d, COH, sigma, p, T
    // parameters: A, m, n, r, E, V
    //
    // note: Units for E and V in the constructor are kJ/mol and cm^3/mol,
    // respectively.
    FlowLawdisGBS(double facA_min, double facA_max,
		   double m_min, double m_max,
           double n_min, double n_max,
		   double r_min, double r_max,
		   double E_min, double E_max,
		   double V_min, double V_max);
    ~FlowLawdisGBS();

    double prepPredict(const State& s) const;
    bool needScaling() const { return true; }
    void setScaling(double val) { A = val; }
    double scaling() const { return A; }
    double predict(const State& s) const;
    double predict(double prep) const;
    double variance(const State& s) const;

private:
    double A;

    Parameter* facA__;
    Parameter* m__; // grain-size exponent
    Parameter* n__; // stress exponent
    Parameter* r__; // water-content exponent
    Parameter* E__; // activation energy (J/mol)
    Parameter* V__; // activation volume (m^3/mol)
};

// end vuong
#endif /* _LABMC_FLOWLAW_H_ */

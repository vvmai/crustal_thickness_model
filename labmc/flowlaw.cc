/*
 * flowlaw.cc
 *
 * Jun Korenaga
 * Summer 2008
 */

#include <cmath>
#include "constants.h"
#include "flowlaw.h"
#include "labmc_nr.h"

// FlowLaw
bool FlowLaw::checkStates(const Array1d<int>& states) const
{
    int ncount=0;
    for (int i=1; i<=required_states.size(); i++){
	int ii=required_states(i);
	for (int j=1; j<=states.size(); j++){
	    if (states(j) == ii){
		ncount++;
		continue;
	    }
	}
    }

    return (ncount == required_states.size()) ? true : false;
}

//FlowLawGen
FlowLawGen::FlowLawGen(	double facA_min, double facA_max,
		   	double n_min, double n_max,
		   	double s_min, double s_max,
		   	double E_min, double E_max)
{
    required_states.resize(1);
    required_states(1) = State::oxygen_fugacity;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"Gen:facA",true);
	unfixed_params.push_back(facA__);
    }else if (facA_min == facA_max){
	facA__ = new Parameter(facA_min,"Gen:facA",true);
    }else{
	nerror++;
    }

    if (n_min < n_max){
	n__ = new Parameter(n_min,n_max,"Gen:n");
	unfixed_params.push_back(n__);
    }else if (n_min == n_max){
	n__ = new Parameter(n_min,"Gen:n");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"Gen:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"Gen:E");
    }else{
	nerror++;
    }

    if (s_min < s_max){
	s__ = new Parameter(s_min,s_max,"Gen:s");
	unfixed_params.push_back(s__);
    }else if (s_min == s_max){
	s__ = new Parameter(s_min,"Gen:s");
    }else{
	nerror++;
    }

    A = 1.0;
    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawGen::~FlowLawGen()
{
    delete facA__;
    delete n__;
    delete E__;
    delete s__;
}

double FlowLawGen::prepPredict(const State& s) const
{
    double val1 = pow(s.oxygenFugacity(), s__->value());
    double val2 = pow(s.stress(), n__->value());
    double val3 = E__->value();
    double val4 = Rgas*s.temperature();

    return val1*val2*exp(-val3/val4);
}

double FlowLawGen::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawGen::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawGen::variance(const State& s) const
{
    double val1 = (s__->value())*s.oxygenFugacityRelError();
    double val2 = (n__->value())*s.stressRelError();
    double T = s.temperature();

    double val3
	= ((E__->value())/(Rgas*T*T))
	*s.temperatureError();


    double var = val1*val1+val2*val2+val3*val3;
    return var;
}


// FlowLawDiffDry
FlowLawDiffDry::FlowLawDiffDry(double facA_min, double facA_max,
			       double m_min, double m_max,
			       double E_min, double E_max,
			       double V_min, double V_max)
{
    required_states.resize(1);
    required_states(1) = State::grain_size;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"DiffDry:facA",true);
	unfixed_params.push_back(facA__);
    }else if (facA_min == facA_max){
	facA__ = new Parameter(facA_min,"DiffDry:facA",true);
    }else{
	nerror++;
    }

    if (m_min < m_max){
	m__ = new Parameter(m_min,m_max,"DiffDry:m");
	unfixed_params.push_back(m__);
    }else if (m_min == m_max){
	m__ = new Parameter(m_min,"DiffDry:m");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"DiffDry:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"DiffDry:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"DiffDry:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"DiffDry:V");
    }else{
	nerror++;
    }

    A = 1.0;
    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawDiffDry::~FlowLawDiffDry()
{
    delete facA__;
    delete m__;
    delete E__;
    delete V__;
}

double FlowLawDiffDry::prepPredict(const State& s) const
{
    double val1 = pow(s.grainSize(), -1.0*(m__->value()));
    double val2 = E__->value()+s.pressure()*V__->value();
    double val3 = Rgas*s.temperature();

    return val1*s.stress()*exp(-val2/val3);
}

double FlowLawDiffDry::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawDiffDry::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawDiffDry::variance(const State& s) const
{
    double val1 = (m__->value())*s.grainSizeRelError();
    double val2 = s.stressRelError();
    double T = s.temperature();
    double val3 = (V__->value()/(Rgas*T))*s.pressureError();
    double val4
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4;
    return var;
}

// FlowLawDisDry
FlowLawDisDry::FlowLawDisDry(double facA_min, double facA_max,
			     double n_min, double n_max,
			     double E_min, double E_max,
			     double V_min, double V_max)
{
    required_states.resize(0);

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"DisDry:facA",true);
	unfixed_params.push_back(facA__);
    }else if (facA_min == facA_max){
	facA__ = new Parameter(facA_min,"DisDry:facA",true);
    }else{
	nerror++;
    }

    if (n_min < n_max){
	n__ = new Parameter(n_min,n_max,"DisDry:n");
	unfixed_params.push_back(n__);
    }else if (n_min == n_max){
	n__ = new Parameter(n_min,"DisDry:n");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"DisDry:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"DisDry:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"DisDry:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"DisDry:V");
    }else{
	nerror++;
    }

    A = 1.0;
    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawDisDry::~FlowLawDisDry()
{
    delete facA__;
    delete n__;
    delete E__;
    delete V__;
}

double FlowLawDisDry::prepPredict(const State& s) const
{
    double val1 = pow(s.stress(), n__->value());
    double val2 = E__->value()+s.pressure()*V__->value();
    double val3 = Rgas*s.temperature();

    return val1*exp(-val2/val3);
}

double FlowLawDisDry::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawDisDry::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawDisDry::variance(const State& s) const
{
    double val1 = (n__->value())*s.stressRelError();
    double T = s.temperature();
    double val2 = (V__->value()/(Rgas*T))*s.pressureError();
    double val3
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3;
    return var;
}

// FlowLawDiffWet
FlowLawDiffWet::FlowLawDiffWet(double facA_min, double facA_max,
			       double m_min, double m_max,
			       double r_min, double r_max,
			       double E_min, double E_max,
			       double V_min, double V_max)
{
    required_states.resize(2);
    required_states(1) = State::grain_size;
    required_states(2) = State::water_content;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"DiffWet:facA",true);
	unfixed_params.push_back(facA__);
    }else if (facA_min == facA_max){
	facA__ = new Parameter(facA_min,"DiffWet:facA",true);
    }else{
	nerror++;
    }

    if (m_min < m_max){
	m__ = new Parameter(m_min,m_max,"DiffWet:m");
	unfixed_params.push_back(m__);
    }else if (m_min == m_max){
	m__ = new Parameter(m_min,"DiffWet:m");
    }else{
	nerror++;
    }

    if (r_min < r_max){
	r__ = new Parameter(r_min,r_max,"DiffWet:r");
	unfixed_params.push_back(r__);
    }else if (r_min == r_max){
	r__ = new Parameter(r_min,"DiffWet:r");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"DiffWet:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"DiffWet:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"DiffWet:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"DiffWet:V");
    }else{
	nerror++;
    }

    A = 1.0;
    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawDiffWet::~FlowLawDiffWet()
{
    delete facA__;
    delete m__;
    delete r__;
    delete E__;
    delete V__;
}

double FlowLawDiffWet::prepPredict(const State& s) const
{
    double val1 = pow(s.grainSize(), -1.0*(m__->value()));
    double val2 = pow(s.waterContent(), r__->value());
    double val3 = E__->value()+s.pressure()*V__->value();
    double val4 = Rgas*s.temperature();

    return val1*val2*s.stress()*exp(-val3/val4);
}

double FlowLawDiffWet::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawDiffWet::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawDiffWet::variance(const State& s) const
{
    double val1 = (m__->value())*s.grainSizeRelError();
    double val2 = (r__->value())*s.waterContentRelError();
    double val3 = s.stressRelError();
    double T = s.temperature();
    double val4 = (V__->value()/(Rgas*T))*s.pressureError();
    double val5
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4+val5*val5;
    return var;
}

// FlowLawDisWet
FlowLawDisWet::FlowLawDisWet(double facA_min, double facA_max,
			     double n_min, double n_max,
			     double r_min, double r_max,
			     double E_min, double E_max,
			     double V_min, double V_max)
{
    required_states.resize(1);
    required_states(1) = State::water_content;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"DisWet:facA",true);
	unfixed_params.push_back(facA__);
    }else if (facA_min == facA_max){
	facA__ = new Parameter(facA_min,"DisWet:facA",true);
    }else{
	nerror++;
    }

    if (r_min < r_max){
	r__ = new Parameter(r_min,r_max,"DisWet:r");
	unfixed_params.push_back(r__);
    }else if (r_min == r_max){
	r__ = new Parameter(r_min,"DisWet:r");
    }else{
	nerror++;
    }

    if (n_min < n_max){
	n__ = new Parameter(n_min,n_max,"DisWet:n");
	unfixed_params.push_back(n__);
    }else if (n_min == n_max){
	n__ = new Parameter(n_min,"DisWet:n");
    }else{
 	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"DisWet:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"DisWet:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"DisWet:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"DisWet:V");
    }else{
	nerror++;
    }

    A = 1.0;
    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawDisWet::~FlowLawDisWet()
{
    delete facA__;
    delete n__;
    delete E__;
    delete V__;
}

double FlowLawDisWet::prepPredict(const State& s) const
{
    double val1 = pow(s.waterContent(), r__->value());
    double val2 = pow(s.stress(), n__->value());
    double val3 = E__->value()+s.pressure()*V__->value();
    double val4 = Rgas*s.temperature();

    return val1*val2*exp(-val3/val4);
}

double FlowLawDisWet::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawDisWet::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawDisWet::variance(const State& s) const
{
    double val1 = (r__->value())*s.waterContentRelError();
    double val2 = (n__->value())*s.stressRelError();
    double T = s.temperature();
    double val3 = (V__->value()/(Rgas*T))*s.pressureError();
    double val4
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4;
    return var;
}

// Vuong flow law disGBS18
FlowLawdisGBS::FlowLawdisGBS(double facA_min, double facA_max,
			       double m_min, double m_max,
                   double n_min, double n_max,
			       double r_min, double r_max,
			       double E_min, double E_max,
			       double V_min, double V_max)
{
    required_states.resize(2);
    required_states(1) = State::grain_size;
    required_states(2) = State::water_content;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"disGBS:facA",true);
	unfixed_params.push_back(facA__);
    }else if (facA_min == facA_max){
	facA__ = new Parameter(facA_min,"disGBS:facA",true);
    }else{
	nerror++;
    }

    if (m_min < m_max){
	m__ = new Parameter(m_min,m_max,"disGBS:m");
	unfixed_params.push_back(m__);
    }else if (m_min == m_max){
	m__ = new Parameter(m_min,"disGBS:m");
    }else{
	nerror++;
    }

    if (n_min < n_max){
	n__ = new Parameter(n_min,n_max,"disGBS:n");
	unfixed_params.push_back(n__);
    }else if (n_min == n_max){
	n__ = new Parameter(n_min,"disGBS:n");
    }else{
 	nerror++;
    }

    if (r_min < r_max){
	r__ = new Parameter(r_min,r_max,"disGBS:r");
	unfixed_params.push_back(r__);
    }else if (r_min == r_max){
	r__ = new Parameter(r_min,"disGBS:r");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"disGBS:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"disGBS:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"disGBS:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"disGBS:V");
    }else{
	nerror++;
    }

    A = 1.0;
    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawdisGBS::~FlowLawdisGBS()
{
    delete facA__;
    delete m__;
    delete n__;
    delete r__;
    delete E__;
    delete V__;
}

double FlowLawdisGBS::prepPredict(const State& s) const
{
    double val1 = pow(s.grainSize(), -1.0*(m__->value()));
    double val2 = pow(s.stress(), n__->value());
    double val3 = pow(s.waterContent(), r__->value());
    double val4 = E__->value()+s.pressure()*V__->value();
    double val5 = Rgas*s.temperature();

    return val1*val2*val3*exp(-val4/val5);
}

double FlowLawdisGBS::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawdisGBS::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawdisGBS::variance(const State& s) const
{
    double val1 = (m__->value())*s.grainSizeRelError();
    double val2 = (n__->value())*s.stressRelError();
    double val3 = (r__->value())*s.waterContentRelError();
    double val4 = s.stressRelError();
    double T = s.temperature();
    double val5 = (V__->value()/(Rgas*T))*s.pressureError();
    double val6
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4+val5*val5+val6*val6;
    return var;
}
// end vuong
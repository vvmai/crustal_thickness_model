/*
 * state.cc
 *
 * Jun Korenaga
 * Summer 2008
 */

#include "error.h"
#include "state.h"

State::State()
{}

void State::setStrainRate(double val, double dval)
{
    edot = val; 
    dedot = dval;
    if (edot<0 || dedot<0) error("invalid strain rate detected.");

    r_dedot = dedot/edot;
}

void State::setStress(double val, double dval)
{
    sigma = val; 
    dsigma = dval;
    if (sigma<0 || dsigma<0) error("invalid stress detected.");

    r_dsigma = dsigma/sigma;
}

void State::setGrainSize(double val, double dval)
{
    d = val; 
    dd = dval;
    if (d<0 || dd<0) error("invalid grain size detected.");

    r_dd = (d == 0.0 ? 0.0 : dd/d); 
}

void State::setTemperature(double val, double dval)
{
    T = val; 
    dT = dval;
    if (T<=0 || dT<0) error("invalid temperature detected.");
}

void State::setPressure(double val, double dval)
{
    p = val; 
    dp = dval;
    if (p<0 || dp<0) error("invalid presure detected.");
}

void State::setWaterContent(double val, double dval)
{
    COH = val; 
    dCOH = dval;
    if (COH<0 || dCOH<0) error("invalid water content detected.");

    r_dCOH = (COH == 0.0 ? 0.0 : dCOH/COH);
}

void State::setOxygenFugacity(double val, double dval)
{
    fO2 = val; 
    dfO2 = dval;
    if (fO2<=0.0 || dfO2<0.0) error("invalid oxygen fugacity detected.");

    r_dfO2 = dfO2/fO2;
}

void State::setMeltFraction(double val, double dval)
{
    phi = val; 
    dphi = dval;
    if (phi<0 || dphi<0) error("invalid melt fraction detected.");
}

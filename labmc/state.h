/*
 * state.h
 *
 * Jun Korenaga
 * Summer 2008
 */

#ifndef _LABMC_STATE_H_
#define _LABMC_STATE_H_

class State {
public:
    State(); 
    ~State(){}

    void setStrainRate(double, double);
    void setStress(double, double);
    void setGrainSize(double, double);
    void setTemperature(double, double);
    void setPressure(double, double);
    void setWaterContent(double, double);
    void setOxygenFugacity(double, double);
    void setMeltFraction(double, double);

    double strainRate() const { return edot; }
    double strainRateError() const { return dedot; }
    double stress() const { return sigma; }
    double stressError() const { return dsigma; }
    double grainSize() const { return d; }
    double grainSizeError() const { return dd; }
    double temperature() const { return T; }
    double temperatureError() const { return dT; }
    double pressure() const { return p; }
    double pressureError() const { return dp; }
    double waterContent() const { return COH; }
    double waterContentError() const { return dCOH; }
    double oxygenFugacity() const { return fO2; }
    double oxygenFugacityError() const { return dfO2; }
    double meltFraction() const { return phi; }
    double meltFractionError() const { return dphi; }

    double strainRateRelError() const { return r_dedot; }
    double grainSizeRelError() const { return r_dd; }
    double stressRelError() const { return r_dsigma; }
    double waterContentRelError() const { return r_dCOH; }
    double oxygenFugacityRelError() const { return r_dfO2; }

    enum {grain_size, water_content, oxygen_fugacity, melt_fraction};

private:
    double edot, dedot; // strain rate (2nd invariant, in 1/s)
    double sigma, dsigma; // stress (2nd invariant, in MPa)
    double d, dd; // grain size (in micron)
    double T, dT; // temperature (in K)
    double p, dp; // pressure (in Pa)
    double COH, dCOH; // water content (in ppm H/Si)
    double fO2, dfO2; // oxygen fugacity (in atm)
    double phi, dphi; // melt fraction

    double r_dedot, r_dd, r_dsigma, r_dCOH, r_dfO2;
};

#endif /* _LABMC_STATE_H_ */


/*
 * parameter.h
 *
 * Jun Korenaga
 * Summer 2008
 */

#ifndef _LABMC_PARAMETER_H_
#define _LABMC_PARAMETER_H_

#include <string>
#include <cstdlib>
#include "labmc_nr.h"

using namespace std;

// parameter class
class Parameter {
public:
    Parameter(double val, const char* str, bool scaling=false);
    Parameter(double min, double max, const char* str, bool scaling=false);

    double value() const { return val_; }
    void setValue(double x) { val_ = x; }
    string name() const { return name_; }
    bool isFixed() const { return is_fixed; }
    bool isScaling() const { return is_scaling; }
    void randomize();

    static void setRanSeed(long seed) { idum = long(-(abs(seed)+1)); }

private:
    double val_;
    double min_;
    double range_;
    bool is_fixed;
    bool is_scaling;
    string name_;
    static long idum; // for ran2()
};

#endif /* _LABMC_PARAMETER_H_ */

/*
 * parameter.cc
 *
 * Jun Korenaga
 * Summer 2008
 */

#include "labmc_nr.h"
#include "parameter.h"

long Parameter::idum = long(-(abs(1)+1));

Parameter::Parameter(double val, const char* str, bool scaling)
{
    val_ = val;
    min_ = val;
    range_ = 0.0;
    name_ = str;
    is_fixed = true;
    is_scaling = scaling;
}

Parameter::Parameter(double min, double max, const char* str, bool scaling)
{
    min_ = min;
    range_ = max-min;
    name_ = str;
    is_fixed = false;
    is_scaling = scaling;
}

void Parameter::randomize()
{
    if (!is_fixed){
	val_ = min_+ran2(&idum)*range_;
    }
}


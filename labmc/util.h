/*
 * util.h
 */

#ifndef _JK_UTIL_H_
#define _JK_UTIL_H_

#include "array.h"
#include <cstdio>
#include <cstdlib>

const int MaxStr = 512;

FILE *openFile(char *fn, char *mode);
int countLines(char *fn);
int set_ivals(char*, Array1d<int>&);
int set_dvals(char*, Array1d<double>&);
long int make_seed();

/* the following routines are from Numerical Recipes */
#define NR_END 1
#define FREE_ARG char*
void gaussj(double **a, int n, double **b, int m);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh); 


#endif /* _JK_UTIL_H_ */

/*
 * util.cc
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <ctime>
#include "util.h"
#include "error.h"

FILE *openFile(char *fn, char *mode)
{
    FILE *tmp;
    
    if ((tmp=fopen(fn, mode)) == NULL){
	fprintf(stderr, "[openFile] can't open file %s\n", fn);
	exit(1);
    }
    return tmp;
}

int countLines(char *fn)
{
    FILE *tmpfp;
    int nline=0;
    char line[MaxStr];

    tmpfp = openFile(fn, (char*) "r");
    while(fgets(line, MaxStr, tmpfp) != NULL){
	nline++;
    }
    fclose(tmpfp);

    return nline;
}

int set_ivals(char* line, Array1d<int>& pp)
{
    const char* SEPCHARS = " ";
    
    pp.resize(0);
    char *p = strtok(line, SEPCHARS);
    while (p != NULL){
	int tmpi;
	sscanf(p, "%d", &tmpi);
	pp.push_back(tmpi);
	p = strtok(NULL, SEPCHARS);
    }
    return pp.size();
}

int set_dvals(char* line, Array1d<double>& pp)
{
    const char* SEPCHARS = " ";

    pp.resize(0);
    char *p = strtok(line, SEPCHARS);
    while (p != NULL){
	double tmpi;
	sscanf(p, "%lf", &tmpi);
	pp.push_back(tmpi);
	p = strtok(NULL, SEPCHARS);
    }
    return pp.size();
}

long int make_seed()
{
    time_t now = time(NULL);
    struct tm* now_tm = gmtime(&now);
    srand(now_tm->tm_sec);
    long int seed = now_tm->tm_sec*rand()
	+ now_tm->tm_min*rand()
	+ now_tm->tm_hour*rand() 
	+ now_tm->tm_mday*rand()
	+ now_tm->tm_mon*rand() 
	+ now_tm->tm_year*rand();

    return seed;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) error("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

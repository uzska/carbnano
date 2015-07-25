#ifndef WALKER_H_
#define WALKER_H_


int RandomWalk(double *Current, gsl_rng *rng, int *N, 
	       double (*iN)[3], double (*fN)[3], int length, int radius,
	       int fineness, int bins);

int calculateBin(double x, double y, double z, int bins);


#endif

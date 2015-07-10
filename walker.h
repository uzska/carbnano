#ifndef WALKER_H_
#define WALKER_H_

/*
 * Updates the position by performing a random walk
 */
int Random_Walk(double *Current, gsl_rng *rng, double i_N[][3], double f_N[][3],
		int *N, double side_length,int bins, int fineness);

/*
 * Returns a value for phi, which depends on the Current y-value,
 * i.e. the boundaries
 */
double get_phi(double y, gsl_rng *rng, double radius, double side_length);



#endif

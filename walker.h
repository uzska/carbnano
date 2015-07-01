#ifndef WALKER_H_
#define WALKER_H_

/*
 * Returns the random walk of a walker, given the locations
 * of the CNs. Matrix has dimensions TIME x 3, each row is
 * a spatial location: x,y,z
 */
int iterate_Random_Walk(double (*Walk)[3], FILE *Nanotube_File, gsl_rng *rng,
			int TIME, int n_Walks, double side_length, int faces);

/*
 * Helper function for iterateRandomWalk. Returns the position
 * vector of one random walk
 */
int Random_Walk(double Next[3], double *Current, gsl_rng *rng, 
		FILE *f, double side_length);

/*
 * Returns a value for phi, which depends on the Current y-value,
 * i.e. the boundaries
 */
double get_phi(double y, gsl_rng *rng, double radius, double side_length);



#endif

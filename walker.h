#ifndef WALKER_H_
#define WALKER_H_

/*
 * Returns the random walk of a walker, given the locations
 * of the CNs. Matrix has dimensions TIME x 3, each row is
 * a spatial location: x,y,z
 */
//int iterate_Random_Walk(double (*Walk)[3], FILE *Nanotube_File, gsl_rng *rng, int iteration, int no,
//			int Times, int Rates, int process_id, int n_Walks, double side_length, int faces);

/*
 * Helper function for iterateRandomWalk. Returns the position
 * vector of one random walk
 */
int Random_Walk(double *Current, gsl_rng *rng, double i_N[][3], double f_N[][3],
		int *N, double side_length,int bins, int fineness);

/*
 * Returns a value for phi, which depends on the Current y-value,
 * i.e. the boundaries
 */
double get_phi(double y, gsl_rng *rng, double radius, double side_length);



#endif

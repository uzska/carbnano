/*
 * walker.c
 */
#include <gsl_rng.h>
#include <gsl_cdf.h>
#include <math.h>

/* 
 * Get New Position on the Nanotube
 */
int Redistribute(double i[3], double f[3], double *Current,
		 int radius, int fineness, int bins, gsl_rng *rng) {
  double p1 = gsl_rng_uniform(rng);
  double p2 = gsl_rng_uniform(rng);
  double p3 = gsl_rng_uniform(rng);
  double off = (1.0*radius) / (1.0*bins*fineness);

  double x1, y1, z1;
  double x2, y2, z2;
  
  if (i[0] < f[0]) {x1 = i[0]; x2 = f[0];} else {x1 = f[0]; x2 = i[0];} 
  if (i[1] < f[1]) {y1 = i[1]; y2 = f[1];} else {y1 = f[1]; y2 = i[1];} 
  if (i[2] < f[2]) {z1 = i[2]; z2 = f[2];} else {z1 = f[2]; z2 = i[2];} 

  Current[0] = x1 - off + p1*(x2 - x1 + 2.0*off);
  Current[1] = y1 - off + p1*(y2 - y1 + 2.0*off); 
  Current[2] = z1 - off + p1*(z2 - z1 + 2.0*off);
  
  return 0;
}


/*
 * Bins are numbered from 0 to n_bins_side^3-1
 * The bins include the 3 lower sides of the cube, 
 * e.g. [0.1,0.2) x [0.2,0.3) x [0.3,0.4)
 *
 * The bins on the boundaries (x=1, y=1, z=1) do contain 
 * their respective upper sides
 */
int calculateBin(double x, double y, double z, int bins) {
  int xx = x * bins; // -> x / (sidelength/n_bin_side)
  int yy = y * bins;
  int zz = z * bins;

  // If a coordinate is on a boundary, we say it's on the upper side
  // of a bin
  if (xx == bins) xx--;
  if (yy == bins) yy--;
  if (zz == bins) zz--;

  return yy*bins*bins + xx + zz*bins;
}

/*
 * Performs a Random Walk:
 * First, check if the walker is already inside a nanotube. 
 * Either redistribute anywhere on the nanotube, or pop out
 * of the nanotube
 *
 * If the walker wants to get into a nanotube, we check the
 * probability. If insufficient, the walker stays where it is.
 * If sufficient, the walker enters the nanotube
 *
 * returns 0 if in the matrix, 1 if in a nanotube
 */
int RandomWalk(double *Current, gsl_rng *rng, int *N,
	       double (*iN)[3], double (*fN)[3], int length, int radius,
	       int fineness, int bins) {
  int in = 0; 
  double Dm = 0.02; // diffusivity
  double sigma = sqrt(2 * Dm);   // standard deviation sqrt(2*D_m*delta_t)
  double fm_cn = 0.5;
  double fcn_m = sigma * fineness * bins * 
    (2.0/(1.0*length)+ 4.0/(1.0*radius)) * fm_cn;

  double x = Current[0]; 
  double y = Current[1];
  double z = Current[2];

  /* check if we're inside a nanotube */
  int lookup = 0;
  if (N) {lookup = N[calculateBin(x,y,z,fineness*bins)];}
  if (0 != lookup) {
    double p = gsl_rng_uniform(rng);
    if (p < 1.0-fcn_m) {
      //redistribute in nanotube
      Redistribute(iN[lookup-1],fN[lookup-1],Current,radius,
		   fineness,bins,rng);
      return 1;
    }
  }

  // Random Walk according to a Normal Distribution,
  // potential new positions x y z
  double s1 = gsl_rng_uniform_pos(rng);
  x += (sigma / sqrt(3)) * gsl_cdf_ugaussian_Qinv(s1);

  double s2 = gsl_rng_uniform_pos(rng);
  z += sigma / sqrt(3) * gsl_cdf_ugaussian_Qinv(s2);
  
  double s3;
  while (1) {
    s3 = gsl_rng_uniform_pos(rng);
    if ((y + sigma / sqrt(3) * gsl_cdf_ugaussian_Qinv(s3) >= 0) &&
	(y + sigma / sqrt(3) * gsl_cdf_ugaussian_Qinv(s3) <= 1)) {
      break;
    } 
  }
  y += sigma / sqrt(3) * gsl_cdf_ugaussian_Qinv(s3);


  // periodic boundary conditions
  while (x > 1.0) { x = x - 1.0;}
  while (z > 1.0) { z = z - 1.0;}
  while (x < 0) { x = x + 1.0;}
  while (z < 0) { z = z + 1.0;}

  /* Check if our new position is located in a nanotube  */
  if (N) {lookup = N[calculateBin(x,y,z,fineness*bins)];}
  if (0 != lookup) {
    double p = gsl_rng_uniform(rng);
    if (p < 1.0-fm_cn) {
      return 0; // don't enter nanotube, walker stays put
    }
    else {in = 1;}
  }

  // updates position, enters nanotube
  Current[0] = x;
  Current[1] = y;
  Current[2] = z;
 
  return in;  
}

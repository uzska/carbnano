/*
 * walker.c
 */
#include <gsl_rng.h>
#include <gsl_sf_result.h>
#include <math.h>

/*
 * Inverse erf
 */ 
double erfinv(double x){
  double z;
  double a  = 0.147;
  double the_sign_of_x;
  if(0==x) {
    the_sign_of_x = 0;
  } else if(x>0){
    the_sign_of_x = 1;
  } else {
    the_sign_of_x = -1;
  }

  // cases for x = 1,-1
  if (x == 1) {
    //infinity
    return INFINITY;
  }
  if (x == -1) {
    //neg infinity
    return -INFINITY;
  }

  if(0 != x) {
    double ln_1minus_x_sqrd = log(1-x*x);
    double ln_1minusxx_by_a = ln_1minus_x_sqrd / a;
    double ln_1minusxx_by_2 = ln_1minus_x_sqrd / 2;
    double ln_etc_by2_plus2 = ln_1minusxx_by_2 + (2/(M_PI * a));
    double first_sqrt = sqrt((ln_etc_by2_plus2*ln_etc_by2_plus2)-ln_1minusxx_by_a);
    double second_sqrt = sqrt(first_sqrt - ln_etc_by2_plus2);
    z = second_sqrt * the_sign_of_x;
  } else { // x is zero                                                                                                                              
    z = 0;
  }
  return z;
}

/* 
 * Get New Position on the Nanotube
 */
int RandomPosBounds(double i[3], double f[3], double *Current,
		    gsl_rng *rng) {
  double p1 = gsl_rng_get(rng) / (1.0*gsl_rng_max(rng));
  double p2 = gsl_rng_get(rng) / (1.0*gsl_rng_max(rng)); 
  double p3 = gsl_rng_get(rng) / (1.0*gsl_rng_max(rng));

  Current[0] = i[0] + p1*(f[0] - i[0]);
  Current[1] = i[1] + p2*(f[1] - i[1]);
  Current[2] = i[2] + p3*(f[2] - i[2]);
  
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
int calculate_Bin(double x, double y, double z, int n_bin_side, double side_length) {
  int xx = x * n_bin_side / side_length; // -> x / (sidelength/n_bin_side)
  int yy = y * n_bin_side / side_length;
  int zz = z * n_bin_side / side_length;

  // If a coordinate is on a boundary, we say it's on the upper side
  // of a bin
  if (xx == n_bin_side) xx--;
  if (yy == n_bin_side) yy--;
  if (zz == n_bin_side) zz--;

  return yy*n_bin_side*n_bin_side + xx + zz*n_bin_side;
}

/*
 * Returns the azemuthal angle phi, and does not allow a walker to exit
 * from the heated sides.
 */
/*
double get_phi(double y, gsl_rng *rng, double radius, double side_length) {
  // The walker is getting close to the boundary y=0, so we bounce it back
  if (y < radius) {
    return gsl_rng_get(rng) / (1.0*gsl_rng_max(rng)) * M_PI;
  }
  // The walker is getting close to the boundary y=1, so we bounce it back
  else if (y > (side_length - radius)) {
    return (gsl_rng_get(rng) / (1.0*gsl_rng_max(rng)) * M_PI) + M_PI;
  }
  else {
    return gsl_rng_get(rng) / (1.0*gsl_rng_max(rng)) * 2 * M_PI;
  }
}
*/

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
int Random_Walk(double *Current, gsl_rng *rng, double i_N[][3], 
		double f_N[][3], int *N, double side_length, 
		int radius, int length,
		int bins, int fineness) {

  int in = 0; 

  double Dm = 0.02; // diffusivity
  double sigma = sqrt(2 * Dm);   // standard deviation sqrt(2*D_m*delta_t)

  // prob. that a walker crosses from the matrix to the nanotube
  // fm_cn = 4 / pC C_m R_bd
  double fm_cn = 0.5;
  // prob. that a walker crosses from the nanotube to the matrix
  // fcn_m = sigma A_c / V_c * fm_cn
  double fcn_m = sigma * fineness * bins * 
    (2.0 / (1.0*length) + 4.0 / (1.0*radius)) * fm_cn;

  // sa / vol = 2*pi radius^2  + 2 pi radius length / pi radius^2 length
  // 2 /length + 2/radius
  // 2 * radius^2 + 4 * length*radius / radius^2 * length
  // 2/length + 4/radius

  double x = Current[0]; double y = Current[1]; double z = Current[2];

  // if in nanotube, check which nanotube we're in
  int lookup = N[calculate_Bin(x,y,z,fineness*bins,side_length)];
  if (lookup != 0) {
    double prob = gsl_rng_get(rng) / (1.0*gsl_rng_max(rng));
    if (prob < fcn_m) {
      // redistribute and return
      RandomPosBounds(i_N[lookup-1],f_N[lookup-1],Current,rng);
      return 1;
    }
  }

  // Random Walk according to a Normal Distribution,
  // potential new positions x y z
  double s1 = gsl_rng_get(rng) / (1.0 *gsl_rng_max(rng));
  x += sigma / sqrt(3) * erfinv(2*s1-1);
  
  double s2 = gsl_rng_get(rng) / (1.0 *gsl_rng_max(rng));
  z += sigma / sqrt(3) * erfinv(2*s2-1);
  
  double s3;
  while (1) {
    s3 = gsl_rng_get(rng) / (1.0 *gsl_rng_max(rng));
    
    if (y + sigma / sqrt(3) * erfinv(2*s3-1) >= 0 &&
	y + sigma / sqrt(3) * erfinv(2*s3-1) <= 1) {
      break;
    } 
  }    
  y += sigma / sqrt(3) * erfinv(2*s3-1);

  // periodic boundary conditions
  if (x > side_length) {x -= side_length;}
  if (x < 0) {x += side_length;}
  if (z > side_length) {z -= side_length;}
  if (z < 0) {z += side_length;}

  // check if walker enters nanotube
  lookup = N[calculate_Bin(x,y,z,fineness*bins,side_length)];
  if (lookup != 0) {
    double prob = gsl_rng_get(rng) / (1.0*gsl_rng_max(rng));
    // exits and walker does not get new position if 
    // probability is lacking
    if (prob < fm_cn) {
      return 0;
    }
    else { in = 1;}
  }

  Current[0] = x;
  Current[1] = y;
  Current[2] = z;

  /*
  // update x, y, z
  Current[0] = radius*sin(theta)*cos(phi) + Current[0];
  Current[1] = radius*sin(theta)*sin(phi) + Current[1];
  Current[2] = radius*cos(theta) + Current[2];
  
  // periodic boundary conditions
  if (Current[0] > side_length) {
    Current[0] -= side_length;
  }
  if (Current[2] > side_length) {
    Current[2] -= side_length;
  }
  if (Current[0] < 0) {
    Current[0] += side_length;
  }
  if (Current[2] < 0) {
    Current[2] += side_length;
  }
  */
  return in;  
}

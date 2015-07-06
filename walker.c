/*
 * walker.c
 */
#include <gsl_rng.h>
#include <math.h>

/*
 * Bins are numbered from 0 to n_bins_side^3-1
 * The bins include the 3 lower sides of the cube, 
 * e.g. [0.1,0.2) x [0.2,0.3) x [0.3,0.4)
 *
 * The bins on the boundaries (x=1, y=1, z=1) do contain 
 * their respective upper sides
 */
int calculate_Bin(double x, double y, double z, int n_bin_side, double side_length) {
  int xx = x * n_bin_side / side_length; // x / (sidelength/n_bin_side)
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

/*
 * Next is the new position of the walker, given it's 
 * old position Current and the location of the nanotubes
 */
int Random_Walk(double *Current, gsl_rng *rng, 
		char *N, double side_length, int bins, int fineness) {


  double theta;
  double phi;
  double radius = 1.0/199.0;
  double nanoBarrier = 0.5;

  // TODO: somehow have to know what nanotube we are in

  // if in nanotube
  if (0) {
    // either redistribute and return or get new position
    double prob = gsl_rng_get(rng) / (1.0*gsl_rng_max(rng));
    if (prob > nanoBarrier) {
      
      return 1;
    }
  }


  // calculate angles  
  theta = (gsl_rng_get(rng) / (1.0*gsl_rng_max(rng))) * M_PI;
  phi = get_phi(Current[1], rng, radius, side_length);
  
  // potential new positions
  double x = radius*sin(theta)*cos(phi) + Current[0];
  double y = radius*sin(theta)*sin(phi) + Current[1];
  double z = radius*cos(theta) + Current[2];
  
  if (N[calculate_Bin(x,y,z,fineness*bins,side_length)] == 'x') {
    // probability that walker goes in nanotube
    double prob = gsl_rng_get(rng) / (1.0*gsl_rng_max(rng));
    // exits and walker does not get new position if 
    // probability is lacking
    if (prob < nanoBarrier) {
      return 1;
    }
  }
  
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

  return 0;  
}

/*
int iterate_Random_Walk(double (*Walk)[3], FILE *Nanotube_File, gsl_rng *rng, int iteration, int n_proc,
			int Times, int Rates, int id, int walks, double side_length, int faces) {

  int i, j, k;
  double Next[3];
  if (Rates == 0) {
    for (i = 0; i < walks; i++) {
      int x = Times-1;
      for (j = 0; j < 2*faces; j++) {
	for (k = 0; k < x; k++) {
	  Random_Walk(Next, *(Walk + k + j*(x+1) + i*2*faces*(x+1)), rng, Nanotube_File, side_length);
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][0] = Next[0];
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][1] = Next[1];
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][2] = Next[2];
	}
      }
    }
  }
  else if (Rates > 0) {
    for (i = 0; i < walks; i++) {
      int x = Times - (id+iteration*n_proc)*walks/Rates - i/Rates - 1;
      for (j = 0; j < 2*faces; j++) {
	for (k = 0; k < x; k++) {
	  Random_Walk(Next, *(Walk + k + j*(x+1) + i*2*faces*(x+1)), rng, Nanotube_File, side_length);
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][0] = Next[0];
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][1] = Next[1];
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][2] = Next[2];
	}
      }
    }
  }
  else {
    for (i = 0; i < walks; i++) {
      int x = Times + Rates*(iteration*n_proc*walks + id*walks) + i*Rates - 1;
      for (j = 0; j < 2*faces; j++) {
	for (k = 0; k < x; k++) {
	  Random_Walk(Next, *(Walk + k + j*(x+1) + i*2*faces*(x+1)), rng, Nanotube_File, side_length);	  
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][0] = Next[0];
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][1] = Next[1];
	  Walk[1 + k + j*(x+1) + i*2*faces*(x+1)][2] = Next[2];
	}
      }
    }
  }

  return 0;
}
*/

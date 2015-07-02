/*
 * walker.c
 */
#include <gsl_rng.h>
#include <gsl_matrix.h>
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
int Random_Walk(double Next[3], double *Current, gsl_rng *rng, 
		FILE *f, double side_length) {


  double theta;
  double phi;
  double radius = 1.0/70.0;

  // look at nanotube arrays and act accordingly
  if (0) {
  }
  // random walk, no restraints on it's new position
  else {
    theta = (gsl_rng_get(rng) / (1.0*gsl_rng_max(rng))) * M_PI;
    phi = get_phi(Current[1], rng, radius, side_length);
  }

  // update x, y, z
  Next[0] = radius*sin(theta)*cos(phi) + Current[0];
  Next[1] = radius*sin(theta)*sin(phi) + Current[1];
  Next[2] = radius*cos(theta) + Current[2];

  // periodic boundary conditions
  if (Next[0] > side_length) {
    Next[0] -= side_length;
  }
  if (Next[2] > side_length) {
    Next[2] -= side_length;
  }
  if (Next[0] < 0) {
    Next[0] += side_length;
  }
  if (Next[2] < 0) {
    Next[2] += side_length;
  }

  return 0;
}

int iterate_Random_Walk(double (*Walk)[3], FILE *Nanotube_File, gsl_rng *rng,
			int TIME, int Times, int Rates, int id, int n_Walks, double side_length, int faces) {

  int i, j, k;
  double Next[3];
  if (Rates = 0) {
    for (i = 0; i < TIME-1; i++) {
      for (j = 0; j < n_Walks; j++) {
	for (k = 0; k < 2*faces; k++) {
	  Random_Walk(Next, *(Walk + i + j*TIME*2*faces + k*TIME), rng, Nanotube_File, side_length);
	  Walk[i + 1 + j*TIME*2*faces + k*TIME][0]= Next[0];
	  Walk[i + 1 + j*TIME*2*faces + k*TIME][1]= Next[1];
	  Walk[i + 1 + j*TIME*2*faces + k*TIME][2]= Next[2];
	  
	  //Walk[i+1+j*TIME][1]= Next[1];
	  //Walk[i+1+j*TIME][2]= Next[2];
	}
      }
    }
  }
  else if (Rates > 0) {
    for (i = 0; i < n_Walks; i++) {
      for (j = 0; j < 2*faces; j++) {
	for (k = 0; k < Times - id*n_Walks/Rates - i/Rates - 1; k++) {

	  int x = Times - id*n_Walks*Rates - i/Rates;
	  Random_Walk(Next, *(Walk + k + j*(Times - id*n_Walks/Rates - i/Rates) + 
			      i*2*faces*(Times - id*n_Walks/Rates - i/Rates)), 
		      rng, Nanotube_File, side_length);

	  Walk[1 + k + j*x + i*2*faces*x][0] = Next[0];
	  Walk[1 + k + j*x + i*2*faces*x][1] = Next[1];
	  Walk[1 + k + j*x + i*2*faces*x][2] = Next[2];
	}
      }
    }
  }
  else {
    for (i = 0; i < n_Walks; i++) {
      for (j = 0; j < 2*faces; j++) {
	for (k = 0; k < Times + id*n_Walks*Rates + i*Rates - 1; k++) {
	  int x = Times + id*n_Walks*Rates + i*Rates;
	  Random_Walk(Next, *(Walk + k + j*x + i*2*faces*x), rng, Nanotube_File, side_length);
	  
	  Walk[1 + k + j*x + i*2*faces*x][0] = Next[0];
	  Walk[1 + k + j*x + i*2*faces*x][1] = Next[1];
	  Walk[1 + k + j*x + i*2*faces*x][2] = Next[2];
	}
      }
    }
  }
		 
  return 0;
}


#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <mpi.h>
#include <gsl_rng.h>

#include "nanotube.h"
#include "walker.h"

/*
 * Least Squares Method to Calculate Slope
 */
int computeSlope(int *T, int bins, double *slope, double *y_int) {
  int i;
  double X_avg = 0;
  double Y_avg = 0;
  for (i = 0; i < bins; i++) {
    X_avg += i; 
    Y_avg += T[i];
  }

  X_avg /= (1.0*bins);
  Y_avg /= (1.0*bins);

  double denom = 0;
  for (i = 0; i < bins; i++) {
    *slope += (i-X_avg) * (T[i] - Y_avg);
    denom += (i-X_avg) * (i-X_avg);
  }

  *slope /= denom;
  *y_int = Y_avg - (*slope)*X_avg;
  
  return 0;
}


/*
 * Returns the max of three numbers
 */
int getMax(int x, int y, int z) {
  if (x >= y && x >= z) {return x;}
  else if (y >= x && y >= z) {return y;}
  else {return z;}
}

/*
 * Starting positions of the walkers.
 */
int initWalkArray(double (*W)[3], double F[][2], int faces) {
  int i;
  for (i = 0; i < 2*faces; i++) {
    W[i][0] = F[i/2][0];
    W[i][1] = i%2;
    W[i][0] = F[i/2][1];
  }
  return 0;
}

/*
 * Offset algorithm to determine how many time steps of the walker's 
 * position should be calculated
 */
int initTime(int Time, int Rate, int id, int walks, int currentWalk) {
  if (Rate > 0) {
    return (Time - id*walks/Rate - currentWalk/Rate);
  }
  else if (Rate < 0) {
    return (Time + id*walks*Rate + currentWalk*Rate);
  }
  else {
    return Time;
  }
}

/* 
 *  Main Loop 
 */
int main(int argc, char *argv[]) {

  /* World Variables */
  double side_length = 1;
  int n_bin_side = 10; //per side 
  int faces = n_bin_side*n_bin_side;
  int DIM = 3;
  
  double interval = side_length / ((double)n_bin_side);
  int fineness = 3; //the much finer the nanotube grid is than the temp bins

  /* Times, Rate, and walks  */
  int Times[] = {20,40};
  int Rate = 1;
  int len_Times = sizeof(Times)/sizeof(Times[0]);

  int walks = 20; // number of different walks each process simulates
  
  /* iterator variables  */
  int i,j,k,p,s;
  
  int m;
  
  /*
   * A list of the initial x and z positions of the walkers. 
   * The list first proceeds across from x=0 to x=1 and then upwards from
   * z=0 to z=1.
   */
  double FacesXZ[faces][2];
  for (i = 0; i < n_bin_side; i++) {
    for (j = 0; j < n_bin_side; j++) {
      FacesXZ[i*n_bin_side + j][0] = 1.0/(2.0*n_bin_side) + (j*1.0)/(1.0*n_bin_side);
      FacesXZ[i*n_bin_side + j][1] = 1.0/(2.0*n_bin_side) + (i*1.0)/(1.0*n_bin_side);
    }
  }
  
  /* Spawn # of processes */
  int process_id;
  int n_processes;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&n_processes);  
  MPI_Comm_rank(MPI_COMM_WORLD,&process_id);


  /* Error Checking, no of processes needed */
  if (Rate > 0) {
    if (n_processes*walks < Times[len_Times-1]*Rate) {
      fprintf(stderr, "Rate > 1: Need n_processes * walks to be >= %d\n", 
	      Times[len_Times-1]*Rate);
      exit(EXIT_FAILURE);
    }
  } 
  else if (Rate < 0) {
    if (n_processes*walks < Times[len_Times-1]/abs(Rate)) {
      fprintf(stderr, "Rate < 1: Need n_processes * walks to be >= %d\n", 
	      Times[len_Times-1]/abs(Rate));
      exit(EXIT_FAILURE);
    }
  }

  /* Total array */
  int (*Total_cs)[len_Times] = NULL;
  if (process_id == 0) {
    Total_cs = malloc(sizeof(*Total_cs) * n_bin_side);
  }
  
  /* Carbon Nanotube Lookup Table */
  int *Nanotubes = calloc(fineness * fineness * 
			  fineness * n_bin_side * 
			  n_bin_side * n_bin_side,
			   sizeof(int));
  int n_tubes = 1;
  double i_Nanotubes[n_tubes][DIM];
  double f_Nanotubes[n_tubes][DIM];

  i_Nanotubes[0][0] = 0;
  i_Nanotubes[0][1] = 0;
  i_Nanotubes[0][2] = 0;

  f_Nanotubes[0][0] = 0;
  f_Nanotubes[0][0] = 1;
  f_Nanotubes[0][0] = 1;
  
  /* fill in bins between initial and final bins */
  for (i = 0; i < n_tubes; i++) {
    // determine the x,y, and z bin numbers
    int x = i_Nanotubes[i][0] * n_bin_side * fineness / side_length;
    int y = i_Nanotubes[i][1] * n_bin_side * fineness / side_length;
    int z = i_Nanotubes[i][2] * n_bin_side * fineness / side_length;

    int xx = f_Nanotubes[i][0] * n_bin_side * fineness / side_length;
    int yy = f_Nanotubes[i][1] * n_bin_side * fineness / side_length;
    int zz = f_Nanotubes[i][2] * n_bin_side * fineness / side_length;
        
    int max = getMax(xx-x,yy-y,zz-z);

    double dx = (xx - x)/max;
    double dy = (yy - y)/max;
    double dz = (zz - z)/max;

    while (x <= f_Nanotubes[i][0]) {
      int b = x + y*(n_bin_side*fineness)*(n_bin_side*fineness) + z*(n_bin_side);
      Nanotubes[b] = i+1;
      x += dx; y += dy; z += dz;
    }
  }

  /* RNG */
  const gsl_rng_type *T = gsl_rng_taus;
  gsl_rng *rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL)+process_id);

  /* Random Walk, holds only the current position  */  
  double (*Walk)[DIM] = malloc(sizeof(*Walk) * 2 * faces);
  if (Walk) {
    initWalkArray(Walk,FacesXZ,faces);
  }
  else {
    fprintf(stderr,"Process %d: Walk array could not be allocated\n", process_id);
    exit(EXIT_FAILURE);
  }
  
  /*
   * Random Walk and walker data collection
   *
   * Each process can generate up to "walks" random walkers. 
   * To determine if a walker needs to be generated, we calculate T,
   * which is how long the walker has alive for, given that we want
   * to measure the system at Times[0], Times[1], ... Times[len_Times-1].
   *
   * If T is nonnegative, we do record the walker's position after T 
   * random walks.
   *
   * We then record the walker's position at the next Time, which means
   * we have to do Times[m+1]-Times[m] more random walks.
   */
  int bin;
  int (*SubTotal_cs)[len_Times] = calloc(n_bin_side,sizeof(*SubTotal_cs));

  for (i = 0; i < walks; i++) {
    // reinitialize Walk array every random walk that we perform
    if (i) {initWalkArray(Walk, FacesXZ, faces);}
    // How long we need to run the simulation for
    int T = initTime(Times[0],Rate,process_id,walks,i);

    // go through each time in the Time array
    for (m = 0; m < len_Times && T >= 0; m++) {
      // random walk for each walker generated on the 2*n_bin_side^2 bins
      for (j = 0; j < 2*faces && T >= 0; j++) {
	for (k = 0; k < T; k++) {
	  Random_Walk(*(Walk+j),rng,i_Nanotubes,f_Nanotubes,
		      Nanotubes,side_length,n_bin_side,fineness);
	}	
	bin = calculate_Bin(Walk[j][0], Walk[j][1], Walk[j][2], 
			    n_bin_side, side_length);

	if (j%2==0) {SubTotal_cs[bin/faces][m]++;} 
	else {SubTotal_cs[bin/faces][m]--;}
      }
      // get new T to determine how much longer we need to run the 
      // simulation until we record the position
      if (m!=len_Times-1) {T = Times[m+1] - Times[m];}
    }
  }

  /* Synchronize processes with a barrier */
  MPI_Barrier(MPI_COMM_WORLD);  

  // free memory
  free(Walk);
    
  MPI_Reduce(SubTotal_cs, Total_cs, n_bin_side*len_Times, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
  // free memory
  free(SubTotal_cs);

  /* print out cross section totals, slopes, and thermal conductivities */
  if (process_id == 0) {
    FILE *f;
    f = fopen("results.txt","w");

    for (k = 0; k < len_Times; k++) {
      FILE *g;                                                                                                      
      char name_g[FILENAME_MAX];
      snprintf(name_g, sizeof(name_g), "output/Totalcs_%d.csv",k);
      g = fopen(name_g,"w");
      
      for (i = 0; i < n_bin_side; i++) {
	if (i != n_bin_side-1) {fprintf(g,"%d,",Total_cs[i][k]);}
	else {fprintf(g,"%d",Total_cs[i][k]);}
      }
      fclose(g);


      // slope, thermal conductivities
      double slope = 0; double y_int;
      computeSlope(&Total_cs[i][0], n_bin_side, &slope, &y_int);
      double flux; 
      if (Rate > 0) {
	flux = Rate*faces;
      }
      else if (Rate < 0) {
	flux = faces/abs(Rate);
      }
      else {
	flux = 0;
      }
      
      double conduct = -1.0*flux / slope;
      fprintf(f,"slope: %g y-int: %g conductivity: %g\n ",slope,y_int,conduct);
    }
    fclose(f);

  }

  // free memory
  free(Total_cs);

  MPI_Finalize();
  return 0;
}

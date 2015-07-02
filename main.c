#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <mpi.h>
#include <gsl_matrix.h>
#include <gsl_vector.h>
#include <gsl_rng.h>

#include "nanotube.h"
#include "walker.h"

int main(int argc, char *argv[]) {

  /*
   * World Variables
   */
  double side_length = 1;
  int n_bin_side = 10; //per side 
  int faces = n_bin_side*n_bin_side;
  int TIME = 10;
  int DIM = 3;

  double interval = side_length / ((double)n_bin_side);
  int fineness = 5; //the much finer the nanotube grid is than the temp bins

  int Times[] = {20}; // Times to measure at
  int Rates[] = {0};
  int len_Times = sizeof(Times)/sizeof(Times[0]);

  int n_Walks = 3; // number of different walks each process simulates
  int iterations = 5; // no of times a process repeats
 
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

  /* 
   * Spawn # of processes
   */
  int process_id;
  int n_processes;
  MPI_Group group_world;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&n_processes);  
  MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
  MPI_Comm_group(MPI_COMM_WORLD,&group_world);

  /*
   * Error Checking, no of processes needed
   */
  if (Rates[len_Times-1] > 0) {
    if (n_processes*n_Walks*iterations != Times[len_Times-1]*Rates[len_Times-1]) {
      fprintf(stderr, "Rate > 1: Need n_processes * n_Walks * iterations to be at exactly %d\n", 
	      Times[len_Times-1]*Rates[len_Times-1]);
      exit(EXIT_FAILURE);
    }
  } 
  else if (Rates[len_Times-1] < 0) {
    if (n_processes*n_Walks*iterations < TIME/abs(Rates[len_Times-1])) {
      fprintf(stderr, "Rate < 1: Need n_processes * n_Walks * iterations to be at least %d\n", TIME/abs(Rates[len_Times-1]));
      exit(EXIT_FAILURE);
    }
  }

  /*
   * Final Total array 
   */
  int (*FinalTotal_cs)[len_Times] = NULL;
  if (process_id == 0) {
    FinalTotal_cs = malloc(sizeof(*FinalTotal_cs) * n_bin_side);
  }
  
  /*
   * Carbon Nanotube Lookup Table
   */
  FILE *Nanotube_File = NULL;

  //char *Nanotubes = malloc(sizeof(int) * fineness * fineness * fineness * n_bins_side * n_bins_side * n_bins_side);
  //int n_tubes;
  //double i_Nanotubes[n_tubes][n_tubes][n_tubes];  
  //double f_Nanotubes[n_tubes][n_tubes][n_tubes];
  //char *Nanotubes = malloc(sizeof(char) * (n_bin_side*fineness) * (n_bin_side*fineness) * (n_bin_side*fineness));

  /*
   * RNG
   */
  const gsl_rng_type *T = gsl_rng_taus;
  gsl_rng *rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL)+process_id);


  /*
   * Multiple cycles to limit the number of nodes used and to recycle memory
   */
  for (p = 0; p < iterations; p++) {
    
    /*
     * allocate memory for Walk array, which records Times + 1 steps (including the start)
     * Walk is a 2-D array, with TIME rows and DIM cols
     * Each row represents the position in (x,y,z) in a different time
     */
    double (*Walk)[DIM];
    if (Rates[len_Times-1] > 0) {
      s=0;
      for (i = 0; i < n_Walks; i++) {
	s += 2*faces * (Times[len_Times-1] + 1
			- n_Walks/Rates[len_Times-1] * (p*n_processes + process_id) 
			- i/Rates[len_Times-1]);
      }
      Walk = malloc(sizeof(*Walk) * s);
    }
    else if (Rates[len_Times-1] < 0) {
      s = 0;
      for (i = 0; i < n_Walks; i++) {
	s = s + (Times[len_Times-1] + process_id*n_Walks*Rates[len_Times-1] + i*Rates[len_Times-1]) * 2 * faces;
      }
      Walk = malloc(sizeof(*Walk) * s);
      printf("%d\n",s);
    }
    else {
      s = (1 + Times[len_Times-1]) * n_Walks * 2 * faces;
      Walk = malloc(sizeof(*Walk) * s);
    }
    
    /* 
     * initialize the starting walkers
     * Give processes their start location
     */
    if (Walk) {
      s = 0;
      if (Rates[len_Times-1] > 0) {
	  for (i = 0; i < n_Walks; i++) {
	    for (j = 0; j < faces*2; j++) {	      
	      Walk[s][0] = FacesXZ[j/2][0];
	      Walk[s][1] = j%2;
	      Walk[s][2] = FacesXZ[j/2][1];
	      s += 1+Times[len_Times-1] -  (p*n_processes + process_id)*n_Walks/Rates[len_Times-1] - i/Rates[len_Times-1];
	    }
	  }
      }
      else if (Rates[len_Times-1] < 0) {
	for (i = 0; i < n_Walks; i++) {
	  for (j = 0; j < faces*2; j++) {
	    Walk[s][0] = FacesXZ[j/2][0];
	    Walk[s][1] = j%2;
	    Walk[s][2] = FacesXZ[j/2][1];
	    s += 1 + Times[len_Times-1] + process_id*n_Walks*Rates[len_Times-1] + i*Rates[len_Times-1];
	  }
	}
      }
      else {
	for (i = 0; i < n_Walks; i++) {
	  for (j = 0; j < faces*2; j++) {
	    Walk[s][0] = FacesXZ[j/2][0];
	    Walk[s][1] = j%2;
	    Walk[s][2] = FacesXZ[j/2][1];
	    s += 1+Times[len_Times-1];
	  }
	}
      }
    }
    else {
      fprintf(stderr,"Process %d: Walk array could not be allocated\n", process_id);
      exit(EXIT_FAILURE);
    }
    
    
    /*
     * Generate a Random Walk for each process
     */
    iterate_Random_Walk(Walk, NULL, rng, p, n_processes, Times[len_Times-1], 
			Rates[len_Times-1], process_id, n_Walks, side_length, faces);

    /*
     * Synchronize processes with a barrier
     */
    MPI_Barrier(MPI_COMM_WORLD);  
    
    int binOfWalker;
    int (*SubTotal_cs)[len_Times] = calloc(n_bin_side,sizeof(*SubTotal_cs));
    //int (*SubTotal)[len_Times] = calloc(n_bin_side*n_bin_side*n_bin_side,sizeof(*SubTotal));

    /*
     * Go through each Time to record walker data
     */
    for (k = 0; k < len_Times; k++) {
      double x,y,z;
      
      /* Rates > 0 */
      if (Rates[k] > 0) {
	int cur = Times[k] - (p*n_processes + process_id)*n_Walks/Rates[k];
	for (i = 0; i < n_Walks && 1; i++) {
	  for (j = 0; j < 2*faces; j++) {
	    x = Walk[cur][0];
	    y = Walk[cur][1];
	    z = Walk[cur][2];
	    binOfWalker = calculate_Bin(x, y, z, n_bin_side, side_length);
	    if (j%2 == 0) {SubTotal_cs[binOfWalker/faces][k]++;} else {SubTotal_cs[binOfWalker/faces][k]--;}
	    cur += 1+Times[k] - (p*n_processes + process_id)*n_Walks/Rates[k] - i/Rates[k];
	  }
	}
      }
      
      /* Rates < 0 */
      else if (Rates[k] < 0) {
	for (i = 0; i < n_Walks && (-p*n_processes*Rates[k]+ -process_id*Rates[k] + -i*Rates[k]) <= Times[k]; i++) {
	  for (j = 0; j < 2*faces; j++) {
	    binOfWalker = calculate_Bin(Walk[i*2*faces*TIME + j*TIME + Times[k] + p*n_processes*Rates[k] + process_id*Rates[k] + i*Rates[k]][0],
					Walk[i*2*faces*TIME + j*TIME + Times[k] + p*n_processes*Rates[k] + process_id*Rates[k] + i*Rates[k]][1],
					Walk[i*2*faces*TIME + j*TIME + Times[k] + p*n_processes*Rates[k] + process_id*Rates[k] + i*Rates[k]][2],
					n_bin_side, side_length);
	    if (j%2 == 0) {SubTotal_cs[binOfWalker/faces][k]++;} else {SubTotal_cs[binOfWalker/faces][k]--;}
	    //if (j%2 == 0) {SubTotal[binOfWalker][k]++;} else {SubTotal[binOfWalker][k]--;}
	  }
	}
      }
      
      /* Rates == 0 */
      else if (Rates[k] == 0) {
	int cur = Times[k];
	for (i = 0; i < n_Walks; i++) {
	  for (j = 0; j < 2*faces; j++) {
	    x = Walk[cur][0];
	    y = Walk[cur][1];
	    z = Walk[cur][2];

	    binOfWalker = calculate_Bin(x,y,z, n_bin_side, side_length);
	    if (j%2 == 0) {SubTotal_cs[binOfWalker/faces][k]++;} else {SubTotal_cs[binOfWalker/faces][k]--;}
	    cur += 1+Times[len_Times-1];
	  }
	}
      }
      
    }
    
    // free memory
    free(Walk);

    int (*Total_cs)[len_Times] = NULL;
    Total_cs = malloc(sizeof(*Total_cs) * n_bin_side);
    
    MPI_Reduce(SubTotal_cs, Total_cs, n_bin_side*len_Times, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (process_id == 0) {
      for (i = 0; i < n_bin_side; i++) {
	for (j = 0; j < len_Times; j++) {
	  FinalTotal_cs[i][j] = FinalTotal_cs[i][j] + Total_cs[i][j];
	}
      }
    }


    // free memory
    //free(SubTotal);
    free(SubTotal_cs);
    free(Total_cs);
  } // end for
  
  if (process_id == 0) {
    // print out cross section totals to a file
    for (k = 0; k < len_Times; k++) {
      FILE *g;                                                                                                      
      char name_g[FILENAME_MAX];
      snprintf(name_g, sizeof(name_g), "output/Totalcs_%d.csv",k);
      g = fopen(name_g,"w");
      
      for (i = 0; i < n_bin_side; i++) {
	if (i != n_bin_side-1) {fprintf(g,"%d,",FinalTotal_cs[i][k]);}
	else {fprintf(g,"%d",FinalTotal_cs[i][k]);}
      }
      fclose(g);
    }
    
  }

  // free memory
  free(FinalTotal_cs);
  
  MPI_Finalize();
  return 0;
}

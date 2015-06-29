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
  int n_bin_side = 45; //per side 
  int TIME = 1000;
  int DIM = 3;
  int faces = n_bin_side*n_bin_side;
  double interval = side_length / ((double)n_bin_side);

  int Times[] = {9,99,199,299,399,499,599,699,799,899,999};
  int Rates[] = {2,2,2,2,2,2,2,2,2,2,2};
  int len_Times = sizeof(Times)/sizeof(Times[0]);

  int n_Walks = 24;

  int i;
  int j;
  int k;
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
    if (n_processes*n_Walks < TIME*Rates[len_Times-1]) {
    //if (process_id == 0) {
	fprintf(stderr, "Rate > 1: Need n_processes * n_Walks to be at least %d\n", TIME*Rates[len_Times-1]);
      //}
      exit(EXIT_FAILURE);
    }
  } 
  else if (Rates[len_Times-1] < 0) {
    if (n_processes*n_Walks < TIME/abs(Rates[len_Times-1])) {
    //if (process_id == 0) {
	fprintf(stderr, "Rate < 1: Need n_processes * n_Walks to be at least %d\n", TIME/abs(Rates[len_Times-1]));
      //}
      exit(EXIT_FAILURE);
    }
  } 
  
  /*
   * Carbon Nanotube Lookup Table
   */
  FILE *Nanotube_File = NULL;
  //int *Nanotubes = malloc(sizeof(int) * );


  /*
   * RNG
   */
  const gsl_rng_type *T = gsl_rng_taus;
  gsl_rng *rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL)+process_id);

  /*
   * Give processes their start location
   *
   * Walk is a 2-D array, with TIME rows and DIM cols
   * Each row represents the position in (x,y,z) in a different time
   */  
  double (*Walk)[DIM] = malloc(sizeof(*Walk) * TIME * n_Walks * 2 * faces);
  if (Walk) {
    for (i = 0; i < n_Walks; i++) {
      // process_id/2 is an int
      // This trick allows us to have a hot and cold walker
      // injected from opposite sides, i.e. (x,0,z) and (x,1,z).
      for (j = 0; j < faces*2; j++) {
	Walk[i*TIME*2*faces+j*TIME][0] = FacesXZ[j/2][0];
	Walk[i*TIME*2*faces+j*TIME][1] = j % 2;
	Walk[i*TIME*2*faces+j*TIME][2] = FacesXZ[j/2][1];
      }
    }
  }
  else {
    //if (process_id == 0) {
      fprintf(stderr,"Walk array could not be allocated\n");
      //}
    exit(EXIT_FAILURE);
  }

  /*
   * Generate a Random Walk for each process
   */
  iterate_Random_Walk(Walk, Nanotube_File, rng, TIME, n_Walks, side_length, faces);
  
  /*
   * Synchronize processes with a barrier
   */
  MPI_Barrier(MPI_COMM_WORLD);  

  int binOfWalker;
  int (*SubTotal_cs)[len_Times] = calloc(n_bin_side,sizeof(*SubTotal_cs));
  //int (*SubTotal)[len_Times] = calloc(n_bin_side*n_bin_side*n_bin_side,sizeof(*SubTotal));

  for (k = 0; k < len_Times; k++) {
    if (Rates[k] > 0) {
      for (i = 0; i < n_Walks && (process_id*n_Walks/Rates[k]+i/Rates[k]) <= Times[k]; i++) {
	for (j = 0; j < 2*faces; j++) {
	  binOfWalker = calculate_Bin(Walk[i*2*faces*TIME + j*TIME + Times[k] - process_id*n_Walks/Rates[k] -i/Rates[k]][0],
				      Walk[i*2*faces*TIME + j*TIME + Times[k] - process_id*n_Walks/Rates[k] -i/Rates[k]][1],
				      Walk[i*2*faces*TIME + j*TIME + Times[k] - process_id*n_Walks/Rates[k] -i/Rates[k]][2],
				      n_bin_side, side_length);
	  
	  if (j%2 == 0) {SubTotal_cs[binOfWalker/faces][k]++;} else {SubTotal_cs[binOfWalker/faces][k]--;}
	  //if (j%2 == 0) {SubTotal[binOfWalker][k]++;} else {SubTotal[binOfWalker][k]--;}
	}
      }
    }
    
    else if (Rates[k] < 0) {
      for (i = 0; i < n_Walks && (process_id*Rates[k] + i*Rates[k]) <= Times[k]; i++) {
	for (j = 0; j < 2*faces; j++) {
	  binOfWalker = calculate_Bin(Walk[i*2*faces*TIME + j*TIME + Times[k] + process_id*Rates[k] + i*Rates[k]][0],
				      Walk[i*2*faces*TIME + j*TIME + Times[k] + process_id*Rates[k] + i*Rates[k]][1],
				      Walk[i*2*faces*TIME + j*TIME + Times[k] + process_id*Rates[k] + i*Rates[k]][2],
				      n_bin_side, side_length);
	  if (j%2 == 0) {SubTotal_cs[binOfWalker/faces][k]++;} else {SubTotal_cs[binOfWalker/faces][k]--;}
	  //if (j%2 == 0) {SubTotal[binOfWalker][k]++;} else {SubTotal[binOfWalker][k]--;}
	}
      }
    }
    
    else if (Rates[k] == 0) {
      for (i = 0; i < n_Walks; i++) {
	for (j = 0; j < 2*faces; j++) {
	  binOfWalker = calculate_Bin(Walk[i*2*faces*TIME + j*TIME + Times[k]][0],
				      Walk[i*2*faces*TIME + j*TIME + Times[k]][1],
				      Walk[i*2*faces*TIME + j*TIME + Times[k]][2],
				      n_bin_side, side_length);
	  if (j%2 == 0) {SubTotal_cs[binOfWalker/faces][k]++;} else {SubTotal_cs[binOfWalker/faces][k]--;}
	  //if (j%2 == 0) {SubTotal[binOfWalker][k]++;} else {SubTotal[binOfWalker][k]--;}
	}
      }
    }
    
  }

  // free memory
  free(Walk);
  
  
  int (*Total_cs)[len_Times] = NULL;
  if (process_id == 0) {
    Total_cs = malloc(sizeof(*Total_cs) * n_bin_side);
  }

  MPI_Reduce(SubTotal_cs, Total_cs, n_bin_side*len_Times, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  // free memory
  //free(SubTotal);
  free(SubTotal_cs);

  if (process_id == 0) {
    // print out cross section totals to a file
    for (k = 0; k < len_Times; k++) {
      FILE *g;                                                                                                      
      char name_g[FILENAME_MAX];
      snprintf(name_g, sizeof(name_g), "output2/Totalcs_%d.csv",k);
      g = fopen(name_g,"w");
      
      for (i = 0; i < n_bin_side; i++) {
	if (i != n_bin_side-1) {fprintf(g,"%d,",Total_cs[i][k]);}
	else {fprintf(g,"%d",Total_cs[i][k]);}
      }
      fclose(g);
    }
  }

  // free memory
  free(Total_cs);
  
  MPI_Finalize();
  return 0;
}

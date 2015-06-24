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
  int n_bin_side = 15; //per side 
  int TIME = 100;
  int DIM = 3;
  int faces = n_bin_side*n_bin_side;
  double interval = side_length / ((double)n_bin_side);

  int Times[] = {0,1,99};//{199,299,399};
  int Rates[] = {3,3,3};
  int len_Times = sizeof(Times)/sizeof(Times[0]);

  int n_Walks = 100;

  int i;
  int j;
  int k;
  int m;

  /*
   * A list of the initial x and z positions of the walkers. 
   * The list first proceeds across from x=0 to x=1 and then upwards from
   * z=0 to z=1.
   */
  double FacesXZ[n_bin_side*n_bin_side][2];
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
      if (process_id == 0) {
	fprintf(stderr, "Rate > 1: Need n_processes * n_Walks to be at least %d\n", TIME*Rates[len_Times-1]);
      }
      exit(EXIT_FAILURE);
    }
  } 
  else if (Rates[len_Times-1] < 0) {
    if (n_processes*n_Walks < TIME/abs(Rates[len_Times-1])) {
      if (process_id == 0) {
	fprintf(stderr, "Rate < 1: Need n_processes * n_Walks to be at least %d\n", TIME/abs(Rates[len_Times-1]));
      }
      exit(EXIT_FAILURE);
    }
  } 

  /*
   * Carbon Nanotube Lookup Table and RNG
   */
  FILE *Nanotube_File = NULL;
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
    if (process_id == 0) {
      fprintf(stderr,"Walk array could not be allocated\n");
    }
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

  /*
   * Gather the Random Walk data from all the processes to find 
   * the walker distribution at different heat flux rates and times
   * 
   * We give new jobs to some of the spawned processes.
   * Specifically, we want to collect the data for len_Times different 
   * times or rates.
   * So, we will reassign work to len_Times processes and stop all the 
   * other processes
   */  
  MPI_Group group_workers;
  MPI_Comm comm_workers;
  int ranks_workers[len_Times]; 
  for (i = 0; i < len_Times; i++) {
    ranks_workers[i] = i;
  }
  // create new group
  MPI_Group_incl(group_world, len_Times, ranks_workers, &group_workers);
  // create new communicator
  MPI_Comm_create(MPI_COMM_WORLD, group_workers, &comm_workers);
  
  // Gather from all processes to root, but Broadcast to working processes only
  double (*AllWalkers)[DIM] = malloc(sizeof(*AllWalkers) * TIME * n_Walks  * n_processes * 2 * faces);
  MPI_Gather(&Walk[0][0], TIME*DIM*n_Walks*2*faces, MPI_DOUBLE, &AllWalkers[0][0], TIME*DIM*n_Walks*2*faces,
	     MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int *Total_cs = NULL;
  int *Total = NULL;
  int binOfWalker;

  // jobs for the working processes
  if (process_id >= len_Times) {
    MPI_Finalize();
  }
  else {
    MPI_Bcast(AllWalkers, TIME*n_Walks*n_processes*DIM*2*faces, MPI_DOUBLE, 0, comm_workers);  

    // arrays that contain the walker distribution
    Total_cs = calloc(n_bin_side,sizeof(int));
    Total = calloc(n_bin_side*n_bin_side*n_bin_side, sizeof(int));

    /*
     * Rate exceeds 1 walker per time step
     * 
     * Processes each generate n_Walks for each bin on y=0 and y=1, i.e. for 2*faces bins.
     * So, we want to collect the data from TIME*Rates[process_id] random walk simulations,
     * meaning, we go determine which bins each of the TIME*Rates[process_id] walkers are 
     * in. 
     *
     * So, we iterate through the processes, then through the n_Walks on the processes, and 
     * finally through all the bins that walkers are coming from. 
     *
     * We have an offset m, which captures the idea that the walkers are injected at a constant
     * rate into the system. m is incremented for each group of Rates[process_id] random walks
     * we analyze and has an initial value of 0.
     */
    if (Rates[process_id] > 0) {
      m = 0;
      for (i = 0; i < n_processes && i*n_Walks < TIME*Rates[process_id]; i++) {
	for (j = 0; j < n_Walks && i*n_Walks+j < TIME*Rates[process_id]; j++) {
	  for (k = 0; k < 2*faces; k++) {
	    binOfWalker = 
	      calculate_Bin(AllWalkers[i*n_Walks*2*faces*TIME + j*TIME*2*faces + k*TIME + Times[process_id] - m][0],
			    AllWalkers[i*n_Walks*2*faces*TIME + j*TIME*2*faces + k*TIME + Times[process_id] - m][1],
			    AllWalkers[i*n_Walks*2*faces*TIME + j*TIME*2*faces + k*TIME + Times[process_id] - m][2],
			    n_bin_side, side_length);
	    // update cross section total
	    if (k%2) {Total[binOfWalker]++;} else {Total[binOfWalker]--;}   	
	    if (k%2) {Total_cs[binOfWalker/faces]++;} else {Total_cs[binOfWalker/faces]--;}
	  }
	  if ((j+1)%Rates[process_id]) {m++;}
	}
      }
    }

    /* 
     * Rate does not exceed 1 walker per time step
     *
     * Here Rates[process_id] < 0, meaning the rate of injection per bin is 
     * 1 / |Rates[process_id]|.
     * 
     * We iterate through as many random walk simulations similarly as the case
     * where Rate > 1.
     */
    else if (Rates[process_id] < 0) {
      m = 0;
      for (i = 0; i < n_processes && i*n_Walks < TIME/abs(Rates[process_id]); i++) {
	for (j = 0; j < n_Walks && i*n_Walks + j < TIME/abs(Rates[process_id]); j++) {
	  for (k = 0; k < 2*faces; k++) {
	    binOfWalker = 
	      calculate_Bin(AllWalkers[i*n_Walks*2*faces*TIME + j*TIME*2*faces + k*TIME + Times[process_id] + m*Rates[process_id]][0],
			    AllWalkers[i*n_Walks*2*faces*TIME + j*TIME*2*faces + k*TIME + Times[process_id] + m*Rates[process_id]][1],
			    AllWalkers[i*n_Walks*2*faces*TIME + j*TIME*2*faces + k*TIME + Times[process_id] + m*Rates[process_id]][2],
			    n_bin_side, side_length);
	    // update cross section total
	    if (k%2) {Total[binOfWalker]++;} else {Total[binOfWalker]--;}   	
	    if (k%2) {Total_cs[binOfWalker/faces]++;} else {Total_cs[binOfWalker/faces]--;}
	  }
	  m++;
	}
      }
    }
    
    /*
     * Rate is 0, simulate thermal equilibrium
     *
     * Here, we are placing n_processes * n_Walks walkers in the center
     * of each bin on the sides y=0 and y=1. Then, we look at the walker
     * distribution at Times[process_id].
     * 
     */
    else {
      for (i = 0; i < n_processes; i++) {
	for (j = 0; j < n_Walks; j++) {
	  for (k = 0; k < 2*faces; k++) {
	    binOfWalker = 
	      calculate_Bin(AllWalkers[i*TIME*n_Walks*2*faces + j*TIME*2*faces + k*TIME + Times[process_id]][0], 
			    AllWalkers[i*TIME*n_Walks*2*faces + j*TIME*2*faces + k*TIME + Times[process_id]][1],
			    AllWalkers[i*TIME*n_Walks*2*faces + j*TIME*2*faces + k*TIME + Times[process_id]][2],
			    n_bin_side, side_length);
	  
	    // update Total at bin accordingly
	    if (k%2) {Total[binOfWalker]++;} else {Total[binOfWalker]--;}   	
	    if (k%2) {Total_cs[binOfWalker/faces]++;} else {Total_cs[binOfWalker/faces]--;}
	  }
	}
      }
    }

    // print out cross section totals to a file
    FILE *g;
    char name_g[FILENAME_MAX];
    snprintf(name_g, sizeof(name_g), "run%d_cs.csv", process_id);
    g = fopen(name_g,"w");

    for (i = 0; i < n_bin_side; i++) {
      if (i != n_bin_side-1) {fprintf(g,"%d,",Total_cs[i]);}
      else {fprintf(g,"%d",Total_cs[i]);}
    }
    fclose(g);
    
    // print out the totals of all the bins to a file
    FILE *f;
    char name[FILENAME_MAX];
    snprintf(name, sizeof(name), "run%d_bins.csv", process_id);
    f = fopen(name, "w");

    for (i = 0; i < n_bin_side*n_bin_side*n_bin_side; i++) {
      if (i != n_bin_side*n_bin_side*n_bin_side-1) 
	fprintf(f, "%d,", Total[i]);
      else 
	fprintf(f, "%d", Total[i]);
    }
    
    fclose(f);
  }
  
  MPI_Finalize();
  return 0;
}

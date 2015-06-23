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
  int TIME = 100;
  int DIM = 3;
  int faces = n_bin_side*n_bin_side;
  double interval = side_length / ((double)n_bin_side);

  int Times[] = {1500,1700,1900};
  int Rates[] = {1,1,1};
  int len_Times = sizeof(Times)/sizeof(Times[0]);

  int i;
  int j;
  int k;

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
  double CenterFaces[n_bin_side];
  for (i = 0; i < n_bin_side; i++) {
    CenterFaces[i] = i*interval + interval/2;
  }

  int X[faces];
  int Z[faces];
  for (i = 0; i < faces; i++) {
    X[i] = i%n_bin_side;
    Z[i] = i/n_bin_side;
  }
  */

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
    if (n_processes < faces * 2 * TIME * Rates[len_Times-1]) {
      if (process_id == 0) {
	fprintf(stderr,"Rate > 1: Need at least %d processes for %d iterations, rate %d, and %d faces per side.\n", 
		faces * 2 * TIME * Rates[len_Times-1],TIME,Rates[len_Times-1],faces);
      }
      exit(EXIT_FAILURE);
    }
  } 
  else if (Rates[len_Times-1] < 0) {
    if (n_processes < faces * 2 * TIME / abs(Rates[len_Times-1])) {
      if (process_id == 0) {
	fprintf(stderr,"Rate < 1: Need at least %d processes for %d iterations, rate 1/%d and %d faces per side..\n", 
		faces * 2 * TIME / abs(Rates[len_Times-1]),TIME,abs(Rates[len_Times-1]),faces);
      }
      exit(EXIT_FAILURE);
    }
  } 
  else {
    if (n_processes % (faces * 2) != 0 ) {
      if (process_id == 0) {
	fprintf(stderr,"Equil: No. of processes must be divisible by %d.\n", faces * 2);
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
  double (*Walk)[DIM] = malloc(sizeof(*Walk) * TIME);
  if (Walk) {
    // process_id/2 is an int
    // This trick allows us to have a hot and cold walker
    // injected from opposite sides, i.e. (x,0,z) and (x,1,z).
    Walk[0][0] = FacesXZ[(process_id/2)%faces][0];
    Walk[0][1] = process_id % 2;
    Walk[0][2] = FacesXZ[(process_id/2)%faces][1];
  }
  else {
    if (process_id == 0) {
      fprintf(stderr,"Walk array could not be allocated\n");
    }
    exit(EXIT_FAILURE);
  }

  //printf("process %d is at <%g,%g,%g>\n",process_id, Walk[0][0], Walk[0][1], Walk[0][2]);

  /*
   * Generate a Random Walk for each process
   */
  iterate_Random_Walk(Walk, Nanotube_File, rng, TIME, side_length);
  
  /*
   * Synchronize processes with a barrier
   */
  MPI_Barrier(MPI_COMM_WORLD);  

  /*
   * Gather the Random Walk data from all the processes to find 
   * the walker distribution at different heat flux rates and times
   * 
   * We give new jobs to some of the spawned processes
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
  double (*AllWalkers)[DIM] = malloc(sizeof(*AllWalkers) * TIME * (n_processes));
  MPI_Gather(&Walk[0][0], TIME*DIM, MPI_DOUBLE, &AllWalkers[0][0], TIME*DIM,
	     MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int *Total_cs = NULL;
  int *Total = NULL;
  int binOfWalker;
  // jobs for the working processes
  if (process_id < len_Times) {
    MPI_Bcast(AllWalkers, TIME*n_processes*DIM, MPI_DOUBLE, 0, comm_workers);  

    Total_cs = calloc(n_bin_side,sizeof(int));
    Total = calloc(n_bin_side*n_bin_side*n_bin_side, sizeof(int));

    //
    // Each working process has a different rate and time 
    // 
    // Rate exceeds 1 walker per time step
    if (Rates[process_id] > 0) {
      for (i = 0, k = 0; i <= Times[process_id]*Rates[process_id]; i++) {
	for (j = 0; j < 2*faces; j++) {
	  binOfWalker = calculate_Bin(AllWalkers[i*TIME*2*faces + j*TIME + Times[process_id] - k][0],
				      AllWalkers[i*TIME*2*faces + j*TIME + Times[process_id] - k][1],
				      AllWalkers[i*TIME*2*faces + j*TIME + Times[process_id] - k][2],
				      n_bin_side, side_length);
	  // update Total at bin accordingly
	  if (j%2) {Total[binOfWalker]++;} else {Total[binOfWalker]--;}
	  if (j%2) {Total_cs[binOfWalker/faces]++;} else {Total_cs[binOfWalker/faces]--;}
	}
	if ((i+1)%Rates[process_id] == 0) {
	  k++;
	}
      } 
    }

    // Rate does not exceed 1 walker per time step
    else if (Rates[process_id] < 0) {	
      for (i = 0; i < Times[process_id] / abs(Rates[process_id]); i++) {
	for (j = 0; i < 2*faces; j++) {	  
	  binOfWalker = calculate_Bin(AllWalkers[i*TIME*2*faces + j*TIME + 
						 Times[process_id] + 
						 i*Rates[process_id]][0],
				      AllWalkers[i*TIME*2*faces + j*TIME + 
						 Times[process_id] + 
						 i*Rates[process_id]][1],
				      AllWalkers[i*TIME*2*faces + j*TIME + 
						 Times[process_id] + 
						 i*Rates[process_id]][2],
				      n_bin_side, side_length);
	  // update Total at bin accordingly
	  if (j%2) {Total[binOfWalker]++;} else {Total[binOfWalker]--;}	  
	  if (j%2) {Total_cs[binOfWalker/faces]++;} else {Total_cs[binOfWalker/faces]--;}
	}
      }
    }
    // Rate is 0, simulate thermal equilibrium
    else {
      for (i = 0; i < n_processes; i++) {
	binOfWalker = calculate_Bin(AllWalkers[i*TIME+Times[process_id]][0], 
				    AllWalkers[i*TIME+Times[process_id]][1],
				    AllWalkers[i*TIME+Times[process_id]][2],
				    n_bin_side, side_length);
	
	// update Total at bin accordingly
	if (i%2) {Total[binOfWalker]++;} else {Total[binOfWalker]--;}   	
	if (i%2) {Total_cs[binOfWalker/faces]++;} else {Total_cs[binOfWalker/faces]--;}
      }
    }

    FILE *g;
    char name_g[FILENAME_MAX];
    snprintf(name_g, sizeof(name_g), "run%d_cs.csv", process_id);
    g = fopen(name_g,"w");

    for (i = 0; i < n_bin_side; i++) {
      if (i != n_bin_side-1) {fprintf(g,"%d,",Total_cs[i]);}
      else {fprintf(g,"%d",Total_cs[i]);}
    }
    fclose(g);
    
    // print out total to a file
    FILE *f;
    char name[FILENAME_MAX];
    snprintf(name, sizeof(name), "%d.csv", process_id);
    f = fopen(name, "w");

    //operations to fill data into file i.txt;
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

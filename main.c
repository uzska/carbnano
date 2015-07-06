#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <mpi.h>
#include <gsl_rng.h>

#include "nanotube.h"
#include "walker.h"

int initWalkArray(double (*W)[3], double F[][2], int faces) {
  int i;
  for (i = 0; i < 2*faces; i++) {
    W[i][0] = F[i/2][0];
    W[i][1] = i%2;
    W[i][0] = F[i/2][1];
  }
  return 0;
}

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

int updateTime(int Time, int Rate, int ith_walk) {
  if (Rate > 0) {
    if (ith_walk % Rate == 0) {return (Time - 1);}
    else {return Time;}
  }
  else if (Rate < 0) {
    return (Time + Rate);
  }
  else {
    return Time;
  }
}

int main(int argc, char *argv[]) {

  /*
   * World Variables
   */
  double side_length = 1;
  int n_bin_side = 200; //per side 
  int faces = n_bin_side*n_bin_side;
  int DIM = 3;
  
  double interval = side_length / ((double)n_bin_side);
  int fineness = 5; //the much finer the nanotube grid is than the temp bins

  int Times[] = {200,400,600,800,1000,1200,1400,1600,1800,2000};
  int Rate = 4;
  int len_Times = sizeof(Times)/sizeof(Times[0]);

  int walks = 500; // number of different walks each process simulates
  
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
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&n_processes);  
  MPI_Comm_rank(MPI_COMM_WORLD,&process_id);


  /*
   * Error Checking, no of processes needed
   */
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

  /*
   * Total array 
   */
  int (*Total_cs)[len_Times] = NULL;
  if (process_id == 0) {
    Total_cs = malloc(sizeof(*Total_cs) * n_bin_side);
  }
  
  /*
   * Carbon Nanotube Lookup Table
   */
  char *Nanotubes = calloc(fineness * fineness * fineness * n_bin_side * n_bin_side * n_bin_side,
			   sizeof(char));
  int n_tubes;
  double i_Nanotubes[n_tubes][DIM];
  double f_Nanotubes[n_tubes][DIM];


  for (i = 0; i < n_tubes; i++) {
    double x = i_Nanotubes[i][0];    
    double y = i_Nanotubes[i][1];    
    double z = i_Nanotubes[i][2];    

    double dx = f_Nanotubes[i][0] - i_Nanotubes[i][0];
    double dy = f_Nanotubes[i][1] - i_Nanotubes[i][1];
    double dz = f_Nanotubes[i][2] - i_Nanotubes[i][2];
    
    double norm = sqrt(dx*dx + dy*dy + dz*dz);
    
    for (j = 0; j <= norm; j++) {
      int b = calculate_Bin(x,y,z,n_bin_side*fineness,1.0);
      Nanotubes[b] = 'x';

      x += dx/norm;
      y += dy/norm;
      z += dz/norm;
    }
  }


  /*
   * RNG
   */
  const gsl_rng_type *T = gsl_rng_taus;
  gsl_rng *rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL)+process_id);

  
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
    if (i) {initWalkArray(Walk, FacesXZ, faces);}

    int T = initTime(Times[0],Rate,process_id,walks,i);

    for (m = 0; m < len_Times && T >= 0; m++) {
      for (j = 0; j < 2*faces && T >= 0; j++) {
	for (k = 0; k < T; k++) {
	  Random_Walk(*(Walk+j),rng,Nanotubes,side_length,n_bin_side,fineness);
	}	
	bin = calculate_Bin(Walk[j][0], Walk[j][1], Walk[j][2], n_bin_side, side_length);
	if (j%2==0) {SubTotal_cs[bin/faces][m]++;} else {SubTotal_cs[bin/faces][m]--;} 
      }
      if (m!=len_Times-1) {T = Times[m+1] - Times[m];}
    }

  }

  /*
   * Synchronize processes with a barrier
   */
  MPI_Barrier(MPI_COMM_WORLD);  

  // free memory
  free(Walk);
    
  MPI_Reduce(SubTotal_cs, Total_cs, n_bin_side*len_Times, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
  // free memory
  free(SubTotal_cs);
  
  if (process_id == 0) {
    // print out cross section totals to a file
    for (k = 0; k < len_Times; k++) {
      FILE *g;                                                                                                      
      char name_g[FILENAME_MAX];
      snprintf(name_g, sizeof(name_g), "output3/Totalcs_%d.csv",k);
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

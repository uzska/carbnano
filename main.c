#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>
#include <gsl_cdf.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "walker.h"
/*
*/

int getMax(int x, int y, int z) {
  if (x >= y && x >= z) {return x;}
  else if (y >= x && y >= z) {return y;}
  else {return z;}
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
 * Main thread
 */
int main(int argc, char *argv[]) {

  int bins = 50;
  int faces = bins*bins;
  int fineness = 3;
  int Times[] = {100,200,300,400,500,600,700,800,900,1000,1100,1200};//{5,50,100,200,300};
  int Rate = 1;
  int walks = 600;
  int len_Times = sizeof(Times)/sizeof(Times[0]);

  int i,j,k,m;

  /* Initial Positions */
  double FacesXZ[faces][2];
  for (i = 0; i < bins; i++) {
    for (j = 0; j < bins; j++) {
      FacesXZ[i*bins + j][0] = 1.0/(2.0*bins) + (j*1.0)/(1.0*bins);
      FacesXZ[i*bins + j][1] = 1.0/(2.0*bins) + (i*1.0)/(1.0*bins);
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
      fprintf(stderr, 
	      "Rate > 1: Need n_processes * walks to be >= %d\n", 
	      Times[len_Times-1]*Rate);
      exit(EXIT_FAILURE);
    }
  } 
  else if (Rate < 0) {
    if (n_processes*walks < Times[len_Times-1]/abs(Rate)) {
      fprintf(stderr, 
	      "Rate < 1: Need n_processes * walks to be >= %d\n", 
	      Times[len_Times-1]/abs(Rate));
      exit(EXIT_FAILURE);
    }
  }


  /* RNG */
  const gsl_rng_type *T = gsl_rng_taus;
  gsl_rng *rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL)+process_id);

  /* Carbon Nanotube Lookup Table */
  int *Nanotubes = calloc(fineness * bins *
			  fineness * bins *
			  fineness * bins, sizeof(int));
  /* Nanotube variables  */
  int ntubes = 12;
  int tubeLen = 15;
  int tubeRad = 2;
  /* End of the Nanotubes */
  //double (*iNano)[3] = malloc(ntubes * sizeof(*iNano));
  //double (*fNano)[3] = malloc(ntubes * sizeof(*fNano));
  int (*iNano)[3] = malloc(ntubes * sizeof(*iNano));
  int (*fNano)[3] = malloc(ntubes * sizeof(*fNano));

  //getNanotubes(iNano,fNano,ntubes,tubeRad,tubeLen,fineness,bins,rng);

  /* manually set nanotube ends */
  if (ntubes) {    
iNano[0][0] = 122; iNano[0][1] = 62; iNano[0][2] = 112;
fNano[0][0] = 122; fNano[0][1] = 77; fNano[0][2] = 112;
iNano[1][0] = 86; iNano[1][1] = 127; iNano[1][2] = 139;
fNano[1][0] = 86; fNano[1][1] = 142; fNano[1][2] = 139;
iNano[2][0] = 21; iNano[2][1] = 104; iNano[2][2] = 77;
fNano[2][0] = 21; fNano[2][1] = 119; fNano[2][2] = 77;
iNano[3][0] = 27; iNano[3][1] = 101; iNano[3][2] = 69;
fNano[3][0] = 27; fNano[3][1] = 116; fNano[3][2] = 69;
iNano[4][0] = 24; iNano[4][1] = 104; iNano[4][2] = 28;
fNano[4][0] = 24; fNano[4][1] = 119; fNano[4][2] = 28;
iNano[5][0] = 39; iNano[5][1] = 117; iNano[5][2] = 125;
fNano[5][0] = 39; fNano[5][1] = 132; fNano[5][2] = 125;
iNano[6][0] = 99; iNano[6][1] = 97; iNano[6][2] = 55;
fNano[6][0] = 99; fNano[6][1] = 112; fNano[6][2] = 55;
iNano[7][0] = 33; iNano[7][1] = 13; iNano[7][2] = 42;
fNano[7][0] = 33; fNano[7][1] = 28; fNano[7][2] = 42;
iNano[8][0] = 119; iNano[8][1] = 81; iNano[8][2] = 25;
fNano[8][0] = 119; fNano[8][1] = 96; fNano[8][2] = 25;
iNano[9][0] = 50; iNano[9][1] = 45; iNano[9][2] = 138;
fNano[9][0] = 50; fNano[9][1] = 60; fNano[9][2] = 138;
iNano[10][0] = 71; iNano[10][1] = 74; iNano[10][2] = 69;
fNano[10][0] = 71; fNano[10][1] = 89; fNano[10][2] = 69;
iNano[11][0] = 9; iNano[11][1] = 96; iNano[11][2] = 70;
fNano[11][0] = 9; fNano[11][1] = 111; fNano[11][2] = 70;

  }
  /* manually set nanotube ends */

  for (i = 0; i < ntubes; i++) {
    // determine the x,y, and z bin numbers 
    /*
    int x = iNano[i][0] * bins * fineness;
    int y = iNano[i][1] * bins * fineness;
    int z = iNano[i][2] * bins * fineness;
    
    int xx = fNano[i][0] * bins * fineness;
    int yy = fNano[i][1] * bins * fineness;
    int zz = fNano[i][2] * bins * fineness;
    */
    double x = iNano[i][0];
    double y = iNano[i][1];
    double z = iNano[i][2];
    
    double xx = fNano[i][0];
    double yy = fNano[i][1];
    double zz = fNano[i][2];

    int max = getMax(xx-x,yy-y,zz-z);
    
    //printf("<%g,%g,%g> <%g,%g,%g>\n",x,y,z,xx,yy,zz);
    
    if (max <= 0) {      
      fprintf(stderr, "nanotube ends are the same");
      exit(EXIT_FAILURE);
    }
    

    double dx = (xx - x)/max;
    double dy = (yy - y)/max;
    double dz = (zz - z)/max;

    while (x <= xx && y <= yy && z <= zz) {
      // go through rad
      for (j = 0; j < tubeRad; j++) {
	for (k = 0; k < tubeRad; k++) {
	  for (m = 0; m < tubeRad; m++) {
	    int b = (x+j) + (y+k)*(bins*fineness)*(bins*fineness) +
	      (z+m)*(bins*fineness);
	    Nanotubes[b] = i+1;
	  }
	}
      }

      for (j = 0; j < tubeRad; j++) {
	for (k = 0; k < tubeRad; k++) {
	  for (m = 0; m < tubeRad; m++) {
	    int b = (x-j) + (y-k)*(bins*fineness)*(bins*fineness) + 
	      (z-m)*(bins*fineness);
	    Nanotubes[b] = i+1;
	  }
	}
      }      
      x += dx; y += dy; z += dz;
    }
  }

  /* Random Walk Array, holds only the current position  */  
  double (*Walk)[3] = malloc(sizeof(*Walk) * 2 * faces);
  /* Total walker Array */
  long long int (*Total)[len_Times] = NULL;
  if (process_id == 0) {
    Total = malloc(sizeof(*Total) * bins);
  }
  /* (sub)Total walker array held by each process */
  //long long int (*SubTotal)[len_Times] = calloc(bins,sizeof(*SubTotal));
  long long int (*SubTotal)[len_Times] = malloc(bins  * sizeof(*SubTotal));
  for (i = 0; i < len_Times; ++i) {
    for (j = 0; j < bins; ++j) {
      SubTotal[j][i] = 0LL;
    }
  }

  /* 
   * Random Walk algorithm--simutlates heat flux by offsetting
   * First move existing walkers in the system, 
   * then inject a new wave at the heated sides
   */
  int Bin;
  for (i = 0; i < walks; ++i) {
    initWalkArray(Walk,FacesXZ,faces);
    int T = initTime(Times[0],Rate,process_id,walks,i);    

    for (m = 0; m < len_Times; ++m) {
      for (j = 0; j < 2*faces && T >= 0; ++j) {
	for (k = 0; k < T; ++k) {
	  RandomWalk(*(Walk+j), rng, Nanotubes, iNano, fNano,
		     tubeLen, tubeRad, fineness, bins);
	}
	Bin = calculateBin(Walk[j][0], Walk[j][1], Walk[j][2], bins); 
	if (j%2==0) {SubTotal[Bin/faces][m] += 1LL;}
	else {SubTotal[Bin/faces][m] -= 1LL;}
      }
      if (m!=len_Times-1) {
	if (T >= 0) {T = Times[m+1] - Times[m];}
	else {
	  if (T + Times[m] >=0) {T += Times[m];}
	}
      }
    }
  }
  /* Spawn walkers, add to the last subtotal array  */
  if (process_id == 0) {
    SubTotal[0][len_Times-1] += (long long) (bins*bins);
    SubTotal[bins-1][len_Times-1] -= (long long) (bins*bins);
  }

  /* Synchronize processes with a barrier */
  MPI_Barrier(MPI_COMM_WORLD);  

  // free memory
  free(Walk);
    
  MPI_Reduce(SubTotal, Total, bins*len_Times,
	     MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
  // free memory
  free(SubTotal);

  if (process_id == 0) {
    for (m = 0; m < len_Times; ++m) {
      FILE *g;
      char name_g[FILENAME_MAX];
      snprintf(name_g, sizeof(name_g), "Total%d.csv", m);
      g = fopen(name_g,"w");
      
      for (i = 0; i < bins; ++i) {
	if (i != bins-1) {fprintf(g,"%lld,",Total[i][m]);}
	else {fprintf(g,"%lld",Total[i][m]);}
      }
      fclose(g);
    }
    // free memory
    free(Total);
  }

  MPI_Finalize();
  return 0;
}

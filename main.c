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
  int Times[] = {800};
  int Rate = 1;
  int walks = 400;
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
  int ntubes = 20;
  int tubeLen = 4;
  int tubeRad = 2;
  /* End of the Nanotubes */
  //double (*iNano)[3] = malloc(ntubes * sizeof(*iNano));
  //double (*fNano)[3] = malloc(ntubes * sizeof(*fNano));
  int (*iNano)[3] = malloc(ntubes * sizeof(*iNano));
  int (*fNano)[3] = malloc(ntubes * sizeof(*fNano));

  //getNanotubes(iNano,fNano,ntubes,tubeRad,tubeLen,fineness,bins,rng);

  /* manually set nanotube ends */
  if (ntubes) {    
    iNano[0][0] = 108; iNano[0][1] = 93; iNano[0][2] = 129;
    fNano[0][0] = 110; fNano[0][1] = 95; fNano[0][2] = 133;
    iNano[1][0] = 83; iNano[1][1] = 105; iNano[1][2] = 134;
    fNano[1][0] = 83; fNano[1][1] = 109; fNano[1][2] = 136;
    iNano[2][0] = 54; iNano[2][1] = 23; iNano[2][2] = 141;
    fNano[2][0] = 57; fNano[2][1] = 23; fNano[2][2] = 144;
    iNano[3][0] = 17; iNano[3][1] = 47; iNano[3][2] = 74;
    fNano[3][0] = 20; fNano[3][1] = 50; fNano[3][2] = 74;
    iNano[4][0] = 73; iNano[4][1] = 84; iNano[4][2] = 34;
    fNano[4][0] = 74; fNano[4][1] = 88; fNano[4][2] = 35;
    iNano[5][0] = 71; iNano[5][1] = 63; iNano[5][2] = 32;
    fNano[5][0] = 72; fNano[5][1] = 66; fNano[5][2] = 35;
    iNano[6][0] = 56; iNano[6][1] = 52; iNano[6][2] = 95;
    fNano[6][0] = 57; fNano[6][1] = 56; fNano[6][2] = 97;
    iNano[7][0] = 135; iNano[7][1] = 69; iNano[7][2] = 91;
    fNano[7][0] = 136; fNano[7][1] = 72; fNano[7][2] = 93;
    iNano[8][0] = 51; iNano[8][1] = 35; iNano[8][2] = 71;
    fNano[8][0] = 52; fNano[8][1] = 36; fNano[8][2] = 74;
    iNano[9][0] = 4; iNano[9][1] = 9; iNano[9][2] = 127;
    fNano[9][0] = 7; fNano[9][1] = 10; fNano[9][2] = 130;
    iNano[10][0] = 133; iNano[10][1] = 57; iNano[10][2] = 56;
    fNano[10][0] = 134; fNano[10][1] = 57; fNano[10][2] = 60;
    iNano[11][0] = 59; iNano[11][1] = 138; iNano[11][2] = 122;
    fNano[11][0] = 63; fNano[11][1] = 139; fNano[11][2] = 123;
    iNano[12][0] = 56; iNano[12][1] = 95; iNano[12][2] = 24;
    fNano[12][0] = 57; fNano[12][1] = 96; fNano[12][2] = 27;
    iNano[13][0] = 15; iNano[13][1] = 58; iNano[13][2] = 45;
    fNano[13][0] = 15; fNano[13][1] = 61; fNano[13][2] = 48;
    iNano[14][0] = 117; iNano[14][1] = 6; iNano[14][2] = 65;
    fNano[14][0] = 120; fNano[14][1] = 9; fNano[14][2] = 66;
    iNano[15][0] = 126; iNano[15][1] = 102; iNano[15][2] = 73;
    fNano[15][0] = 127; fNano[15][1] = 103; fNano[15][2] = 76;
    iNano[16][0] = 32; iNano[16][1] = 42; iNano[16][2] = 80;
    fNano[16][0] = 35; fNano[16][1] = 45; fNano[16][2] = 81;
    iNano[17][0] = 127; iNano[17][1] = 76; iNano[17][2] = 72;
    fNano[17][0] = 127; fNano[17][1] = 76; fNano[17][2] = 76;
    iNano[18][0] = 118; iNano[18][1] = 124; iNano[18][2] = 94;
    fNano[18][0] = 120; fNano[18][1] = 127; fNano[18][2] = 96;
    iNano[19][0] = 137; iNano[19][1] = 80; iNano[19][2] = 40;
    fNano[19][0] = 140; fNano[19][1] = 82; fNano[19][2] = 43;
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
      for (j = 1; j < tubeRad; j++) {
	for (k = 1; k < tubeRad; k++) {
	  for (m = 1; m < tubeRad; m++) {
	    int b = (x+j) + (y+k)*(bins*fineness)*(bins*fineness) +
	      (z+m)*(bins*fineness);
	    Nanotubes[b] = i+1;
	  }
	}
      }

      for (j = 1; j < tubeRad; j++) {
	for (k = 1; k < tubeRad; k++) {
	  for (m = 1; m < tubeRad; m++) {
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

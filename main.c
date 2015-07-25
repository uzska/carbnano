#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>
#include <gsl_cdf.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "walker.h"
/*
int outofBounds(double xi, double yi, double zi, 
		double xf, double yf, double zf, 
		int bins, int fineness, int rad) {

  double off = (rad*1.0) / (bins*fineness);
  if ((xi-off) < 0 || (xi + off) > 1.0) {return 1;}
  if ((yi-off) < 0 || (yi + off) > 1.0) {return 1;}
  if ((zi-off) < 0 || (zi + off) > 1.0) {return 1;}  

  if ((xf-off) < 0 || (xf + off) > 1.0) {return 1;}
  if ((yf-off) < 0 || (yf + off) > 1.0) {return 1;}
  if ((zf-off) < 0 || (zf + off) > 1.0) {return 1;}

  return 0;
}

int intersects(int x1, int y1, int z1,
	       int x2, int y2, int z2,
	       double Itubes[][3], double Ftubes[][3],
	       int tubes, int bins, int fineness, int rad) {
  int i;
  int intersect = 0;
  for (i = 0; i < tubes; i++) {
    if ((x1 >= (int)(Itubes[i][0]*bins*fineness - rad) && 
	x1 <= (int)(Ftubes[i][0]*bins*fineness + rad)) ||
	(x2 >= (int)(Itubes[i][0]*bins*fineness - rad) && 
	 x2 <= (int)(Ftubes[i][0]*bins*fineness + rad)) ) {
      intersect = 1;
    } else {return 0;}
    if ((y1 >= (int)(Itubes[i][1]*bins*fineness - rad) && 
	y1 <= (int)(Ftubes[i][1]*bins*fineness + rad)) ||
	(y2 >= (int)(Itubes[i][1]*bins*fineness - rad) && 
	 y2 <= (int)(Ftubes[i][1]*bins*fineness + rad)) ) {
      intersect = 1;
    } else {return 0;}
    if ((z1 >= (int)(Itubes[i][2]*bins*fineness - rad) && 
	z1 <= (int)(Ftubes[i][2]*bins*fineness + rad)) ||
	(z2 >= (int)(Itubes[i][2]*bins*fineness - rad) && 
	 z2 <= (int)(Ftubes[i][2]*bins*fineness + rad)) ) {
      intersect = 1;
    } else {return 0;}
  }
  return intersect;
}

int getNanotubes(double Itubes[][3], double Ftubes[][3], int ntubes, 
		 int rad, int len, int fineness, int bins, gsl_rng *rng) {
  int i; int j;
  double dist = (len*1.0)/(bins*fineness);

  double phi; 
  double theta;
  double x1,y1,z1;
  double x2,y2,z2;
  int xx1,yy1,zz1;
  int xx2,yy2,zz2;

  for (i = 0; i < ntubes; i++) {
    
    // get angles/orientations

    // random
    theta = gsl_rng_uniform(rng) * M_PI;
    phi = gsl_rng_uniform(rng) * 2*M_PI;

    // horiz
    //theta = M_PI/2;
    //phi = M_PI/2;

    //vertical
    //theta = 0;
    //phi = 0;

    // pick new random nanotube start locations
    x1 = gsl_rng_uniform(rng); xx1 = x1 * fineness * bins;
    y1 = gsl_rng_uniform(rng); yy1 = y1 * fineness * bins;
    z1 = gsl_rng_uniform(rng); zz1 = z1 * fineness * bins;

    // set final positions
    x2 = x1 + dist*sin(theta)*cos(phi); xx2 = x2 * fineness * bins;
    y2 = y1 + dist*sin(theta)*sin(phi); yy2 = y2 * fineness * bins;
    z2 = z1 + dist*cos(theta);          zz2 = z2 * fineness * bins;

    // bounds, reset nanotube positions if necessary
    while (outofBounds(x1,y1,z1,x2,y2,z2,bins,fineness,rad) ||
	   intersects(xx1,yy1,zz1,xx2,yy2,zz2,
		      Itubes,Ftubes,i,bins,fineness,rad)) {
      x1 = gsl_rng_uniform(rng);
      y1 = gsl_rng_uniform(rng);
      z1 = gsl_rng_uniform(rng);
      x2 = x1 + dist*sin(theta)*cos(phi);
      y2 = y1 + dist*sin(theta)*sin(phi);
      z2 = z1 + dist*cos(theta);
    }
    
    if (x1 < x2) {Itubes[i][0] = x1; Ftubes[i][0] = x2;}
    else {Ftubes[i][0] = x1; Itubes[i][0] = x2;}
    if (y1 < y2) {Itubes[i][1] = y1; Ftubes[i][1] = y2;}
    else {Ftubes[i][1] = y1; Itubes[i][1] = y2;}
    if (z1 < z2) {Itubes[i][2] = z1; Ftubes[i][2] = z2;}
    else {Ftubes[i][2] = z1; Itubes[i][2] = z2;}

  }
  return 0;
}
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
  int Times[] = {200,400,600,800,1000,1200,1400,1600,1800,2000};
  int Rate = 1;
  int walks = 125;
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
  double (*iNano)[3] = malloc(ntubes * sizeof(*iNano));
  double (*fNano)[3] = malloc(ntubes * sizeof(*fNano));
  //getNanotubes(iNano,fNano,ntubes,tubeRad,tubeLen,fineness,bins,rng);

  /* manually set nanotube ends */

  for (i = 0; i < ntubes; i++) {
    // determine the x,y, and z bin numbers
    int x = iNano[i][0] * bins * fineness;
    int y = iNano[i][1] * bins * fineness;
    int z = iNano[i][2] * bins * fineness;
    
    int xx = fNano[i][0] * bins * fineness;
    int yy = fNano[i][1] * bins * fineness;
    int zz = fNano[i][2] * bins * fineness;
    
    int max = getMax(xx-x,yy-y,zz-z);
    
    //printf("<%d,%d,%d> <%d,%d,%d>\n",x,y,z,xx,yy,zz);
    
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
  long long int (*SubTotal)[len_Times] = calloc(bins,sizeof(*SubTotal));

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
  }

  // free memory
  free(Total);
  
  MPI_Finalize();
  return 0;
}

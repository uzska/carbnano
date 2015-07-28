#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>
#include <math.h>

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
	       int Itubes[][3], int Ftubes[][3],
	       int tubes, int bins, int fineness, int rad) {
  int i;
  int intersect = 0;
  for (i = 0; i < tubes; i++) {
    /*
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
    */
    if ((x1 >= (int)(Itubes[i][0] - rad) && 
	x1 <= (int)(Ftubes[i][0] + rad)) ||
	(x2 >= (int)(Itubes[i][0] - rad) && 
	 x2 <= (int)(Ftubes[i][0] + rad)) ) {
      intersect = 1;
    } else {return 0;}
    if ((y1 >= (int)(Itubes[i][1] - rad) && 
	y1 <= (int)(Ftubes[i][1] + rad)) ||
	(y2 >= (int)(Itubes[i][1] - rad) && 
	 y2 <= (int)(Ftubes[i][1] + rad)) ) {
      intersect = 1;
    } else {return 0;}
    if ((z1 >= (int)(Itubes[i][2] - rad) && 
	z1 <= (int)(Ftubes[i][2] + rad)) ||
	(z2 >= (int)(Itubes[i][2] - rad) && 
	 z2 <= (int)(Ftubes[i][2] + rad)) ) {
      intersect = 1;
    } else {return 0;}

  }
  return intersect;
}

int getNanotubes(FILE *f, int Itubes[][3], int Ftubes[][3], int ntubes, 
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
    //theta = gsl_rng_uniform(rng) * M_PI;
    //phi = gsl_rng_uniform(rng) * 2*M_PI;

    // horiz
    theta = M_PI/2;
    phi = M_PI/2;

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
	   (getMax(abs(xx1-xx2),abs(yy1-yy2),abs(zz1-zz2)) == 0) ||
	   intersects(xx1,yy1,zz1,xx2,yy2,zz2,
		      Itubes,Ftubes,i,bins,fineness,rad)) {
      // set initial positions
      x1 = gsl_rng_uniform(rng); xx1 = x1 * fineness * bins;
      y1 = gsl_rng_uniform(rng); yy1 = y1 * fineness * bins;
      z1 = gsl_rng_uniform(rng); zz1 = z1 * fineness * bins;
      // set final positions
      x2 = x1 + dist*sin(theta)*cos(phi); xx2 = x2 * fineness * bins;
      y2 = y1 + dist*sin(theta)*sin(phi); yy2 = y2 * fineness * bins;
      z2 = z1 + dist*cos(theta);          zz2 = z2 * fineness * bins;
    }
    
    if (xx1 < xx2) {Itubes[i][0] = xx1; Ftubes[i][0] = xx2;}
    else {Ftubes[i][0] = xx1; Itubes[i][0] = xx2;}
    if (yy1 < yy2) {Itubes[i][1] = yy1; Ftubes[i][1] = yy2;}
    else {Ftubes[i][1] = yy1; Itubes[i][1] = yy2;}
    if (zz1 < zz2) {Itubes[i][2] = zz1; Ftubes[i][2] = zz2;}
    else {Ftubes[i][2] = zz1; Itubes[i][2] = zz2;}

    fprintf(f,"iNano[%d][0] = %d; iNano[%d][1] = %d; iNano[%d][2] = %d;\n",
	    i,Itubes[i][0],i,Itubes[i][1],i,Itubes[i][2]);
    fprintf(f,"fNano[%d][0] = %d; fNano[%d][1] = %d; fNano[%d][2] = %d;\n",
	    i,Ftubes[i][0],i,Ftubes[i][1],i,Ftubes[i][2]);

  }
  fclose(f);
  return 0;
}

int getMax(int x, int y, int z) {
  if (x >= y && x >= z) {return x;}
  else if (y >= x && y >= z) {return y;}
  else {return z;}
}


/*
 * Main thread
 */
int main(int argc, char *argv[]) {

  FILE *f = fopen("nanoEnds.txt","w");
  int bins = 50;
  int fineness = 3;

  const gsl_rng_type *T = gsl_rng_taus;
  gsl_rng *rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, time(NULL));  

  /* Nanotube variables  */
  int ntubes = 12;
  int tubeLen = 15;
  int tubeRad = 2;
  /* End of the Nanotubes */
  int (*iNano)[3] = malloc(ntubes * sizeof(*iNano));
  int (*fNano)[3] = malloc(ntubes * sizeof(*fNano));
  getNanotubes(f,iNano,fNano,ntubes,tubeRad,tubeLen,fineness,bins,rng);

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chardefs.h"
#include "mathdefs.h"
#include "misc.c"
#include "misc_math.c"
#include "fitsio.c"
//----------------------------------------------------------------------

#define TRUE      1                      // Define some handy constants  
#define FALSE     0                      // Define some handy constants  
#define fwhm      107.22                 // FWHM of the primary beam
#define SingleFieldNoise 0.28            // Single field centre noise in mJy/beam


float gasdev(float sigma);
int addgauss(float x, float y, float A);


int main(void)
{
  void exit(), readfits_map(), writefits_map();
  char infilename[20];
  int  i, j, k, n1, n2, n1center, n2center;
  float xmax, ymax;
  float mean_noise_level;
  float count, max, min, var;
  float xwidth, ywidth, area, num_pointings;
  float x, y;
  double weight;
  float cutoff, cellsize;
  float *data;
  char line[80];
  FILE  *infile, *pointingfile;
  header_param_list hpar;

  /*  read in the number of fields across the map, the field separation
      and the primary beam cutoff radius  */

  strcpy(infilename,"df.fits");

  readfits_map(infilename,&data,&hpar);
  srand(255);

  printf("\n Read fits file %s",infilename);
  printf("\n dimension = %d by %d \n", hpar.naxis[0],hpar.naxis[1]);

  n1 = hpar.naxis[0];
  n2 = hpar.naxis[1];


  // now calculate the matrix of noise values at each grid position

  n1center = n1/2;
  n2center = n2/2;
  data[n1center+n2center*n1] = 1.0;

 for(j=0;j<n2;j++) {
    for(i=0;i<n1;i++) {
 //     data[i+j*n1] = gasdev(data[i+j*n1]);
      printf("\n %f, %f",data[i+j*n1],gasdev(data[i+j*n1]));
    }
  }


  writefits_map ("out.fits", data, &hpar);
  free(data);



}

float gasdev(float sigma)
{
  static int iset=0;
  static float gset;
  float fac, rsq, v1, v2;

  do {
         v1 = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
	 v2 = 2.0*(rand()/(RAND_MAX+1.0)) - 1.0;
	 rsq = v1*v1 + v2*v2;
	 } while (rsq >= 1.0 || rsq == 0.0);
     fac = sqrt(-2.0*log(rsq)/rsq);
     return v2*fac*sigma;
}



int addgauss(float x, float y, float A)
{

//    return(exp(-2.772589*(offset/fwhm)*(offset/fwhm)));

}
















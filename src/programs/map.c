#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/chardefs.h"
#include "programs/mathdefs.h"
#include "programs/misc.c"
#include "programs/misc_math.c"
#include "programs/coord_utils.c"
#include "programs/fitsio.c"
//----------------------------------------------------------------------

int main(void)
{

  void exit(), readfits_map(), writefits_map();
  int  i, j, k, n1, n2;
  float xmax, ymax;
  int  N, m, numcenters;
  int rah, ram, decd, decm;
  float RAcen, DECcen, RArange, DECrange;
  float ramin = 6.03*15.0;
  float ramax = 8.23*15.0;
  float decmax = 12.2;
  float decmin = 10.8;
  float x, y;
  float RA, DEC, AST;
  float I, Q, U, V;
  float cellsize = 1.0;             // cell size of map in arc minutes
  float *dataI, *dataQ, *dataU, *dataV, *weight;  
  FILE  *infile;
  header_param_list hpar;





  cellsize = cellsize/60.0;                          //convert to degrees
  RAcen = (ramax + ramin)/2.0;
  DECcen = (decmax + decmin)/2.0;
  RArange =  ramax - ramin;
  DECrange = decmax - decmin;
  n1 = (int)(RArange/cellsize) + 1;
  n2 = (int)(DECrange/cellsize) + 1;

  printf("\n Map centre: %9.5f %8.4f, Map size: %5d x %5d \n", RAcen, DECcen, n1, n2);
  printf("Cell size: %7.3f\n", cellsize);

  dataI  = (float *) malloc (n1 * n2 * sizeof (float));
  dataQ  = (float *) malloc (n1 * n2 * sizeof (float)); 
  dataU  = (float *) malloc (n1 * n2 * sizeof (float)); 
  dataV  = (float *) malloc (n1 * n2 * sizeof (float));
  weight = (float *) malloc (n1 * n2 * sizeof (float));


  for(j=0;j<n2;j++) {
    for(i=0;i<n1;i++) {
      dataI[i+j*n1] = 0.0;
      dataQ[i+j*n1] = 0.0;
      dataU[i+j*n1] = 0.0;
      dataV[i+j*n1] = 0.0;
      weight[i+j*n1] = 0.0;
    }
  }

// ----------------first file ---------------------------------
  if( (infile = fopen("53103/fluxtime.dat","r") ) == NULL )
     { printf("can't open input file \n\r"); }

  printf("\n starting to read input file 53103");
  k = 0;  
  fscanf(infile,"%f %f %f %f %f %f %f", &RA, &DEC,&AST,&I,&Q,&U,&V);
  while(!feof(infile)) {
//    printf("\n %4d %8.5f %8.5f", k, RA, DEC);
      k++;
      RA = RA*15.0;
  
      x = n1-(RA-ramin)/cellsize;
      y = (DEC-decmin)/cellsize;

      i = (int)x;
      j = (int)y;
 //     printf("\n %4d %4d", i, j);
      if ( (i >= 0) && (i < n1) ) {
        if( (j >= 0) && (j < n2) ) {
          dataI[i+j*n1] = dataI[i+j*n1] + I;
          dataQ[i+j*n1] = dataQ[i+j*n1] + Q;          
          dataU[i+j*n1] = dataU[i+j*n1] + U;          
          dataV[i+j*n1] = dataV[i+j*n1] + V;
          weight[i+j*n1] = weight[i+j*n1] + 1.0;
	}
      } 
    fscanf(infile,"%f %f %f %f %f %f %f", &RA, &DEC,&AST,&I,&Q,&U,&V);
  } 

 printf("\n read %5d data points",k);
 fflush(stdout);

// ----------------second file ---------------------------------
 if( (infile = fopen("53105/fluxtime.dat","r") ) == NULL )
     { printf("can't open input file \n\r"); }

  printf("\n starting to read input file 53105");
  k = 0;  
  fscanf(infile,"%f %f %f %f %f %f %f", &RA, &DEC,&AST,&I,&Q,&U,&V);
  while(!feof(infile)) {
//    printf("\n %4d %8.5f %8.5f", k, RA, DEC);
      k++;
      RA = RA*15.0;
  
      x = n1-(RA-ramin)/cellsize;
      y = (DEC-decmin)/cellsize;

      i = (int)x;
      j = (int)y;
 //     printf("\n %4d %4d", i, j);
      if ( (i >= 0) && (i < n1) ) {
        if( (j >= 0) && (j < n2) ) {
          dataI[i+j*n1] = dataI[i+j*n1] + I;
          dataQ[i+j*n1] = dataQ[i+j*n1] + Q;          
          dataU[i+j*n1] = dataU[i+j*n1] + U;          
          dataV[i+j*n1] = dataV[i+j*n1] + V;
          weight[i+j*n1] = weight[i+j*n1] + 1.0;
	}
      } 
    fscanf(infile,"%f %f %f %f %f %f %f", &RA, &DEC,&AST,&I,&Q,&U,&V);
  } 

 printf("\n read %5d data points",k);
 fflush(stdout);
//---------------------------------------------------------

// ----------------third file ---------------------------------
 if( (infile = fopen("53107/fluxtime.dat","r") ) == NULL )
     { printf("can't open input file \n\r"); }

  printf("\n starting to read input file 53107");
  k = 0;  
  fscanf(infile,"%f %f %f %f %f %f %f", &RA, &DEC,&AST,&I,&Q,&U,&V);
  while(!feof(infile)) {
//    printf("\n %4d %8.5f %8.5f", k, RA, DEC);
      k++;
      RA = RA*15.0;
  
      x = n1-(RA-ramin)/cellsize;
      y = (DEC-decmin)/cellsize;

      i = (int)x;
      j = (int)y;
 //     printf("\n %4d %4d", i, j);
      if ( (i >= 0) && (i < n1) ) {
        if( (j >= 0) && (j < n2) ) {
          dataI[i+j*n1] = dataI[i+j*n1] + I;
          dataQ[i+j*n1] = dataQ[i+j*n1] + Q;          
          dataU[i+j*n1] = dataU[i+j*n1] + U;          
          dataV[i+j*n1] = dataV[i+j*n1] + V;
          weight[i+j*n1] = weight[i+j*n1] + 1.0;
	}
      } 
    fscanf(infile,"%f %f %f %f %f %f %f", &RA, &DEC,&AST,&I,&Q,&U,&V);
  } 

 printf("\n read %5d data points",k);
 fflush(stdout);
//---------------------------------------------------------


// ----------------fourth file ---------------------------------
 if( (infile = fopen("53108/fluxtime.dat","r") ) == NULL )
     { printf("can't open input file \n\r"); }

  printf("\n starting to read input file 53108");
  k = 0;  
  fscanf(infile,"%f %f %f %f %f %f %f", &RA, &DEC,&AST,&I,&Q,&U,&V);
  while(!feof(infile)) {
//    printf("\n %4d %8.5f %8.5f", k, RA, DEC);
      k++;
      RA = RA*15.0;
  
      x = n1-(RA-ramin)/cellsize;
      y = (DEC-decmin)/cellsize;

      i = (int)x;
      j = (int)y;
 //     printf("\n %4d %4d", i, j);
      if ( (i >= 0) && (i < n1) ) {
        if( (j >= 0) && (j < n2) ) {
          dataI[i+j*n1] = dataI[i+j*n1] + I;
          dataQ[i+j*n1] = dataQ[i+j*n1] + Q;          
          dataU[i+j*n1] = dataU[i+j*n1] + U;          
          dataV[i+j*n1] = dataV[i+j*n1] + V;
          weight[i+j*n1] = weight[i+j*n1] + 1.0;
	}
      } 
    fscanf(infile,"%f %f %f %f %f %f %f", &RA, &DEC,&AST,&I,&Q,&U,&V);
  } 

 printf("\n read %5d data points",k);
 fflush(stdout);
//---------------------------------------------------------



//------------------- divide sums by weights ----------------------
  for(j=0;j<n2;j++) {
    for(i=0;i<n1;i++) {
      if (weight[i+j*n1] > 0.0) {
        dataI[i+j*n1] = dataI[i+j*n1]/weight[i+j*n1];
        dataQ[i+j*n1] = dataQ[i+j*n1]/weight[i+j*n1];
        dataU[i+j*n1] = dataU[i+j*n1]/weight[i+j*n1];
        dataV[i+j*n1] = dataV[i+j*n1]/weight[i+j*n1];
      }
    }
  }


  fclose(infile); 

  // Write the noisemap as a fits file

// #if 0
  init_header_param_list (&hpar);  /* initialize parameter records */
  hpar.bitpix = -32;
  hpar.num_axes = 2;
  hpar.naxis[0] = n1;
  hpar.naxis[1] = n2;
  sprintf (hpar.ctype[0], "RA---TAN");
  sprintf (hpar.ctype[1], "DEC---TAN");
  hpar.crval[0] = RAcen;		/* hours */
  hpar.crval[1] = DECcen;               /* degrees */
  hpar.crpix[0] = 0.5 + n1 / 2.0;	/* image center in pixels */
  hpar.crpix[1] = 0.5 + n2 / 2.0;
  hpar.cdelt[0] = -cellsize;                     /* degrees */
  hpar.cdelt[1] = cellsize;                     /* degrees */
  sprintf (hpar.bunit, "Kelvin");
  sprintf (hpar.object, "GALFACTS Test Region");
  sprintf (hpar.telescope, "Arecibo");
  writefits_map ("I.fits", dataI, &hpar);
  writefits_map ("Q.fits", dataQ, &hpar);
  writefits_map ("U.fits", dataU, &hpar);
  writefits_map ("V.fits", dataV, &hpar);
  free(dataI); 
  free(dataQ);
  free(dataU);
  free(dataV);
// #endif
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chardefs.h"
#include "mathdefs.h"
#include "misc.c"
#include "misc_math.c"
#include "coord_utils.c"
#include "fitsio.c"
//----------------------------------------------------------------------

#define TRUE      1                      // Define some handy constants  
#define FALSE     0                      // Define some handy constants  
#define fwhm      107.22                 // FWHM of the primary beam
#define SingleFieldNoise 0.19            // Single field centre noise in mJy/beam

double primary_beam(float offset);
void degtohms(float ra0, int *h_ptr, int *m_ptr , float *s_ptr);
void degtodms(float dec0, int *d_ptr, int *m_ptr, float *s_ptr);
double arclength(float ra1, float dec1, float ra2, float dec2);


int main(void)
{

  void exit(), readfits_map(), writefits_map();
  int  i, j, k, n1, n2;
  float xmax, ymax;
  int  N, m, numcenters;
  int rah, ram, decd, decm;
  float ras, decs;
  float xc[1000], yc[1000], x, y, distance, mean_noise_level, rms_noise_var;
  float xgrid[1000], ygrid[1000];
  float count, max, min, var;
  float RA, DEC, dra, ddec, rowDEC;
  float latitude, longitude;
  double xrad, yrad ;    // RA and DEC or positions in radians for input to conversion routines
  double lat, lon;
  float diameter, area;
  double dum1, dum2, dum3;
  double weight;
  float cutoff, s, cellsize;
  int numlayers;
  float *data;  
  char line[80];
  FILE  *infile, *pointingfile, *galannfile;
  header_param_list hpar;

  /*  read in the number of fields across the map, the field separation
      and the primary beam cutoff radius  */

  if( (infile = fopen("deepfield.par","r") ) == NULL )
     { printf("can't open input file \n\r"); }

  if( (pointingfile = fopen("deepfieldcel.ann","w") ) == NULL )
     { printf("can't open output file \n\r"); }

 if( (galannfile = fopen("deepfieldgal.ann","w") ) == NULL )
     { printf("can't open output file \n\r"); }


  fgets(line,80,infile);
  sscanf(line,"%d",&numlayers);         // number of layers of hexagons
  fgets(line,80,infile);
  sscanf(line,"%f",&s);                 // separation of field centres in arcminutes
  fgets(line,80,infile);
  sscanf(line,"%f",&cutoff);            // cut off distance for primary beam
  fgets(line,80,infile);
  sscanf(line,"%f",&cellsize);          // cellsize in arcminutes
  fgets(line,80,infile);
  sscanf(line,"%f %f",&RA, &DEC);       // RA and DEC of field centre

  RA = RA * 15.0;                 // convert RA to decimal degrees
 
  fprintf(pointingfile,"COLOUR GREEN \n");
  fprintf(galannfile,"COLOUR RED \n");

  printf("\n numlayers= %d, separation = %7.2f, cutoff= %5.1f, cellsize = %5.1f",
	 numlayers, s, cutoff, cellsize);

  // num_centers is the number of pointing centres within the square region, 
  // assuming a close-packed hexagonal grid of fields separated by separation.
  // The region has num_fields along each edge.  The example below is for 
  // num_fields = 4;  The total num_centers is then 16 + 9 = 25; The spacing
  // of fields along the vertical edge (y axis) is sqrt(3)*s and the  
  // spacing along the horizontal edge is equal to s.

  //           
  //        X   X   X
  //      X   X   X   X
  //     X  X   X   X   X 
  //      X   X   X   X
  //        X   X   X
  //    

// For n layers (i.e. n hexagons) there are N=2*(n-1)+1 fields across the centre row
// And there will be N rows, n-1 above the centre row and n-1 below. If k is the
// index of the row above or below the centre (k = 1, to n-1), then the kth row has
// N-k fields.           

// note for a field spacing is 's' the spacing between rows is sqrt(3/4)*s, and 
// for row k the first field position is shifted by k*s/2 in the x direction. 
 
  N = 2*(numlayers-1)+1;                         // number of fields in centre row 
  numcenters = 3*(numlayers-1)*numlayers + 1;    // total number of fields
  diameter = N*s/60.0;
  area = (PI/4.0)*(diameter)*(diameter); 
 
  printf("\n diameter = %7.1f deg, Deep Field Area = %7.2f sq deg. \n", 
           diameter, area);
  printf("\n number of rows in centre row = %d",N);
  printf("\n number of field centres in map  = %d \n",numcenters);
 
//  Load array of coordinates of the field centers (xc,yc).  These coordinates
//  are in grid units.  The scale of the grids on the sky is given by cellsize.
   
  m = 0;                                 // index counter for fields
  ddec = s/60.0*sqrt(3.0/4.0);               // vertical (dec) distance between rows
  dra = (s/60.0)*(1.0/15.0*cos(DEC*PI/180.0));  // separation in RA

// fill centre row
  for(j=0;j<N;j++) {
   dra = (s/60.0)*(1.0/cos(DEC*PI/180.0));
   xc[m] = (float)(j-N/2)*dra + RA;
   yc[m] = DEC;
   xrad = xc[m]*D2RAD;
   yrad = yc[m]*D2RAD;
   eq2000_to_gal(xrad, yrad, &lon, &lat);
   latitude = lat*RAD2D;
   longitude = lon*RAD2D;
   degtohms(xc[m],&rah,&ram,&ras);
   degtodms(yc[m],&decd,&decm,&decs);
   printf("\n xc[%d], yc[%d] = %3d %02d %4.1f, %3d %02d %4.1f",m,m,rah,ram,ras,decd,decm,decs);
//   printf("\n long = %6.2f, lat = %6.2f",longitude, latitude);
   fprintf(pointingfile,"CROSS W %7.2f %7.2f 0.1  0.1  \n",xc[m], yc[m]);
   fprintf(galannfile,"CROSS W %7.2f %7.2f 0.2  0.2  \n",longitude, latitude);
   m++;
   }
  printf("\n");

// rows above the centre
  for(k=1;k<=numlayers-1;k++) {
    for(j=0; j<N-k; j++) { 
      yc[m] = DEC + (float)k * ddec;
      dra = (s/60.0)*(1.0/cos(yc[m]*PI/180.0));      
      xc[m] = RA + (float)(j-N/2)*dra + (float)k*dra/2.0;
      xrad = xc[m]*D2RAD;
      yrad = yc[m]*D2RAD;
      eq2000_to_gal(xrad, yrad, &lon, &lat);
      latitude = lat*RAD2D;
      longitude = lon*RAD2D;
      fprintf(pointingfile,"CROSS W %7.2f %7.2f 0.04  0.04  \n",xc[m], yc[m]);
      degtohms(xc[m],&rah,&ram,&ras);
      degtodms(yc[m],&decd,&decm,&decs);
      printf("\n xc[%d], yc[%d] = %3d %02d %4.1f, %3d %02d %4.1f",m,m,rah,ram,ras,decd,decm,decs);
      fprintf(galannfile,"CROSS W %7.2f %7.2f 0.2  0.2  \n",longitude, latitude);
      m++;
    }
    printf("\n");
   }

// rows below the centre
  for(k=1;k<=numlayers-1;k++) {
    for(j=0; j<N-k; j++) { 
      yc[m] = DEC -(float)k * ddec; 
      dra = (s/60.0)*(1.0/cos(yc[m]*PI/180.0));        
      xc[m] = RA + (float)(j-N/2)*dra + (float)k*dra/2.0;
      xrad = xc[m]*D2RAD;
      yrad = yc[m]*D2RAD;
      eq2000_to_gal(xrad, yrad, &lon, &lat);
      latitude = lat*RAD2D;
      longitude = lon*RAD2D;
      fprintf(pointingfile,"CROSS W %7.2f %7.2f 0.04  0.04  \n",xc[m], yc[m]);
      degtohms(xc[m],&rah,&ram,&ras);
      degtodms(yc[m],&decd,&decm,&decs);
      printf("\n xc[%d], yc[%d] = %3d %02d %4.1f, %3d %02d %4.1f",m,m,rah,ram,ras,decd,decm,decs);
      fprintf(galannfile,"CROSS W %7.2f %7.2f 0.2  0.2  \n",longitude, latitude);
      m++;
    }
  printf("\n");
  }

 fclose(pointingfile);
 printf("\n total number of fields = %3d \n", m);

  /*  Noise statistics are calculate for the region within the outer ring of 
      field centres. Therefore the lower left hand corner of the map is set 
      at the pointing centre of the field at lower left. */

  ymax = 1.2*diameter*60.0/cellsize;
  xmax = ymax;
//  xmax = ymax*(1.0/cos(DEC*PI/180.0));
   // field dimensions in pixels 
  n1 = xmax-1;
  n2 = ymax-1;
  data = (float *) malloc (n1 * n2 * sizeof (float));
  printf("\n pixel size  = %5.1f",cellsize);  
 
  printf("\n field width = %d, field height = %d pixels",n1, n2);
  
 // now calculate the matrix of noise values at each grid position 

  mean_noise_level = 0.0; max = 0.0; min = 10.0; count = 0.0;
  for(j=0;j<n2;j++) {
    for(i=0;i<n1;i++) {
      data[i+j*n1] = 0.0;
      x = RA-(((float)n1/2.0-(float)i)*cellsize/60);
      y = DEC-(((float)n2/2.0-(float)j)*cellsize/60);
	for (k=0;k<numcenters;k++) {
	  distance = 60.0*arclength(x,y, xc[k], yc[k]);
	  if(distance < cutoff) {
	    weight = primary_beam(distance);
	    data[i+j*n1] = data[i+j*n1] + weight*weight;
	  }
	}
      data[i+j*n1] = SingleFieldNoise/sqrt(data[i+j*n1]);
      count++;
      if(data[i+j*n1] > max) max = data[i+j*n1];
      if(data[i+j*n1] < min) min = data[i+j*n1];
      mean_noise_level = mean_noise_level + data[i+j*n1];
    }
  }
  mean_noise_level = mean_noise_level/count;
  var =( max - min) / mean_noise_level;
  printf("\n  mean noise level = %f\n", mean_noise_level);
  printf("\n  max = %f, min = %f, fractional range = %f \n", max, min, var );

  // now calculate rms of noise level about the mean 

  rms_noise_var = 0.0;
  count = 0.0;
  for(j=0;j<n2;j++) {
    for(i=0;i<n1;i++) {
      rms_noise_var = rms_noise_var + 
         (mean_noise_level - data[i+j*n1]) * (mean_noise_level - data[i+j*n1]);
      count++;
    }
  }
   rms_noise_var = sqrt(rms_noise_var/count);
  printf("\n  rms noise variation = %f\n\n", rms_noise_var);

  // Write the noisemap as a fits file

  init_header_param_list (&hpar);  /* initialize parameter records */
  hpar.bitpix = -32;
  hpar.num_axes = 2;
  hpar.naxis[0] = n1;
  hpar.naxis[1] = n2;
  sprintf (hpar.ctype[0], "Right Ascension");
  sprintf (hpar.ctype[1], "Declination");
  hpar.crval[0] = RA;			/* degrees */
  hpar.crval[1] = DEC;                  /* degrees */
  hpar.crpix[0] = 0.5 + n1 / 2.0;	/* image center in pixels */
  hpar.crpix[1] = 0.5 + n2 / 2.0;
  hpar.cdelt[0] = cellsize/60.0;                     /* degrees */
  hpar.cdelt[1] = cellsize/60.0;                     /* degrees */
  sprintf (hpar.bunit, "mJy/beam");
  sprintf (hpar.object, "DRAO Deep Field Simulation");
  sprintf (hpar.telescope, "DRAO-ST");
  writefits_map ("df.fits", data, &hpar);
  free(data);

}

void degtohms(float ra0, int *h_ptr, int *m_ptr, float *s_ptr)
{
  // converts input right ascension, ra0, in decimal degrees
  // to hours, minutes and seconds.  

  int  ra_h, ra_m;
  float rahours, ramin, rasec;

  rahours = ra0/15.0;
  ra_h = rahours;
  *h_ptr = ra_h;
  ramin = 60.0*(rahours - ra_h);
  ra_m = ramin;
  *m_ptr = ra_m;
  rasec = 60.0*(ramin - ra_m);
   *s_ptr = rasec;

}
  
void degtodms(float dec0, int *d_ptr, int *m_ptr, float *s_ptr)
{
  // converts input right ascension, ra0, in decimal degrees
  // to hours, minutes and seconds.  

  int  dec_d, dec_m;
  float decdeg, decmin, decsec;

  decdeg = dec0;
  dec_d = decdeg;
  *d_ptr = dec_d;
  decmin = 60.0*(decdeg - dec_d);
  dec_m = decmin;
  *m_ptr = dec_m;
  decsec = 60.0*(decmin - dec_m);
   *s_ptr = decsec;

}

double arclength(float ra1, float dec1, float ra2, float dec2)
 {
  double sep, sep_approx;
  double d1,d2,r1,r2;
  double convrt = PI/180;
 // convert angles from degrees to radians
  d1 = dec1* convrt;
  r1 = ra1 * convrt;
  d2 = dec2* convrt;
  r2 = ra2 * convrt;

  sep_approx = sqrt (cos(d1)*cos(d2)*(r2-r1)*(r2-r1) + (d2-d1)*(d2-d1));
  if (fabs (cos (sep_approx)) > cos (MAX_APPROX_ANGLE)) {
    /** use approximation -- angle very close to 0 or 180 **/
    sep = sep_approx;
  } else {
    /** use full forumlation **/
    sep = acos (sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(r2-r1));
  }
  sep = sep/convrt;	// convert back to degrees
  return sep;
}


double primary_beam(float offset)
{

    return(exp(-2.772589*(offset/fwhm)*(offset/fwhm)));

}

















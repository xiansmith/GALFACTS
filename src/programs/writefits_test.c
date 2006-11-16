/* WRITEFITS_TEST.C */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "/u/gibson/src/util/chardefs.h"
#include "/u/gibson/src/util/mathdefs.h"
#include "/u/gibson/src/util/misc.c"
#include "/u/gibson/src/util/misc_math.c"
#include "/u/gibson/src/fits/fitsio.c"

/*****************************************************************************
 **									    **
 **				M A I N					    **
 **									    **
 *****************************************************************************/
main (argc, argv)
  int argc;		/* # command line args, INCLUDING command itself */
  char *argv[];		/* arg strings; argv[0] is command string */
{
  void exit(), readfits_map(), writefits_map();
  int i, j, n1, n2;
  float *data;
  header_param_list hpar;

  n1 = 100;
  n2 = 100;
  data = (float *) malloc (n1 * n2 * sizeof (float));
  for (j=0; j<n2; j++) {
    for (i=0; i<n1; i++) {
      data[i+j*n1] = 1.0*i;
    }
  }
  init_header_param_list (&hpar);  /* initialize parameter records */
  hpar.bitpix = -32;
  hpar.num_axes = 2;
  hpar.naxis[0] = n1;
  hpar.naxis[1] = n2;
  sprintf (hpar.ctype[0], "X offset (arcmin)");
  sprintf (hpar.ctype[1], "Y offset (arcmin)");
  hpar.crval[0] = 0.0;			/* arcmin */
  hpar.crval[1] = 0.0; 
  hpar.crpix[0] = 0.5 + n1 / 2.0;	/* image center in pixels */
  hpar.crpix[1] = 0.5 + n2 / 2.0;
  hpar.cdelt[0] = 1.0;			/* arcmin */
  hpar.cdelt[1] = 1.0;
  sprintf (hpar.bunit, "mJy/beam");
  sprintf (hpar.object, "DRAO Deep Field Simulation");
  sprintf (hpar.telescope, "DRAO-ST");
  writefits_map ("writefits_test.fits", data, &hpar);
  free (data);
}

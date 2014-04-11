#include "io.h"
#include "fitsio.h"
#include "mathdefs.h"
#include "misc_math.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>



/*****************************************************************************
 **									    **
 **		P R I N T H E X   --  output a string in hex format	    **
 **									    **
 *****************************************************************************/
void printhex (startaddr, numbytes, endline)
    char *startaddr;
    int numbytes, endline;
{
    int k;

    for (k=0; k<numbytes; k++) {
        if ( (0xff & (int)startaddr[k]) < 0x10) printf ("0");
	printf ("%X", (0xff & (int)startaddr[k]) );
    }
    if (endline==1) printf ("\n");
}

/*****************************************************************************
 **									    **
 **		I N I T _ S U B S E T _ P A R A M _ L I S T		    **
 **									    **
 *****************************************************************************/
void init_subset_param_list (spar_ptr)
  subset_param_list *spar_ptr;
{
  int i;
  subset_param_list spar;

  spar = *spar_ptr;

  spar.num_axes = 0;
  for (i=0; i<MAX_NUM_AXES; i++) {
    spar.start[i]  = 1;
    spar.stop[i]   = 1;
    spar.stride[i] = 1;
  }

  *spar_ptr = spar;
}

/*****************************************************************************
 **									    **
 **		I N I T _ H E A D E R _ P A R A M _ L I S T		    **
 **									    **
 *****************************************************************************/
void init_header_param_list (hpar_ptr)
  header_param_list *hpar_ptr;
{
  int i, j;

  (*hpar_ptr).simple = TRUE;
  (*hpar_ptr).blocked = FALSE;
  (*hpar_ptr).extend = FALSE;
  (*hpar_ptr).bitpix = 0;
  (*hpar_ptr).num_axes = 0;
  for (i=0; i<MAX_NUM_AXES; i++) {
    (*hpar_ptr).naxis[i] = 1;
    (*hpar_ptr).crval[i] = 0.0;
    (*hpar_ptr).crpix[i] = 0.0;
    (*hpar_ptr).cdelt[i] = 0.0;
    (*hpar_ptr).crota[i] = 0.0;
    for (j=0; j<LINELEN; j++) (*hpar_ptr).ctype[i][j] = 0;
  }
  (*hpar_ptr).bscale = 1.0;
  (*hpar_ptr).bzero = 0.0;
  (*hpar_ptr).blank = -32768;
  (*hpar_ptr).datamin = 0.0;
  (*hpar_ptr).datamax = 0.0;
  (*hpar_ptr).obsfreq = 0.0;
  (*hpar_ptr).exptime = 0.0;
  (*hpar_ptr).epoch = 0.0;
  (*hpar_ptr).equinox = 0.0;
  (*hpar_ptr).ncomment = 0;
  (*hpar_ptr).nhistory = 0;
  for (i=0; i<LINELEN; i++) {
    (*hpar_ptr).bunit[i] = 0;
    (*hpar_ptr).date[i] = 0;
    (*hpar_ptr).object[i] = 0;
    (*hpar_ptr).origin[i] = 0;
    (*hpar_ptr).telescope[i] = 0;
    (*hpar_ptr).instrument[i] = 0;
    (*hpar_ptr).observer[i] = 0;
    (*hpar_ptr).filename[i] = 0;
  }
  for (i=0; i<MAX_COMMENT; i++) 
    for (j=0; j<LINELEN; j++) (*hpar_ptr).comment[i][j] = 0;
  for (i=0; i<MAX_HISTORY; i++) 
    for (j=0; j<LINELEN; j++) (*hpar_ptr).history[i][j] = 0;
}

/*****************************************************************************
 **									    **
 **		C O P Y _ H E A D E R _ P A R A M _ L I S T		    **
 **									    **
 *****************************************************************************/
void copy_header_param_list (hpar2_ptr, hpar1_ptr)
  header_param_list *hpar2_ptr, *hpar1_ptr;
{
  int i, j;

  (*hpar2_ptr).simple   = (*hpar1_ptr).simple;
  (*hpar2_ptr).blocked  = (*hpar1_ptr).blocked;
  (*hpar2_ptr).extend   = (*hpar1_ptr).extend;
  (*hpar2_ptr).bitpix   = (*hpar1_ptr).bitpix;
  (*hpar2_ptr).num_axes = (*hpar1_ptr).num_axes;
  for (i=0; i<MAX_NUM_AXES; i++) {
    for (j=0; j<LINELEN; j++) 
	(*hpar2_ptr).ctype[i][j] = (*hpar1_ptr).ctype[i][j];
    (*hpar2_ptr).naxis[i]    = (*hpar1_ptr).naxis[i];
    (*hpar2_ptr).crval[i]    = (*hpar1_ptr).crval[i];
    (*hpar2_ptr).crpix[i]    = (*hpar1_ptr).crpix[i];
    (*hpar2_ptr).cdelt[i]    = (*hpar1_ptr).cdelt[i];
    (*hpar2_ptr).crota[i]    = (*hpar1_ptr).crota[i];
  }
  (*hpar2_ptr).bscale        = (*hpar1_ptr).bscale;
  (*hpar2_ptr).bzero         = (*hpar1_ptr).bzero;
  (*hpar2_ptr).blank         = (*hpar1_ptr).blank;
  (*hpar2_ptr).datamin       = (*hpar1_ptr).datamin;
  (*hpar2_ptr).datamax       = (*hpar1_ptr).datamax;
  (*hpar2_ptr).obsfreq       = (*hpar1_ptr).obsfreq;
  (*hpar2_ptr).exptime       = (*hpar1_ptr).exptime;
  (*hpar2_ptr).epoch         = (*hpar1_ptr).epoch;
  (*hpar2_ptr).equinox       = (*hpar1_ptr).equinox;
  (*hpar2_ptr).ncomment      = (*hpar1_ptr).ncomment;
  (*hpar2_ptr).nhistory      = (*hpar1_ptr).nhistory;
  for (j=0; j<LINELEN; j++) {
    (*hpar2_ptr).bunit[j]      = (*hpar1_ptr).bunit[j];
    (*hpar2_ptr).date[j]       = (*hpar1_ptr).date[j];
    (*hpar2_ptr).object[j]     = (*hpar1_ptr).object[j];
    (*hpar2_ptr).origin[j]     = (*hpar1_ptr).origin[j];
    (*hpar2_ptr).telescope[j]  = (*hpar1_ptr).telescope[j];
    (*hpar2_ptr).instrument[j] = (*hpar1_ptr).instrument[j];
    (*hpar2_ptr).observer[j]   = (*hpar1_ptr).observer[j];
    (*hpar2_ptr).filename[j]   = (*hpar1_ptr).filename[j];
    for (i=0; i<MAX_COMMENT; i++) 
	(*hpar2_ptr).comment[i][j] = (*hpar1_ptr).comment[i][j];
    for (i=0; i<MAX_HISTORY; i++) 
	(*hpar2_ptr).history[i][j] = (*hpar1_ptr).history[i][j];
  }
}

/*****************************************************************************
 **									    **
 **		W C S _ F R O M _ P I X                      		    **
 **									    **
 **	Given a header parameter structure, axis index, and pixel	    **
 **	coordinate, return the world coordinate value for that position.    **
 **	The axis and pixel coord both start counting at 1, not 0!	    **
 **									    **
 **	NOTE: this assumes no projection transformation is required!	    **
 **									    **
 *****************************************************************************/
double wcs_from_pix (hpar_ptr, axis, pixel)
  header_param_list *hpar_ptr;
  int axis;
  double pixel;
{
  double coord;

  coord = (*hpar_ptr).crval[axis-1] 
	  + (pixel-(*hpar_ptr).crpix[axis-1]) * (*hpar_ptr).cdelt[axis-1];
  return coord;
}  

/*****************************************************************************
 **									    **
 **		P I X _ F R O M _ W C S                      		    **
 **									    **
 **	Given a header parameter structure, axis index, and world	    **
 **	coordinate, return the pixel coordinate value for that position.    **
 **	The axis and pixel coord both start counting at 1, not 0!	    **
 **									    **
 **	NOTE: this assumes no projection transformation is required!	    **
 **									    **
 *****************************************************************************/
double pix_from_wcs (hpar_ptr, axis, coord)
  header_param_list *hpar_ptr;
  int axis;
  double coord;
{
  double pixel;

  pixel = (*hpar_ptr).crpix[axis-1] 
	  + (coord-(*hpar_ptr).crval[axis-1]) / (*hpar_ptr).cdelt[axis-1];

  return pixel;
}  

/*****************************************************************************
 **									    **
 **		S E T _ B Z E R O _ B S C A L E              		    **
 **									    **
 **	Set FITS header BSCALE and BZERO parameters to values appropriate   **
 **	for a 16-bit representation covering the min-max range requested.   **
 **									    **
 *****************************************************************************/
void set_bzero_bscale (hpar_ptr, min, max)
  double min, max;
  header_param_list *hpar_ptr;
{
  (*hpar_ptr).bscale = (max - min) / 65534.0;
  (*hpar_ptr).bzero  = (max + min) / 2.0;
}

/*****************************************************************************
 **									    **
 **		S E T _ B Z E R O _ B S C A L E _ D A T A R A N G E	    **
 **									    **
 **	Shell routine to measure the range of data values explicitly and    **
 **	feed these to set_bzero_bscale.					    **
 **									    **
 *****************************************************************************/
void set_bzero_bscale_datarange (hpar_ptr, data)
  float *data;
  header_param_list *hpar_ptr;
{
  int i, number_of_pixels, MINMAX_INIT;
  float min, max;

  number_of_pixels = 1;
  for (i=0; i<(*hpar_ptr).num_axes; i++) {
    number_of_pixels *= (*hpar_ptr).naxis[i];
  }
  MINMAX_INIT = FALSE;
  for (i=0; i<number_of_pixels; i++) {
    if (data[i] < BLANK_THRESH) {
      if (MINMAX_INIT == FALSE) {
	min = data[i];
	max = data[i];
	MINMAX_INIT = TRUE;
      } else {
	if (min > data[i]) min = data[i];
	if (max < data[i]) max = data[i];
      }
    }
  }
  set_bzero_bscale (hpar_ptr, min, max);
}

/*****************************************************************************
 **									    **
 **		I 2 _ T O _ R 4                                    	    **
 **									    **
 **	Return float equivalent of scaled short integer value for specified **
 **	image bscale & bzero header parameter values.			    **
 **									    **
 *****************************************************************************/
static float i2_to_r4 (hpar_ptr, short_val)
  short short_val;
  header_param_list *hpar_ptr;
{
  float float_val;

  if ((*hpar_ptr).blank != 0 && short_val == (*hpar_ptr).blank) {
    float_val = BLANK_PIXEL;
  } else {
    float_val = (*hpar_ptr).bscale * short_val + (*hpar_ptr).bzero;
  }
  return float_val;
}

/*****************************************************************************
 **									    **
 **		R 4 _ T O _ I 2                                    	    **
 **									    **
 **	Return scaled short integer equivalent of float value for specified **
 **	image bscale & bzero header parameter values.			    **
 **									    **
 **	Ensure that I2 scaling range is not exceeded.			    **
 **									    **
 *****************************************************************************/
static short r4_to_i2 (hpar_ptr, float_val, float_min, float_max)
  float float_val, float_min, float_max;
  header_param_list *hpar_ptr;
{
  short short_val;
  float float_scaled;

  if ((*hpar_ptr).blank != 0 && float_val >= BLANK_THRESH) {
    short_val = (*hpar_ptr).blank;
  } else {
    if (float_val < float_min) float_val = float_min;
    if (float_val > float_max) float_val = float_max;
    float_scaled = (float_val - (*hpar_ptr).bzero) / (*hpar_ptr).bscale;
    if (float_scaled >= 0.0) {
      short_val = (short) floor (float_scaled);
    } else {
      short_val = (short) ceil (float_scaled);
    }
  }
  return short_val;
}

/*****************************************************************************
 **									    **
 **		F I N D _ R 4 _ M I N _ M A X                      	    **
 **									    **
 **	Find allowed range of r4 values which map onto a scaled i2 set      **
 **	defined by a supplied header parameter set.               	    **
 **									    **
 *****************************************************************************/
static void find_r4_min_max (hpar_ptr, float_min_ptr, float_max_ptr)
  float *float_min_ptr, *float_max_ptr;
  header_param_list *hpar_ptr;
{
  if ((*hpar_ptr).blank != 0) {
    *float_min_ptr = (float) ((*hpar_ptr).bscale * ((*hpar_ptr).blank + 1) 
							+ (*hpar_ptr).bzero);
  } else {
    *float_min_ptr = (float) ((*hpar_ptr).bscale * (-32768) 
							+ (*hpar_ptr).bzero);
  }
  *float_max_ptr = (float) ((*hpar_ptr).bscale * 32767 + (*hpar_ptr).bzero);
}

/*****************************************************************************
 **									    **
 **		L I M I T _ R 4 _ V A L                            	    **
 **									    **
 **	Force a given float value to lie within the inclusive range         **
 **	defined by two other float values.                        	    **
 **									    **
 *****************************************************************************/
static void limit_r4_val (float_val_ptr, float_min, float_max)
  float *float_val_ptr, float_min, float_max;
{
  if (*float_val_ptr < BLANK_THRESH) {
    if (*float_val_ptr < float_min) *float_val_ptr = float_min;
    if (*float_val_ptr > float_max) *float_val_ptr = float_max;
  }
}

/*****************************************************************************
 **									    **
 **		G E T _ F I T S _ M I N _ M A X                    	    **
 **									    **
 **	Open FITS file, read header & all data values in stream mode to     **
 **	obtain min & max values, close file & return values found.	    **
 **									    **
 *****************************************************************************/
static void get_fits_min_max (fname, min_ptr, max_ptr)
  char fname[];
  double *min_ptr, *max_ptr;
{
  int i, j, k, n1, n2, n3, MINMAX_INIT;
  float min, max, *inplane, data_val;
  header_param_list hpar;
  FILE *infile;

  open_read (&infile, fname);
  readfits_header (infile, &hpar);
  n1 = hpar.naxis[0];
  n2 = hpar.naxis[1];
  n3 = hpar.naxis[2];
  inplane = (float *) malloc (n1 * n2 * sizeof (float));
  MINMAX_INIT = FALSE;
  for (k=0; k<n3; k++) {
    printf ("\tReading plane %d of %d.\n%c[1A", k+1, n3, 27);
    readfits_plane (infile, inplane, &hpar);
    for (j=0; j<n2; j++) {
      for (i=0; i<n1; i++) {
	data_val = inplane[i+j*n1];
	if (data_val < BLANK_THRESH) {
	  if (MINMAX_INIT == FALSE) {
	    min = data_val;
	    max = data_val;
	    MINMAX_INIT = TRUE;
	  } else {
	    if (min > data_val) min = data_val;
	    if (max < data_val) max = data_val;
	  }
	}
      }
    }
  }
  free (inplane);
  fclose (infile);
  *min_ptr = (double) min;
  *max_ptr = (double) max;
}

/*****************************************************************************
 **									    **
 **		G E T _ F I T S _ M I N _ M A X _ G R S             	    **
 **									    **
 **	Open FITS file, read header & all data values in stream mode to     **
 **	obtain min & max values, close file & return values found.	    **
 **									    **
 **	Exclude pixels with values outside a reasonable range as BLANKs.    **
 **	This has been added to accomodate Galactic Ring Survey cubes.	    **
 **									    **
 *****************************************************************************/
static void get_fits_min_max_grs (fname, min_ptr, max_ptr)
  char fname[];
  double *min_ptr, *max_ptr;
{
  int i, j, k, n1, n2, n3, MINMAX_INIT;
  float min, max, *inplane, data_val;
  header_param_list hpar;
  FILE *infile;

  open_read (&infile, fname);
  readfits_header (infile, &hpar);
  n1 = hpar.naxis[0];
  n2 = hpar.naxis[1];
  n3 = hpar.naxis[2];
  inplane = (float *) malloc (n1 * n2 * sizeof (float));
  MINMAX_INIT = FALSE;
  for (k=0; k<n3; k++) {
    printf ("\tReading plane %d of %d.\n%c[1A", k+1, n3, 27);
    readfits_plane (infile, inplane, &hpar);
    for (j=0; j<n2; j++) {
      for (i=0; i<n1; i++) {
	data_val = inplane[i+j*n1];
	if (data_val < BLANK_THRESH) {
	  if (data_val < -100.0 || data_val > 100.0) continue;
	  if (MINMAX_INIT == FALSE) {
	    min = data_val;
	    max = data_val;
	    MINMAX_INIT = TRUE;
	  } else {
	    if (min > data_val) min = data_val;
	    if (max < data_val) max = data_val;
	  }
	}
      }
    }
  }
  free (inplane);
  fclose (infile);
  *min_ptr = (double) min;
  *max_ptr = (double) max;
}

/*****************************************************************************
 **									    **
 **		     G E T _ F I T S _ S T R _ P A R           		    **
 **									    **
 **	Given a FITS 80-character header line, extract a string parameter   **
 **	enclosed in single right quotes, e.g., 'GLAT-CAR'.  Make the best   **
 **	of it if only the first quote is there.  Do not include the quotes  **
 **	in the final string.  Also remove leading and trailing spaces.	    **
 **									    **
 *****************************************************************************/
void get_fits_str_par (line, str)
  char line[], str[];
{
  int i, j, k, len;

  /* find 1st instance of (') character; if none, return empty string */
  str[0] = 0;
  i = 0;
  while (i < LINELEN && line[i] != '\'') i++;
  if (i < LINELEN) {
    /* find 2nd instance of (') character */
    j = i+1;
    while (j < LINELEN && line[j] != '\'') j++;
    if (j >= LINELEN) {
      /* if 2nd ' missing, look for (/) comment character */
      j = i+1;
      while (j < LINELEN && line[j] != '/') j++;
      if (j >= LINELEN) {
	/* if no comment character either, use last non-space character */
	j = LINELEN-1;
	while (j > i && line[j] == ' ') j--;
	j++;
      }
    }
    for (k=i+1; k<j; k++) str[k-i-1] = line[k];
    str[j-i-1] = 0;
    len = j-i-1;

    /* remove leading spaces */
    k = 0;
    while (k < len && str[k] == ' ') k++;
    for (i=k; i<len; i++) str[i-k] = str[i];
    str[len-k] = 0;
    len = len-k;

    /* remove trailing spaces */
    j = len - 1;
    while (j > 0 && str[j] == ' ') j--;
    str[j+1] = 0;
  }
}

/*****************************************************************************
 **									    **
 **		     R E A D F I T S _ H E A D E R			    **
 **									    **
 **	Similar to read_fits_header except the file is already open.	    **
 **									    **
 *****************************************************************************/
void readfits_header (file_ptr, hpar_ptr)
  FILE *file_ptr;
  header_param_list *hpar_ptr;
{
  char buf[LINELEN+1], flagstr[LINELEN+1];
  int i, j, axis, bytecount, linecount, blockcount, END_FOUND, 
	datarawmin, datarawmax;

  init_header_param_list (hpar_ptr);	/* initializations */
  (*hpar_ptr).ncomment = 0;
  (*hpar_ptr).nhistory = 0;
  bytecount = 0;
  linecount = 0;
  blockcount = 0;
  END_FOUND = FALSE;

  do {
    for (i=0; i<LINESPERBLOCK; i++) {
      fread (buf, sizeof(char), LINELEN, file_ptr); 
      buf[LINELEN] = 0;  /* end of string character */
      if (strncmp (buf, "SIMPLE  ", 8)==0) {
	sscanf(&buf[10], "%s", flagstr);
	if (flagstr[0]=='T') (*hpar_ptr).simple = TRUE;
	if (flagstr[0]=='F') (*hpar_ptr).simple = FALSE;
      }
      if (strncmp (buf, "BLOCKED ", 8)==0) {
	sscanf(&buf[10], "%s", flagstr);
	if (flagstr[0]=='T') (*hpar_ptr).blocked = TRUE;
	if (flagstr[0]=='F') (*hpar_ptr).blocked = FALSE;
      }
      if (strncmp (buf, "EXTEND  ", 8)==0) {
	sscanf(&buf[10], "%s", flagstr);
	if (flagstr[0]=='T') (*hpar_ptr).extend = TRUE;
	if (flagstr[0]=='F') (*hpar_ptr).extend = FALSE;
      }
      if (strncmp(buf,"BITPIX  ",8)==0) 
			sscanf(&buf[10], "%d" , &(*hpar_ptr).bitpix);
      if (strncmp(buf,"BUNIT   ",8)==0) 
			get_fits_str_par (buf, (*hpar_ptr).bunit);
      if (strncmp(buf,"BSCALE  ",8)==0) 
			sscanf(&buf[10], "%le", &(*hpar_ptr).bscale);
      if (strncmp(buf,"BZERO   ",8)==0) 
			sscanf(&buf[10], "%le", &(*hpar_ptr).bzero);
      if (strncmp(buf,"BLANK   ",8)==0) 
			sscanf(&buf[10], "%d" , &(*hpar_ptr).blank);
      if (strncmp(buf,"DATAMAX ",8)==0) 
			sscanf(&buf[10], "%le", &(*hpar_ptr).datamax);
      if (strncmp(buf,"DATAMIN ",8)==0) 
			sscanf(&buf[10], "%le", &(*hpar_ptr).datamin);
      if (strncmp(buf,"OBSFREQ ",8)==0) 
			sscanf(&buf[10], "%le", &(*hpar_ptr).obsfreq);
      if (strncmp(buf,"EXPTIME ",8)==0) 
			sscanf(&buf[10], "%le", &(*hpar_ptr).exptime);
      if (strncmp(buf,"EPOCH   ",8)==0) 
			sscanf(&buf[10], "%le", &(*hpar_ptr).epoch);
      if (strncmp(buf,"EQUINOX ",8)==0) 
			sscanf(&buf[10], "%le", &(*hpar_ptr).equinox);
      if (strncmp(buf,"DATE"    ,4)==0) 
			get_fits_str_par (buf, (*hpar_ptr).date);
      if (strncmp(buf,"OBJECT  ",8)==0) 
			get_fits_str_par (buf, (*hpar_ptr).object);
      if (strncmp(buf,"ORIGIN  ",8)==0) 
			get_fits_str_par (buf, (*hpar_ptr).origin);
      if (strncmp(buf,"TELESCOP",8)==0) 
			get_fits_str_par (buf, (*hpar_ptr).telescope);
      if (strncmp(buf,"INSTRUME",8)==0) 
			get_fits_str_par (buf,(*hpar_ptr).instrument);
      if (strncmp(buf,"OBSERVER",8)==0) 
			get_fits_str_par (buf, (*hpar_ptr).observer);
      if (strncmp(buf,"FILENAME",8)==0) 
			get_fits_str_par (buf, (*hpar_ptr).filename);
      if (strncmp(buf,"NAXIS   ",8)==0) 
			sscanf(&buf[10],"%d",&(*hpar_ptr).num_axes);
      if (strncmp(buf,"NAXIS"   ,5)==0 && buf[5] >= '1' && buf[5] <= '7') {
	sscanf (&buf[5], "%d", &axis);
	if (axis <= (*hpar_ptr).num_axes) 
			sscanf (&buf[10], "%d", &(*hpar_ptr).naxis[axis-1]);
      }
      if (strncmp(buf,"CTYPE"   ,5)==0 && buf[5] >= '1' && buf[5] <= '7') {
	sscanf (&buf[5], "%d", &axis);
	if (axis <= (*hpar_ptr).num_axes) 
			get_fits_str_par (buf, (*hpar_ptr).ctype[axis-1]);
      }
      if (strncmp(buf,"CRVAL"   ,5)==0 && buf[5] >= '1' && buf[5] <= '7') {
	sscanf (&buf[5], "%d", &axis);
	if (axis <= (*hpar_ptr).num_axes) 
			sscanf (&buf[10],"%le",&(*hpar_ptr).crval[axis-1]);
      }
      if (strncmp(buf,"CRPIX"   ,5)==0 && buf[5] >= '1' && buf[5] <= '7') {
	sscanf (&buf[5], "%d", &axis);
	if (axis <= (*hpar_ptr).num_axes) 
			sscanf (&buf[10],"%le",&(*hpar_ptr).crpix[axis-1]);
      }
      if (strncmp(buf,"CDELT"   ,5)==0 && buf[5] >= '1' && buf[5] <= '7') {
	sscanf (&buf[5], "%d", &axis);
	if (axis <= (*hpar_ptr).num_axes) 
			sscanf (&buf[10],"%le",&(*hpar_ptr).cdelt[axis-1]);
      }
      if (strncmp(buf,"CROTA"   ,5)==0 && buf[5] >= '1' && buf[5] <= '7') {
	sscanf (&buf[5], "%d", &axis);
	if (axis <= (*hpar_ptr).num_axes) 
			sscanf (&buf[10],"%le",&(*hpar_ptr).crota[axis-1]);
      }
      if (strncmp(buf,"COMMENT" ,7)==0) {
	if ((*hpar_ptr).ncomment >= MAX_COMMENT) {
	  fprintf (stderr, "*** error in READFITS_HEADER: MAX COMMENT COUNT");
	  fprintf (stderr, " OF %d EXCEEDED! ***\n", MAX_COMMENT);
	  exit (1);
	}
	for (j=8; j<LINELEN; j++) 
			(*hpar_ptr).comment[(*hpar_ptr).ncomment][j-8] =buf[j];
	(*hpar_ptr).ncomment++;
      }
      if (strncmp(buf,"HISTORY" ,7)==0) {
	if ((*hpar_ptr).nhistory >= MAX_HISTORY) {
	  fprintf (stderr, "*** error in READFITS_HEADER: MAX HISTORY COUNT");
	  fprintf (stderr, " OF %d EXCEEDED! ***\n", MAX_HISTORY);
	  exit (1);
	}
	for (j=8; j<LINELEN; j++) 
			(*hpar_ptr).history[(*hpar_ptr).nhistory][j-8] =buf[j];
	(*hpar_ptr).nhistory++;
      }
      if (strncmp(buf,"END",3)==0) END_FOUND = TRUE;
    }
    blockcount++;
  } while (END_FOUND == FALSE);

  if (VERBOSE==TRUE) {
    /* output file parameters & header parameters of interest */
    printf ("total %d bytes, %d lines and %d blocks in FITS header\n", 
		bytecount, linecount, blockcount);
    printf ("bitpix = %d", (*hpar_ptr).bitpix);
    for (axis=0; axis<(*hpar_ptr).num_axes; axis++) {
      printf ("  naxis%d = %d", axis, (*hpar_ptr).naxis[axis]);
    }
    printf ("\n");
    printf ("  bscale = %e", (*hpar_ptr).bscale);
    printf ("  bzero = %e",  (*hpar_ptr).bzero);
    printf ("  bunit = %s",  (*hpar_ptr).bunit);
    printf ("  blank = %d",  (*hpar_ptr).blank);
    printf ("\n");
    for (axis=0; axis<(*hpar_ptr).num_axes; axis++) {
      printf ("  crval%d = %e", axis, (*hpar_ptr).crval[axis]);
      printf ("  crpix%d = %e", axis, (*hpar_ptr).crpix[axis]);
      printf ("  cdelt%d = %e", axis, (*hpar_ptr).cdelt[axis]);
      printf ("  crota%d = %e", axis, (*hpar_ptr).crota[axis]);
      printf ("\n");
    }
    printf ("FITS given datamax = %e  datamin = %e\n", 
				(*hpar_ptr).datamax, (*hpar_ptr).datamin);
    datarawmin = floor ( ((*hpar_ptr).datamin 
				- (*hpar_ptr).bzero) / (*hpar_ptr).bscale );
    datarawmax = ceil  ( ((*hpar_ptr).datamax 
				- (*hpar_ptr).bzero) / (*hpar_ptr).bscale );
    printf ("predicted rawmax = %d = %X  rawmin = %d = %X\n",
		datarawmax, datarawmax, datarawmin, datarawmin);
  }
}

/*****************************************************************************
 **									    **
 **		    W R I T E F I T S _ H E A D E R         		    **
 **									    **
 *****************************************************************************/
void writefits_header (file_ptr, hpar_ptr)
  FILE *file_ptr;
  header_param_list *hpar_ptr;
{

  int i, j, k, hblocks, hlines, used_lines;

  if ((*hpar_ptr).num_axes > 4) {
    fprintf (stderr, 
	"*** error in WRITEFITS_HEADER: %d axes not supported! ***\n", 
	(*hpar_ptr).num_axes);
    exit(1);
  }
  if ((*hpar_ptr).bitpix == 16) {
    if ((*hpar_ptr).blank != -32768) {
      if (WARN_BLANK_CHANGE == TRUE) {
	fprintf (stderr, "  WRITEFITS warning: forcing BLANK = -32768\n");
      }
      (*hpar_ptr).blank = -32768;
    }
  }
  file_ptr = file_ptr;    
  hblocks = 1;
  hlines = hblocks * LINESPERBLOCK;

		/*   1234567890123456789012345678901234567890 */
  fprintf (file_ptr, "SIMPLE  =                    T / ");
  fprintf (file_ptr, " Standard FITS file                            ");
  fprintf (file_ptr, "BITPIX  =                  %3d / ", (*hpar_ptr).bitpix);
  switch ((*hpar_ptr).bitpix) {
    case 16 : 
      fprintf (file_ptr, " IEEE 2-byte scaled integer                    ");
      break;
    case -32 : 
      fprintf (file_ptr, " IEEE 4-byte float                             ");
      break;
    default :
      fprintf (file_ptr, " *** Unknown pixel type!  Hang on!! ***        ");
      fprintf (stderr, "%%%% Warning: bad BITPIX in WRITEFITS_HEADER %%%%\n");
  }
  fprintf (file_ptr, "NAXIS   =                    %1d / ", 
							(*hpar_ptr).num_axes);
  fprintf (file_ptr, " Number of image dimensions                    ");
  fprintf (file_ptr, "NAXIS1  =              %7d / ", (*hpar_ptr).naxis[0]);
  fprintf (file_ptr, " Size of 1st dimension in pixels               ");
  fprintf (file_ptr, "NAXIS2  =              %7d / ", (*hpar_ptr).naxis[1]);
  fprintf (file_ptr, " Size of 2nd dimension in pixels               ");
  if ((*hpar_ptr).num_axes >= 3) {
    fprintf (file_ptr, "NAXIS3  =              %7d / ", (*hpar_ptr).naxis[2]);
    fprintf (file_ptr, " Size of 3rd dimension in pixels               ");
  }
  if ((*hpar_ptr).num_axes >= 4) {
    fprintf (file_ptr, "NAXIS4  =              %7d / ", (*hpar_ptr).naxis[3]);
    fprintf (file_ptr, " Size of 4th dimension in pixels               ");
  }
  /* find length of OBJECT string; truncate if too long */
  i = 0;   while (i < 53 && (*hpar_ptr).object[i] != 0) i++;
  (*hpar_ptr).object[i] = 0;
  fprintf (file_ptr, "OBJECT  = '%s'", (*hpar_ptr).object);
  /* pad remainder of OBJECT string area with spaces */
  for (j=0; j<53-i; j++) fprintf (file_ptr, " ");
  fprintf (file_ptr, " / ");
  fprintf (file_ptr, " Object name");

  for (k=0; k<(*hpar_ptr).num_axes; k++) {
    /* ensure that CTYPE strings are also the proper length */
    i = 0;  while (i < 18 && (*hpar_ptr).ctype[k][i] != 0) i++;
    (*hpar_ptr).ctype[k][i] = 0;
    fprintf (file_ptr, "CTYPE%1d  = '%s'", k+1, (*hpar_ptr).ctype[k]);
    for (j=0; j<18-i; j++) fprintf (file_ptr, " "); fprintf (file_ptr, " / ");
    if (k==0) 
      fprintf (file_ptr, " 1st axis type                                 ");
    if (k==1) 
      fprintf (file_ptr, " 2nd axis type                                 ");
    if (k==2) 
      fprintf (file_ptr, " 3rd axis type                                 ");
    if (k==3) 
      fprintf (file_ptr, " 4th axis type                                 ");
    fprintf (file_ptr, "CRVAL%1d  = %20.6f / ", k+1, (*hpar_ptr).crval[k]);
    fprintf (file_ptr, " Reference pixel value                         ");
    fprintf (file_ptr, "CRPIX%1d  = %20.2f / ", k+1, (*hpar_ptr).crpix[k]);
    fprintf (file_ptr, " Reference pixel                               ");
    fprintf (file_ptr, "CDELT%1d  = %20.7f / ", k+1, (*hpar_ptr).cdelt[k]);
    fprintf (file_ptr, " Pixel size in world coordinate units          ");
    fprintf (file_ptr, "CROTA%1d  = %20.4f / ", k+1, (*hpar_ptr).crota[k]);
    fprintf (file_ptr, " Axis rotation in degrees                      ");
  }
  fprintf (file_ptr, "EQUINOX = %20.2f / ", (*hpar_ptr).equinox);
  fprintf (file_ptr, " Equinox of coordinates (if any)               ");
  /* find length of BUNIT string; truncate if too long */
  i = 0;   while (i < 38 && (*hpar_ptr).bunit[i] != 0) i++;
  (*hpar_ptr).bunit[i] = 0;
  fprintf (file_ptr, "BUNIT   = '%s'", (*hpar_ptr).bunit);
  /* pad remainder of BUNIT string area with spaces */
  for (j=0; j<38-i; j++) fprintf (file_ptr, " ");
  fprintf (file_ptr, " / ");
  fprintf (file_ptr, " Units of pixel data values");
  if ((*hpar_ptr).bitpix == 16) {
    fprintf (file_ptr, "BSCALE  = %20.7E / ", (*hpar_ptr).bscale);
    fprintf (file_ptr, " Real pixel value = RAW * BSCALE + BZERO       ");
    fprintf (file_ptr, "BZERO   = %20.7E / ", (*hpar_ptr).bzero);
    fprintf (file_ptr, "                                               ");
    fprintf (file_ptr, "BLANK   = %20d / ", (*hpar_ptr).blank);
    fprintf (file_ptr, " Raw pixel value indicating no data            ");
  }
  for (i=0; i<(*hpar_ptr).ncomment; i++) {
    fprintf (file_ptr, "COMMENT %s", (*hpar_ptr).comment[i]);
    for (j=8+strlen((*hpar_ptr).comment[i]); j<LINELEN; j++) {
      fprintf (file_ptr, " ");
    }
  }
  for (i=0; i<(*hpar_ptr).nhistory; i++) {
    fprintf (file_ptr, "HISTORY %s", (*hpar_ptr).history[i]);
    for (j=8+strlen((*hpar_ptr).history[i]); j<LINELEN; j++) {
      fprintf (file_ptr, " ");
    }
  }
  fprintf (file_ptr, "END                              ");
  fprintf (file_ptr, "                                               ");

  /* pad out rest of header block */
  used_lines = 7+6*(*hpar_ptr).num_axes 
		+ (*hpar_ptr).ncomment + (*hpar_ptr).nhistory;
  if ((*hpar_ptr).bitpix == 16) used_lines += 3;
/*
  for (i=used_lines; i<hlines; i++) {
*/
  for (i = 0; i < hlines - (used_lines % hlines); i++) {
    fprintf (file_ptr, "                                 ");
    fprintf (file_ptr, "                                               ");
  }
}

/*****************************************************************************
 **									    **
 **			 R E A D F I T S _ P L A N E		            **
 **									    **
 **	Read a single "plane" of FITS data, specified by a number of 	    **
 **	pixels and the bits per pixel.  The file is presumed to already     **
 **	be open and pointing to the first pixel of the plane.		    **
 **									    **
 **	In addition, the 2-D array to be written to is presumed to already  **
 **	be allocated and ready for use.					    **
 **									    **
 *****************************************************************************/
void readfits_plane (file_ptr, data, hpar_ptr)
  float *data;
  FILE *file_ptr;
  header_param_list *hpar_ptr;
{

  char *pixbytes, byteblock[BLOCKSIZE];
  short shortval;
  int i, j, k, BLANK_USED, 
	blank, bytes_per_pixel, datablocks, pixels_per_block, pixels_in_block, 
	bytes_in_block, pixels_in_last, bytes_in_last, empties, npix, bitpix, 
	rawval, rawmin, rawmax, expon, ffcount;
  float datum_f;
  double bscale, bzero, compmax, compmin, compsum, compmean;

  bitpix = (*hpar_ptr).bitpix;
  npix = (*hpar_ptr).naxis[0] * (*hpar_ptr).naxis[1];
  bzero = (*hpar_ptr).bzero;
  bscale = (*hpar_ptr).bscale;
  blank = (*hpar_ptr).blank;
  bytes_per_pixel = floor (1.0 * abs(bitpix) / 8.0);
  datablocks = (int) ceil (1.0 * npix * bytes_per_pixel / (1.0 * BLOCKSIZE));
  pixels_per_block = BLOCKSIZE / bytes_per_pixel;
  pixels_in_last = npix - pixels_per_block * (datablocks - 1);
  bytes_in_last = pixels_in_last * bytes_per_pixel;
  if (blank == 0) {
    BLANK_USED = FALSE;
  } else {
    BLANK_USED = TRUE;
  }
  if (VERBOSE==TRUE) {
    printf ("BITPIX = %d\n", bitpix);
    printf ("BZERO  = %e\n", bzero);
    printf ("BSCALE = %e\n", bscale);
    printf ("BLANK  = %d\n", blank);
  }

  empties = 0;
  compmax = -1.e15;
  compmin = 1.e15;
  compsum = 0.0;
  compmean = 0.0;

  /* Integer I/O:
       FITS stores bytes consecutively in order of significance, from most to
       least.  Machines, in memory, do not.  This loop below reads in each set
       of bytes from the FITS file and ARITHMETICALLY places them one at a time
       into the waiting integer.  This way the compiler actually decides where
       in memory the proper bytes go, and the algorithm is machine-independent.
       The & 0xff operation ensures that when each byte is converted to a 
       temporary integer, there are no leading sign bits -- all bytes in the
       temp int are 00 except the LS byte.  This is NOT done for the MS byte
       when it is read in, because the sign bits are needed if the bytes do 
       not fill the waiting integer.  WARNING: this mapping is only valid 
       for bytes_per_pixel <= 4 = sizeof(int). */

  /* Real I/O:
       IRAF will write 4-byte reals (floats) if BITPIX = -32 is specified in
       wfits, and 8-byte reals (doubles) if BITPIX = -64.  Only the 4-byte
       reals are implemented here.  The IEEE format is

       	     sign  exponent          mantissa  
	      bit   8 bits     23 bits (24 implicit)

         	S EEE,E EEE,E MMM,MMMM,MMMM,MMMM,MMMM,MMMM

	where the exponent (range 00-FF) has a built-in positive bias of 127,
	so the true exponent is found by subtracting 7F.  The sign bit is 0
	for positive mantissa, 1 for negative.  The mantissa's most significant
	bit is always 1, and therefore not stored.  Two examples are

	            1.0 = 0 011,1 111,1 000,0000,0000,0000,0000,0000 = 3F800000
	  	    (in base 2) exp = 0 and mant = 1.0000...

	-3/16 = -0.1875 = 1 011,1 110,0 100,0000,0000,0000,0000,0000 = BE400000
	 	   (in base 2) exp = -3 and mant = 1.1000...

	The same technique described for the integers is used for the floats to
	let the compiler decide how to order bytes on the machine.  Two items
	then remain.  On the VAX, Stephan says that the exponent bias is off 
	by 2 from IEEE standard, but the DECstations, Alphas, Suns & Indigos
	all have it as it should be.  Finally, the rawval integer bytes must
	be mapped onto a float (NOT a double!) by recasting with pointers; 
	simple recasting won't do, and arithmetic assignment definitely won't.
	Then a double is arithmetically assigned the value of the float. */

  ffcount = 0;
  /*linear (x,y) order: (0,0),(1,0),...(naxis1,0),(0,1),...(naxis1,naxis2) */
  for (j=0; j<datablocks; j++) {
    if (j < datablocks-1) {
      bytes_in_block = BLOCKSIZE;
      pixels_in_block = pixels_per_block;
    } else {
      /* do not load bytes beyond number of pixels requested */
      bytes_in_block = bytes_in_last;
      pixels_in_block = pixels_in_last;
    }
    fread (byteblock, sizeof(char), bytes_in_block, file_ptr);
    for (i=0; i<pixels_in_block; i++) {
      pixbytes = &byteblock[i*bytes_per_pixel];
      if (VERBOSE==TRUE) {     /* output 1st 5 raw pixel values in hex */
        if (ffcount < 5) {
	  printhex (pixbytes, abs(bytes_per_pixel), 0); printf(" ");
	  ffcount++;
        }
      }
      if (bitpix==16) {
        /* KLUDGE -- force proper treatment of 16-bit signed integers */
        shortval = (short)pixbytes[0];         /* MS byte in LS pos */
        for (k=1; k<bytes_per_pixel; k++) {    /* add in any extra bytes*/
	  shortval <<= 8;	       	       /* shift everybody left 1 byte*/
	  shortval |= (((short) pixbytes[k]) & 0xff); /* next byte in LS pos */
        }
        rawval = (int) shortval;
      } else {
        rawval = (int)pixbytes[0];             /* MS byte in LS pos */
        for (k=1; k<bytes_per_pixel; k++) {    /* add in any extra bytes*/
	  rawval <<= 8;	       	               /* shift everybody left 1 byte*/
 	  rawval |= (((int) pixbytes[k]) & 0xff); /* next byte in LS pos */
        }
      }
      if (bitpix==8) {
	rawval = (int) pixbytes[0];
	if (rawval < 0) rawval += 256;  /* for unsigned values >= 0x80 */
      }
/*
printhex (pixbytes, 1, 0);
printf(" %4d %4d %9.4f\n", (int) pixbytes[0], rawval, bscale * rawval + bzero);
*/
      if (VERBOSE==TRUE && ffcount < 5) printf (" %6d ", rawval);
      if (i==0 && j==0) rawmax = rawmin = rawval;
      switch (bitpix) {
        case 8: case 16: case 32:     /* unsigned, short, & long int */
    	  if (BLANK_USED==TRUE && rawval==blank) {
	    datum_f = BLANK_PIXEL;
	    empties++;
	  } else {
	    if (rawval > rawmax) rawmax = rawval;
	    if (rawval < rawmin) rawmin = rawval;
	    datum_f = bscale * rawval + bzero;
	  }
          if (VERBOSE==TRUE && ffcount < 5) printf (" %e", datum_f);
	  break;
        case -32:	/* single-precision real (IEEE-format float) */
	  if (VMS==TRUE) {
	    expon = (rawval & 0x7f800000) >> 23; /* extract exponent */
	    expon += 2;  /* remove bias, reinsert corrected exponent */
	    rawval = (rawval & 0x807fffff) | (expon << 23);
	  }
	  datum_f = * ((float *) &rawval); /* map onto float var */
	  if (isnan(datum_f)) datum_f = BLANK_PIXEL;
	  break;
        default: 
	  printf ("unexpected BITPIX value %d in case statement\n", bitpix);
	  break;
      } /* end switch */
      if (datum_f > compmax) compmax = datum_f;
      if (datum_f < compmin) compmin = datum_f;
      compsum += datum_f;
      data[i+j*pixels_per_block] = datum_f;
      if (VERBOSE==TRUE && ffcount < 5) printf ("  data[%4d] = %e\n", 
		i+j*pixels_per_block, data[i+j*pixels_per_block]);
    }
  }
  compmean = compsum / (npix - empties);

  if (VERBOSE==TRUE) {
    /* output computed data parameters */
    if (bitpix > 0) {
      printf ("%d empty pixels\n", empties);
      printf ("raw max = %d = %X  raw min = %d = %X\n",
						rawmax,rawmax,rawmin,rawmin);
    }
    printf("computed max = %e  min = %e  mean = %e\n",
						compmax,compmin,compmean);
  }
}

/*****************************************************************************
 **									    **
 **			W R I T E F I T S _ P L A N E			    **
 **									    **
 **	Write a single "plane" of FITS data, specified by a number of 	    **
 **	pixels and the bits per pixel.  The file is presumed to already     **
 **	be open and pointing to the location of the first pixel to be 	    **
 **	written.  This routine expects a **FLOAT** array of data values.    **
 **									    **
 *****************************************************************************/
void writefits_plane (file_ptr, data, hpar_ptr)
  float *data;
  FILE *file_ptr;
  header_param_list *hpar_ptr;
{

  char 	pixbytes[4],	/* 4-byte hex string for fits output */
    	byteblock[BLOCKSIZE];
  short shortval;    /* shortval is a short containing the machine-format hex*/
  int 	rawval,	     /* rawval is an int containing the machine-format hex */
	i, j, k, expon, bitpix, npix, blank, 
	bytes_per_pixel, pixels_per_block, pixels_in_last, bytes_in_last,
	bytes_in_block, pixels_in_block, datablocks;
  float datum_f,	/* float to map onto rawval */
	bzero, bscale, float_scaled;

  bitpix = (*hpar_ptr).bitpix;
  npix = (*hpar_ptr).naxis[0] * (*hpar_ptr).naxis[1];
  bzero = (*hpar_ptr).bzero;
  bscale = (*hpar_ptr).bscale;
  blank = (*hpar_ptr).blank;
  bytes_per_pixel = floor (1.0 * abs (bitpix) / 8.0);
  datablocks = (int) ceil (1.0 * npix * bytes_per_pixel / (1.0 * BLOCKSIZE));
  pixels_per_block = BLOCKSIZE / bytes_per_pixel;
  pixels_in_last = npix - pixels_per_block * (datablocks - 1);
  bytes_in_last = pixels_in_last * bytes_per_pixel;

  /* fits conversion kludgy, less general than reader version */
  for (j=0; j<datablocks; j++) {
    if (j < datablocks-1) {
      bytes_in_block = BLOCKSIZE;
      pixels_in_block = pixels_per_block;
    } else {
      /* do not load bytes beyond number of pixels requested */
      bytes_in_block = bytes_in_last;
      pixels_in_block = pixels_in_last;
    }
    for (i=0; i<pixels_in_block; i++) {
      datum_f = data[i+j*pixels_per_block];
      switch (bitpix) {
	case  16 :
	  if (datum_f >= BLANK_THRESH) {
	    shortval = blank;
	  } else {
	    float_scaled = (datum_f - bzero) / bscale;
	    if (float_scaled >= 0.0) {
	      rawval = (int) floor (float_scaled);
	    } else {
	      rawval = (int) ceil (float_scaled);
	    }
	    if (rawval <      blank + 1 ) rawval =       blank + 1 ;
	    if (rawval > abs (blank + 1)) rawval =  abs (blank + 1);
	    shortval = (short) rawval;
	  }
	  pixbytes[0] = (char) ((shortval & 0xff00) >> 8);
	  pixbytes[1] = (char) ((shortval & 0x00ff));
	  break;
	case -32 :
	  /* Save IEEE float values if not blank, NaN, -inf, or +inf.
	     The BLANK threshhold is set by a constant defined elsewhere.
	     On standard IRIX and linux C compilers, the largest "finite"
	     float value appears to be +/- 3.4028234664e+38.  The level 
	     of +/- 1.0e+37 hardwired here backs off from this a little
	     in case some other compilers don't go as far.  -inf and +inf
	     values will fail one of the threshhold tests, and NaN values
	     will fail all of them -- the math library works that way. */
	  if (   datum_f < BLANK_THRESH	
	      && datum_f > -1.0e+37 && datum_f < 1.0e+37) {
	    rawval = * ((int *) &datum_f);
	    if (VMS==TRUE) {  /* correct the exponent bias -- see FITS reader*/
	      expon = (rawval & 0x7f800000) >> 23;
	      expon -= 2;
	      rawval = (rawval & 0x807fffff) | (expon << 23);
	    }
	    pixbytes[0] = (char) ((rawval & 0xff000000) >> 24);
	    pixbytes[1] = (char) ((rawval & 0x00ff0000) >> 16);
	    pixbytes[2] = (char) ((rawval & 0x0000ff00) >>  8);
	    pixbytes[3] = (char)  (rawval & 0x000000ff);
	  } else {
		/* IEEE NaN, commonly interpreted as FITS BLANK for BITPIX=-32,
		   e.g., by Richard Gooch's Karma software. */
	    for (k=0; k<4; k++) pixbytes[k] = 0xFF;
	  }
	  break;
	default :
	  fprintf (stderr, 
		"error in WRITEFITS_PLANE: bitpix = %d unsupported.\n",bitpix);
	  exit (1);
      }
      if ((VERBOSE==TRUE) && (i<5)) printhex (pixbytes, bytes_per_pixel, 1);
      for (k=0; k<bytes_per_pixel; k++) {
	byteblock[k+i*bytes_per_pixel] = pixbytes[k];
      }
    }
    fwrite (byteblock, sizeof(char), bytes_in_block, file_ptr); 
  }
}

/*****************************************************************************
 **									    **
 **			W R I T E F I T S _ P A D    			    **
 **									    **
 **	Write a specified number of NULL bytes to an open FITS file.  	    **
 **	This is to fill up the remainder of the last data block.            **
 **	The calling routine must determine how many bytes to pad out. 	    **
 **									    **
 *****************************************************************************/
void writefits_pad (file_ptr, nbytes)
  int nbytes;
  FILE *file_ptr;
{
  int i;

  for (i=0; i<nbytes; i++) fputc (0, file_ptr);
}

/*****************************************************************************
 **									    **
 **			W R I T E F I T S _ P A D _ E N D   	       	    **
 **									    **
 **	Write the necessary number of NULL bytes to an open FITS file. 	    **
 **	This is to fill up the remainder of the last data block.            **
 **	The routine presumes all the data has been written, and it	    **
 **	calculates the number of bytes to pad based on this assumption.     **
 **									    **
 *****************************************************************************/
void writefits_pad_end (file_ptr, hpar_ptr)
  header_param_list *hpar_ptr;
  FILE *file_ptr;
{
  int i, data_bytes, dblocks, pad_bytes;

  data_bytes = abs ((*hpar_ptr).bitpix / 8);
  for (i=0; i<(*hpar_ptr).num_axes; i++) data_bytes *= (*hpar_ptr).naxis[i];
  dblocks    = (int) ceil (1.0 * data_bytes / BLOCKSIZE);
  pad_bytes  = dblocks * BLOCKSIZE - data_bytes;
  writefits_pad (file_ptr, pad_bytes);
}

/*****************************************************************************
 **									    **
 **		R E A D F I T S _ P L A N E _ A L L O C A T E	            **
 **									    **
 **	Read a single "plane" of FITS data, specified by a number of 	    **
 **	pixels and the bits per pixel.  The file is presumed to already     **
 **	be open and pointing to the first pixel of the plane.		    **
 **									    **
 **	This is a "shell" routine which calls READFITS_PLANE.               **
 **	It allocates a 2-D array to hold the FITS data prior to this call,  **
 **     and returns a pointer to the array when finished.                   **
 **									    **
 *****************************************************************************/
void readfits_plane_allocate (file_ptr, data_ptr, hpar_ptr)
  float **data_ptr;
  FILE *file_ptr;
  header_param_list *hpar_ptr;
{

  (*data_ptr) = (float *) malloc ((*hpar_ptr).naxis[0] * (*hpar_ptr).naxis[1] 
							* sizeof (float));
  readfits_plane (file_ptr, *data_ptr, hpar_ptr);
}

/*****************************************************************************
 **									    **
 **			 R E A D F I T S _ M A P		            **
 **									    **
 **	Read 2-dimensional FITS image (of floats) and its header.           **
 **                                                                         **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void readfits_map (fname, data_ptr, hpar_ptr)
  char *fname;
  float **data_ptr;
  header_param_list *hpar_ptr;
{
  FILE *infile;

  open_read (&infile, fname);
  readfits_header (infile, hpar_ptr);
  if ((*hpar_ptr).num_axes < 2) {
    fprintf (stderr, "*** error in READFITS_MAP: num_axes = %d < 2 ***\n",
							(*hpar_ptr).num_axes);
    exit (1);
  }
  readfits_plane_allocate (infile, data_ptr, hpar_ptr);
  fclose (infile);
}

/*****************************************************************************
 **									    **
 **		      R E A D F I T S _ A V G _ M A P		            **
 **									    **
 **	Read list of 2-D FITS images in sequence, computing a weighted      **
 **     average of the set on the fly, and returning the average map.       **
 **									    **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void readfits_avg_map (num_in, fname_list, weight_list, avgmap_ptr, hpar_ptr)
  char fname_list[MAX_NAVG][MAX_PARSE_LENGTH+1];
  int num_in;
  float **avgmap_ptr, weight_list[];
  header_param_list *hpar_ptr;
{
  char object_old[LINELEN];
  int i, j, k, n1, n2;
  float *indata, *tweight_map;
  FILE *infile;

  for (k=0; k<num_in; k++) {
    printf ("  Averaging map %s with weight %.4f\n", 
						fname_list[k], weight_list[k]);
    open_read (&infile, fname_list[k]);
    readfits_header (infile, hpar_ptr);
    if ((*hpar_ptr).num_axes < 2) {
      fprintf (stderr, 
		"*** error in READFITS_AVG_MAP: num_axes = %d < 2 ***\n",
							(*hpar_ptr).num_axes);
      exit (1);
    }
    if (k==0) {
      n1 = (*hpar_ptr).naxis[0];
      n2 = (*hpar_ptr).naxis[1];
      (*avgmap_ptr) = (float *) malloc (n1 * n2 * sizeof (float));
      tweight_map = (float *) malloc (n1 * n2 * sizeof (float));
      for (i=0; i<n1*n2; i++) {
	(*avgmap_ptr)[i] = BLANK_PIXEL;
	tweight_map[i] = 0.0;
      }
      indata = (float *) malloc (n1 * n2 * sizeof (float));
    } else {
      if (n1 != (*hpar_ptr).naxis[0] || n2 != (*hpar_ptr).naxis[1]) {
	fprintf (stderr, 
	       "*** error in READFITS_AVG_MAP: inconsistent dimensions ***\n");
	exit (1);
      }
    }
    readfits_plane (infile, indata, hpar_ptr);
    fclose (infile);
    for (j=0; j<n2; j++) {
      for (i=0; i<n1; i++) {
	if (indata[i+j*n1] < BLANK_THRESH) {
	  if ((*avgmap_ptr)[i+j*n1] >= BLANK_THRESH) {
	    (*avgmap_ptr)[i+j*n1] = weight_list[k] * indata[i+j*n1];
	  } else {
	    (*avgmap_ptr)[i+j*n1] += weight_list[k] * indata[i+j*n1];
	  }
	  tweight_map[i+j*n1] += weight_list[k];
	}
      }
    }
  }
  free (indata);
  for (i=0; i<n1*n2; i++) {
    if (tweight_map[i] != 0.0) (*avgmap_ptr)[i] /= tweight_map[i];
  }
  free (tweight_map);
  strcpy (object_old, (*hpar_ptr).object);
  sprintf ((*hpar_ptr).object, "Weighted Inverse Average %s", object_old);
}

/*****************************************************************************
 **									    **
 **		      R E A D F I T S _ A V G _ M A P 2		            **
 **									    **
 **	Read list of 2-D FITS images in sequence, computing a weighted      **
 **     average of the set on the fly, and returning the average map.       **
 **	This routine differs from READFITS_AVG_MAP in that it uses a	    **
 **	companion list of 2-D images as an extra form of weighting in	    **
 **	addition to the "standard" set of scalar weights.  Such weighting   **
 **	is useful when combining non-intensity images, e.g., maps of dust   **
 **	temperature, which should average fairly well with 100um intensity  **
 **	maps as extra weights.						    **
 **									    **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void readfits_avg_map2 (num_in, fname_list, weight_list, weightmap_list, 
							avgmap_ptr, hpar_ptr)
  char fname_list[MAX_NAVG][MAX_PARSE_LENGTH+1], 
	weightmap_list[MAX_NAVG][MAX_PARSE_LENGTH+1];
  int num_in;
  float **avgmap_ptr, weight_list[];
  header_param_list *hpar_ptr;
{
  char object_old[LINELEN];
  int i, j, k, n1, n2;
  float *indata, *wdata, *tweight_map, weight;
  FILE *infile, *wfile;
  header_param_list hpar_w;

  for (k=0; k<num_in; k++) {
    printf ("  Averaging map %s with scalar weight %.4f and weight map %s\n", 
			fname_list[k], weight_list[k], weightmap_list[k]);
    open_read (&infile, fname_list[k]);
    readfits_header (infile, hpar_ptr);
    if ((*hpar_ptr).num_axes < 2) {
      fprintf (stderr, 
		"*** error in READFITS_AVG_MAP2: img num_axes = %d < 2 ***\n",
							(*hpar_ptr).num_axes);
      exit (1);
    }
    open_read (&wfile, weightmap_list[k]);
    readfits_header (wfile, &hpar_w);
    if ((*hpar_ptr).num_axes < 2) {
      fprintf (stderr, 
		"*** error in READFITS_AVG_MAP2: wght num_axes = %d < 2 ***\n",
							hpar_w.num_axes);
      exit (1);
    }
    if (k==0) {
      n1 = (*hpar_ptr).naxis[0];
      n2 = (*hpar_ptr).naxis[1];
      (*avgmap_ptr) = (float *) malloc (n1 * n2 * sizeof (float));
      tweight_map = (float *) malloc (n1 * n2 * sizeof (float));
      for (i=0; i<n1*n2; i++) {
	(*avgmap_ptr)[i] = BLANK_PIXEL;
	tweight_map[i] = 0.0;
      }
      indata = (float *) malloc (n1 * n2 * sizeof (float));
      wdata  = (float *) malloc (n1 * n2 * sizeof (float));
    } else {
      if (   n1 != (*hpar_ptr).naxis[0] || n2 != (*hpar_ptr).naxis[1]
	  || n1 != hpar_w.naxis[0]      || n2 != hpar_w.naxis[1]) {
	fprintf (stderr, 
	      "*** error in READFITS_AVG_MAP2: inconsistent dimensions ***\n");
	exit (1);
      }
    }
    readfits_plane (infile, indata, hpar_ptr);
    fclose (infile);
    readfits_plane (wfile, wdata, &hpar_w);
    fclose (wfile);
    for (j=0; j<n2; j++) {
      for (i=0; i<n1; i++) {
	if (indata[i+j*n1] < BLANK_THRESH) {
	  if (wdata[i+j*n1] < BLANK_THRESH) {
	    weight = weight_list[k] * wdata[i+j*n1];
	  } else {
	    weight = weight_list[k];
	  }
	  if ((*avgmap_ptr)[i+j*n1] >= BLANK_THRESH) {
	    (*avgmap_ptr)[i+j*n1] = weight * indata[i+j*n1];
	  } else {
	    (*avgmap_ptr)[i+j*n1] += weight * indata[i+j*n1];
	  }
	  tweight_map[i+j*n1] += weight;
	}
      }
    }
  }
  free (indata);
  free (wdata);
  for (i=0; i<n1*n2; i++) {
    if (tweight_map[i] != 0.0) (*avgmap_ptr)[i] /= tweight_map[i];
  }
  free (tweight_map);
  strcpy (object_old, (*hpar_ptr).object);
  sprintf ((*hpar_ptr).object, "Weighted Inverse Average %s", object_old);
}

/*****************************************************************************
 **									    **
 **		  R E A D F I T S _ I N V _ A V G _ M A P	            **
 **									    **
 **	Read list of 2-D FITS images in sequence, computing a weighted      **
 **     INVERSE average of the set on the fly, and return the average map.  **
 **									    **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void readfits_inv_avg_map (num_in, fname_list, weight_list, avgmap_ptr, 
								hpar_ptr)
  char fname_list[MAX_NAVG][MAX_PARSE_LENGTH+1];
  int num_in;
  float **avgmap_ptr, weight_list[];
  header_param_list *hpar_ptr;
{
  char object_old[LINELEN];
  int i, j, k, n1, n2;
  float *indata, *tweight_map;
  FILE *infile;

  for (k=0; k<num_in; k++) {
    printf ("  Inverse-averaging map %s with weight %.4f\n", 
						fname_list[k], weight_list[k]);
    open_read (&infile, fname_list[k]);
    readfits_header (infile, hpar_ptr);
    if ((*hpar_ptr).num_axes < 2) {
      fprintf (stderr, 
		"*** error in READFITS_INV_AVG_MAP: num_axes = %d < 2 ***\n",
							(*hpar_ptr).num_axes);
      exit (1);
    }
    if (k==0) {
      n1 = (*hpar_ptr).naxis[0];
      n2 = (*hpar_ptr).naxis[1];
      (*avgmap_ptr) = (float *) malloc (n1 * n2 * sizeof (float));
      tweight_map = (float *) malloc (n1 * n2 * sizeof (float));
      for (i=0; i<n1*n2; i++) {
	(*avgmap_ptr)[i] = BLANK_PIXEL;
	tweight_map[i] = 0.0;
      }
      indata = (float *) malloc (n1 * n2 * sizeof (float));
    } else {
      if (n1 != (*hpar_ptr).naxis[0] || n2 != (*hpar_ptr).naxis[1]) {
	fprintf (stderr, 
	   "*** error in READFITS_INV_AVG_MAP: inconsistent dimensions ***\n");
	exit (1);
      }
    }
    readfits_plane (infile, indata, hpar_ptr);
    fclose (infile);
    for (j=0; j<n2; j++) {
      for (i=0; i<n1; i++) {
	if (indata[i+j*n1] < BLANK_THRESH && indata[i+j*n1] != 0.0) {
	  if ((*avgmap_ptr)[i+j*n1] >= BLANK_THRESH) {
	    (*avgmap_ptr)[i+j*n1] = weight_list[k] / indata[i+j*n1];
	  } else {
	    (*avgmap_ptr)[i+j*n1] += weight_list[k] / indata[i+j*n1];
	  }
	  tweight_map[i+j*n1] += weight_list[k];
	}
      }
    }
  }
  free (indata);
  for (i=0; i<n1*n2; i++) {
    if (tweight_map[i] != 0.0) {
      (*avgmap_ptr)[i] = 1.0 / ((*avgmap_ptr)[i] / tweight_map[i]);
    }
  }
  free (tweight_map);
  strcpy (object_old, (*hpar_ptr).object);
  sprintf ((*hpar_ptr).object, "Weighted Inverse Average %s", object_old);
}

/*****************************************************************************
 **									    **
 **		   R E A D F I T S _ I N V _ A V G _ M A P 2	            **
 **									    **
 **	Read list of 2-D FITS images in sequence, computing a weighted      **
 **     INVERSE average of the set on the fly, and return the average map.  **
 **	This routine differs from READFITS_INV_AVG_MAP in that it uses a    **
 **	companion list of 2-D images as an extra form of weighting in	    **
 **	addition to the "standard" set of scalar weights.  Such weighting   **
 **	is useful when combining non-intensity images, e.g., maps of dust   **
 **	temperature, which should average fairly well with 100um intensity  **
 **	maps as extra weights.						    **
 **									    **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void readfits_inv_avg_map2 (num_in, fname_list, weight_list, weightmap_list, 
							avgmap_ptr, hpar_ptr)
  char fname_list[MAX_NAVG][MAX_PARSE_LENGTH+1], 
	weightmap_list[MAX_NAVG][MAX_PARSE_LENGTH+1];
  int num_in;
  float **avgmap_ptr, weight_list[];
  header_param_list *hpar_ptr;
{
  char object_old[LINELEN];
  int i, j, k, n1, n2;
  float *indata, *wdata, *tweight_map, weight;
  FILE *infile, *wfile;
  header_param_list hpar_w;

  for (k=0; k<num_in; k++) {
    printf ("  Averaging map %s with scalar weight %.4f and weight map %s\n", 
			fname_list[k], weight_list[k], weightmap_list[k]);
    open_read (&infile, fname_list[k]);
    readfits_header (infile, hpar_ptr);
    if ((*hpar_ptr).num_axes < 2) {
      fprintf (stderr, 
	"*** error in READFITS_INV_AVG_MAP2: img num_axes = %d < 2 ***\n",
							(*hpar_ptr).num_axes);
      exit (1);
    }
    open_read (&wfile, weightmap_list[k]);
    readfits_header (wfile, &hpar_w);
    if ((*hpar_ptr).num_axes < 2) {
      fprintf (stderr, 
	"*** error in READFITS_INV_AVG_MAP2: wght num_axes = %d < 2 ***\n",
							hpar_w.num_axes);
      exit (1);
    }
    if (k==0) {
      n1 = (*hpar_ptr).naxis[0];
      n2 = (*hpar_ptr).naxis[1];
      (*avgmap_ptr) = (float *) malloc (n1 * n2 * sizeof (float));
      tweight_map = (float *) malloc (n1 * n2 * sizeof (float));
      for (i=0; i<n1*n2; i++) {
	(*avgmap_ptr)[i] = BLANK_PIXEL;
	tweight_map[i] = 0.0;
      }
      indata = (float *) malloc (n1 * n2 * sizeof (float));
      wdata  = (float *) malloc (n1 * n2 * sizeof (float));
    } else {
      if (   n1 != (*hpar_ptr).naxis[0] || n2 != (*hpar_ptr).naxis[1]
	  || n1 != hpar_w.naxis[0]      || n2 != hpar_w.naxis[1]) {
	fprintf (stderr, 
	  "*** error in READFITS_INV_AVG_MAP2: inconsistent dimensions ***\n");
	exit (1);
      }
    }
    readfits_plane (infile, indata, hpar_ptr);
    fclose (infile);
    readfits_plane (wfile, wdata, &hpar_w);
    fclose (wfile);
    for (j=0; j<n2; j++) {
      for (i=0; i<n1; i++) {
	if (indata[i+j*n1] < BLANK_THRESH && indata[i+j*n1] != 0.0) {
	  if (wdata[i+j*n1] < BLANK_THRESH && wdata[i+j*n1] != 0.0) {
	    weight = weight_list[k] / wdata[i+j*n1];
	  } else {
	    weight = weight_list[k];
	  }
	  if ((*avgmap_ptr)[i+j*n1] >= BLANK_THRESH) {
	    (*avgmap_ptr)[i+j*n1] = weight / indata[i+j*n1];
	  } else {
	    (*avgmap_ptr)[i+j*n1] += weight / indata[i+j*n1];
	  }
	  tweight_map[i+j*n1] += weight;
	}
      }
    }
  }
  free (indata);
  free (wdata);
  for (i=0; i<n1*n2; i++) {
    if (tweight_map[i] != 0.0) {
      (*avgmap_ptr)[i] = 1.0 / ((*avgmap_ptr)[i] / tweight_map[i]);
    }
  }
  free (tweight_map);
  strcpy (object_old, (*hpar_ptr).object);
  sprintf ((*hpar_ptr).object, "Weighted Inverse Average %s", object_old);
}

/*****************************************************************************
 **									    **
 **			 W R I T E F I T S _ M A P		            **
 **									    **
 **	Write 2-dimensional FITS image and its header.                      **
 **                                                                         **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void writefits_map (fname, data, hpar_ptr)
  const char *fname;
  float *data;
  header_param_list *hpar_ptr;
{
  int data_bytes, dblocks, pad_bytes;
  FILE *outfile;

  if ((*hpar_ptr).num_axes != 2) {
    fprintf (stderr, "*** WARNING in WRITEFITS_MAP: num_axes = %d != 2 ***\n",
							(*hpar_ptr).num_axes);
  }
  open_write (&outfile, fname);
  writefits_header (outfile, hpar_ptr);
  writefits_plane (outfile, data, hpar_ptr);
  data_bytes = (*hpar_ptr).naxis[0] * (*hpar_ptr).naxis[1] 
						* abs ((*hpar_ptr).bitpix / 8);
  dblocks    = (int) ceil (1.0 * data_bytes / BLOCKSIZE);
  pad_bytes  = dblocks * BLOCKSIZE - data_bytes;
  writefits_pad (outfile, pad_bytes);
  fclose (outfile);
}

/*****************************************************************************
 **									    **
 **			 R E A D F I T S _ C U B E		            **
 **									    **
 **	Read 3-dimensional FITS image (of floats) and its header.           **
 **                                                                         **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void readfits_cube (fname, data_ptr, hpar_ptr, verbose_flag)
  const char *fname;
  int verbose_flag;
  float **data_ptr;
  header_param_list *hpar_ptr;
{
  int i;
  FILE *infile;

  open_read (&infile, fname);
  readfits_header (infile, hpar_ptr);
  if ((*hpar_ptr).num_axes < 3) {
    fprintf (stderr, "*** error in READFITS_CUBE: num_axes = %d < 3 ***\n",
							(*hpar_ptr).num_axes);
    exit (1);
  }
  *data_ptr = (float *) malloc (  (*hpar_ptr).naxis[0] * (*hpar_ptr).naxis[1] 
				* (*hpar_ptr).naxis[2] * sizeof (float));
  for (i=0; i<(*hpar_ptr).naxis[2]; i++) {
    if (verbose_flag == TRUE) {
      printf ("\tReading plane %d of %d.\n%c[1A", i+1,(*hpar_ptr).naxis[2],27);
    }
    readfits_plane (infile, &((*data_ptr)[i
		* (*hpar_ptr).naxis[0] * (*hpar_ptr).naxis[1]]), hpar_ptr);
  }
  fclose (infile);
}

/*****************************************************************************
 **									    **
 **			 W R I T E F I T S _ C U B E		            **
 **									    **
 **	Write 3-dimensional FITS image and its header.                      **
 **                                                                         **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void writefits_cube (fname, data, hpar_ptr, verbose_flag)
  const char *fname;
  int verbose_flag;
  float *data;
  header_param_list *hpar_ptr;
{
  int i, data_bytes, dblocks, pad_bytes;
  FILE *outfile;

  if ((*hpar_ptr).num_axes < 3) {
    fprintf (stderr, "*** error in WRITEFITS_CUBE: num_axes = %d < 3 ***\n",
							(*hpar_ptr).num_axes);
    exit (1);
  }
  if ((*hpar_ptr).num_axes > 3) {
    fprintf (stderr, "*** WARNING in WRITEFITS_CUBE: num_axes = %d > 3 ***\n",
							(*hpar_ptr).num_axes);
  }
  open_write (&outfile, fname);
  writefits_header (outfile, hpar_ptr);
  for (i=0; i<(*hpar_ptr).naxis[2]; i++) {
    if (verbose_flag == TRUE) {
      printf ("\tWriting plane %d of %d.\n%c[1A", i+1,(*hpar_ptr).naxis[2],27);
  }
    writefits_plane (outfile, &(data[i
		* (*hpar_ptr).naxis[0] * (*hpar_ptr).naxis[1]]), hpar_ptr);
  }
  data_bytes = (*hpar_ptr).naxis[0] * (*hpar_ptr).naxis[1] 
			* (*hpar_ptr).naxis[2] * abs ((*hpar_ptr).bitpix / 8);
  dblocks    = (int) ceil (1.0 * data_bytes / BLOCKSIZE);
  pad_bytes  = dblocks * BLOCKSIZE - data_bytes;
  writefits_pad (outfile, pad_bytes);
  fclose (outfile);
}

/*****************************************************************************
 **									    **
 **			 R E A D F I T S _ 4 C U B E		            **
 **									    **
 **	Read 4-dimensional FITS image (of floats) and its header.           **
 **                                                                         **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void readfits_4cube (fname, data_ptr, hpar_ptr)
  const char *fname;
  float **data_ptr;
  header_param_list *hpar_ptr;
{
  int i, j, n1, n2, n3, n4;
  FILE *infile;

  open_read (&infile, fname);
  readfits_header (infile, hpar_ptr);
  if ((*hpar_ptr).num_axes < 4) {
    fprintf (stderr, "*** error in READFITS_4CUBE: num_axes = %d < 4 ***\n",
							(*hpar_ptr).num_axes);
    exit (1);
  }
  n1 = (*hpar_ptr).naxis[0];  n2 = (*hpar_ptr).naxis[1];
  n3 = (*hpar_ptr).naxis[2];  n4 = (*hpar_ptr).naxis[3];
  *data_ptr = (float *) malloc (n1 * n2 * n3 * n4 * sizeof (float));
  for (j=0; j<n4; j++) {
    for (i=0; i<n3; i++) {
      printf ("\tReading plane [%d,%d] of [%d,%d].\n%c[1A", 
						i+1, j+1, n3, n4, 27);
      readfits_plane (infile, &((*data_ptr)[i*n1*n2+j*n1*n2*n3]), hpar_ptr);
    }
  }
  fclose (infile);
}

/*****************************************************************************
 **									    **
 **			 W R I T E F I T S _ 4 C U B E		            **
 **									    **
 **	Write 4-dimensional FITS image and its header.                      **
 **                                                                         **
 **	This is a "shell" routine which calls several others.               **
 **									    **
 *****************************************************************************/
void writefits_4cube (fname, data, hpar_ptr)
  const char *fname;
  float *data;
  header_param_list *hpar_ptr;
{
  int i, j, n1, n2, n3, n4, data_bytes, dblocks, pad_bytes;
  FILE *outfile;

  if ((*hpar_ptr).num_axes < 4) {
    fprintf (stderr, "*** error in WRITEFITS_4CUBE: num_axes = %d < 4 ***\n",
							(*hpar_ptr).num_axes);
    exit (1);
  }
  if ((*hpar_ptr).num_axes > 4) {
    fprintf (stderr, "*** WARNING in WRITEFITS_4CUBE: num_axes = %d > 4 ***\n",
							(*hpar_ptr).num_axes);
  }
  n1 = (*hpar_ptr).naxis[0];  n2 = (*hpar_ptr).naxis[1];
  n3 = (*hpar_ptr).naxis[2];  n4 = (*hpar_ptr).naxis[3];
  open_write (&outfile, fname);
  writefits_header (outfile, hpar_ptr);
  for (j=0; j<n4; j++) {
    for (i=0; i<n3; i++) {
      printf ("\tWriting plane [%d,%d] of [%d,%d].\n%c[1A", 
					i+1, j+1, n3, n4, 27);
      writefits_plane (outfile, &(data[i*n1*n2+j*n1*n2*n3]), hpar_ptr);
    }
  }
  data_bytes = n1*n2*n3*n4 * abs ((*hpar_ptr).bitpix / 8);
  dblocks    = (int) ceil (1.0 * data_bytes / BLOCKSIZE);
  pad_bytes  = dblocks * BLOCKSIZE - data_bytes;
  writefits_pad (outfile, pad_bytes);
  fclose (outfile);
}

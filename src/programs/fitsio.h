/*############################################################################
 #############################################################################
 ##                                                                         ##
 ##			       F I T S I O . h                                          ##
 ##                                                                         ##
 ##	Source code for a few FITS I/O routines, to be                          ##
 ##	#included in various data reduction codes.                              ##
 ##                                                                         ##
 ##	consolidated Monday 29 May 1995 sjg                                     ##
 ## converted to a library with header February 2006 jsd                    ##
 ##                                                                         ##
 #############################################################################
 ############################################################################*/


#include <stdio.h>
#include "mathdefs.h"
#include "io.h"

#ifndef _FITSIO_H
#define _FITSIO_H

#define NIL            		0	/* empty pointer value */
#define BLOCKSIZE   	     2880	/* READFITS-WRITEFITS: bytes */
#define LINELEN      	       80	/* READFITS-WRITEFITS: bytes */
#define LINESPERBLOCK	       36	/* READFITS-WRITEFITS: 2880/80 */
#define VMS         	    FALSE	/* READFITS-WRITEFITS: VAX/VMS C? */
#define VERBOSE		    FALSE	/* READFITS: debugging */
#define MAX_NUM_AXES		7	/* maximum number of axes supported */
#define MAX_COMMENT	     1000	/* maximum number of comment lines */
#define MAX_HISTORY	     1000	/* maximum number of history lines */
#define BLANK_PIXEL		1.0e36  /* internal storage value (float) */
#define BLANK_THRESH		5.0e35  /* threshhold for flagging as blank */
#define WARN_BLANK_CHANGE   FALSE	/* announce change of BLANK value? */
#define MAX_NAVG	     1024	/* max # of maps to average together */
#define MAX_PARSE_LENGTH 	1024	/* maximum line length for parsing */


#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )


typedef struct subset_param_list {
  int num_axes, start[MAX_NUM_AXES], stop[MAX_NUM_AXES], stride[MAX_NUM_AXES];
} subset_param_list;

typedef struct header_param_list {
  char ctype[MAX_NUM_AXES][LINELEN], bunit[LINELEN],
	date[LINELEN], object[LINELEN], origin[LINELEN], 
	telescope[LINELEN], instrument[LINELEN],
	observer[LINELEN], filename[LINELEN], 
	comment[MAX_COMMENT][LINELEN], history[MAX_HISTORY][LINELEN];
  int simple, blocked, extend,
	bitpix, num_axes, naxis[MAX_NUM_AXES], blank, ncomment, nhistory;
  double bscale, bzero, 
	datamin, datamax, obsfreq, exptime, epoch, equinox,
	crval[MAX_NUM_AXES], crpix[MAX_NUM_AXES], cdelt[MAX_NUM_AXES], 
	crota[MAX_NUM_AXES];
} header_param_list;


void printhex (char *startaddr, int numbytes, int endline);

void init_subset_param_list (subset_param_list *spar_ptr);

void init_header_param_list (header_param_list *hpar_ptr);

void copy_header_param_list (header_param_list *hpar2_ptr, header_param_list *hpar1_ptr);

double wcs_from_pix (header_param_list *hpar_ptr, int axis, double pixel);

double pix_from_wcs (header_param_list *hpar_ptr, int axis, double coord);

void set_bzero_bscale (header_param_list *hpar_ptr, double min, double max);

void set_bzero_bscale_datarange (header_param_list *hpar_ptr, float *data);

void get_fits_str_par (char *line, char *str);

void readfits_header (FILE *file_ptr, header_param_list *hpar_ptr);

void writefits_header (FILE *file_ptr, header_param_list *hpar_ptr);

void readfits_plane (FILE *file_ptr, float* data, header_param_list *hpar_ptr);

void writefits_plane (FILE *file_ptr, float* data, header_param_list *hpar_ptr);

void writefits_pad (FILE *file_ptr, int nbytes);

void writefits_pad_end (FILE *file_ptr, header_param_list *hpar_ptr);

void readfits_plane_allocate (FILE *file_ptr, float **data_ptr, header_param_list *hpar_ptr);

void readfits_map (char *fname, float **data_ptr, header_param_list *hpar_ptr);

void readfits_avg_map (int num_in, char fname_list[MAX_NAVG][MAX_PARSE_LENGTH+1], float weight_list[], float **avgmap_ptr, header_param_list *hpar_ptr);

void readfits_avg_map2 (int num_in, char fname_list[MAX_NAVG][MAX_PARSE_LENGTH+1], float weight_list[], char weightmap_list[MAX_NAVG][MAX_PARSE_LENGTH+1], float **avgmap_ptr, header_param_list *hpar_ptr);

void writefits_map (const char *fname, float *data, header_param_list *hpar_ptr);

void readfits_cube (const char *fname, float **data_ptr, header_param_list *hpar_ptr, int verbose_flag);

void writefits_cube (const char *fname, float *data, header_param_list *hpar_ptr, int verbose_flag);

void readfits_4cube (const char *fname, float **data_ptr, header_param_list *hpar_ptr);

void writefits_4cube (const char *fname, float *data, header_param_list *hpar_ptr);

#endif

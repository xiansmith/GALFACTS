/*
 * Concatenates fits maps into a fits cube.
 * 
 * The values to be writtin into the fits header for the third (Z) axis
 * must be specified at the command line.
 * The value of crpix is assumed to be zero.
 * 
 * file1.fits ... fileN.fits - list of maps (eg *.fits)
 * cube.fits - the name of the cube file to create
 * title - the title of the fits cube
 * Z ctype - the axis label
 * Z crval - the value for the first plane
 * Z cdelt - the increment for subsequent planes (can be negative)
 * 
 * Jeff Dever
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsio.h"


static void print_usage(char * prog)
{
	printf("Usage: %s <file1.fits> <file2.fits> ... <fileN.fits> <cube.fits> <title> <Z ctype> <Z crval> <Z cdelt>\n", prog);
}

int main(int argc, char *argv[]) 
{
	int i;
	char * cube_name;
	char * avg_name;
	FILE * cube_file;
	header_param_list cube_hpar;
	char * cube_title;
	FILE * tmpfile;
	int num_files;
	char * zctype;
	float zcrval;
	float zcdelt;

	/* process the arguments */
	if (argc < 6) {
		return EXIT_FAILURE;
	}
	num_files = argc - 6;
	cube_name = argv[argc-5];
	cube_title = argv[argc-4];
	zctype = argv[argc-3];
	zcrval = strtof(argv[argc-2], NULL);
	zcdelt = strtof(argv[argc-1], NULL);

	/* print arguments after processing so the user knows */
	printf("cube_title: %s\n", cube_title);
	printf("cube_name: %s\n", cube_name);
	printf("num_files: %i\n", num_files);
	printf("zctype: %s\n", zctype);
	printf("zcrval: %f\n", zcrval);
	printf("zcdelt: %f\n", zcdelt);

	/* quick check to sanitize the argument processing */
	if (!isfinite(zcrval) || !isfinite(zcdelt)) {
		printf("ERROR: invalid arguments\n");
		print_usage(argv[0]);
		return EXIT_FAILURE;
	} 

	/* do a correction so the start plane actually starts at zcrval */
	zcrval -= zcdelt;

	/* open the cube to create from the files */
	cube_file = fopen(cube_name, "w");
	if (cube_file == NULL) {
		printf("ERROR: unable to open '%s' for writing\n", cube_name);
		return EXIT_FAILURE;
	}
	
	/* read the header from the first file to use as the template */
	tmpfile = fopen(argv[1], "r");
	readfits_header (tmpfile, &cube_hpar);
	fclose(tmpfile);
	cube_hpar.num_axes = 3;
	cube_hpar.naxis[2] = num_files;
	cube_hpar.crval[2] = zcrval;
	cube_hpar.crpix[2] = 0.0;
	cube_hpar.cdelt[2] = zcdelt;
	sprintf (cube_hpar.ctype[2], zctype);

	strcpy(cube_hpar.object, cube_title);
	writefits_header (cube_file, &cube_hpar);

	//read the maps and write a planes of the cube
	for (i=1; i<=num_files; i++) 
	{
		char * fname = argv[i];
		float * data;
		header_param_list hpar;
		readfits_map (fname, &data, &hpar);
		writefits_plane (cube_file, data, &cube_hpar);
		free(data);
	}
	
	/* close and clean up */
	writefits_pad_end (cube_file, &cube_hpar);
	fclose(cube_file);
	
	return EXIT_SUCCESS;
}



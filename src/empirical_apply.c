#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "programs/chardefs.h"
#include "programs/mathdefs.h"
#include "programs/misc.c"
#include "programs/misc_math.c"
#include "programs/coord_utils.c"
#include "programs/fitsio.c"


#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )


void create_outfilename(char * outfilename, char * infilename)
{
	char * ext;

	ext = strrchr(infilename, '.');
	if (ext == '\0') {
		printf("ERROR: no '.' in '%s'\n", infilename);
		return;
	}
	outfilename[0] = '\0';
	strncat(outfilename, infilename, ext-infilename);
	strcat(outfilename, "_empirical");
	strcat(outfilename, ext);
}

int main(int argc, char * argv[])
{

	header_param_list Ihpar, Qhpar, Uhpar, Vhpar;
	header_param_list Iouthpar, Qouthpar, Uouthpar, Vouthpar;
	FILE *Ifile, *Qfile, *Ufile, *Vfile, *datfile;
	FILE *Ioutfile, *Qoutfile, *Uoutfile, *Voutfile;
	char Ioutfilename[128], Qoutfilename[128], Uoutfilename[128], Voutfilename[128];
	char *Ifilename, *Qfilename, *Ufilename, *Vfilename;
	int num_planes, num_points;
	int i;
	char *datfilename;
	float *Iplane, *Qplane, *Uplane, *Vplane;

	if (argc != 6) {
		printf("usage: %s <I.fits> <Q.fits> <U.fits> <V.fits> <empirical.dat>\n", argv[0]);
		return EXIT_FAILURE;
	}
	
	Ifilename = argv[1];
	Qfilename = argv[2];
	Ufilename = argv[3];
	Vfilename = argv[4];
	datfilename = argv[5];

	//open files for reading
	Ifile = fopen(Ifilename, "r");
	Qfile = fopen(Qfilename, "r");
	Ufile = fopen(Ufilename, "r");
	Vfile = fopen(Vfilename, "r");
	datfile = fopen(datfilename, "r");
	if (Ifile == NULL || Qfile == NULL || Ufile == NULL || Vfile == NULL || datfile == NULL) {
		printf("ERROR: unable to open input file '%s' for reading\n", Ifilename);
		return EXIT_FAILURE;
	}

	//open files for writing
	create_outfilename(Ioutfilename, Ifilename); 
	Ioutfile = fopen(Ioutfilename, "w");
	create_outfilename(Qoutfilename, Qfilename); 
	Qoutfile = fopen(Qoutfilename, "w");
	create_outfilename(Uoutfilename, Ufilename); 
	Uoutfile = fopen(Uoutfilename, "w");
	create_outfilename(Voutfilename, Vfilename); 
	Voutfile = fopen(Voutfilename, "w");
	if (Ioutfile == NULL || Qoutfile == NULL || Uoutfile == NULL || Voutfile == NULL) {
		printf("ERROR: unable to open output file '%s' for writing\n", Ioutfilename);
		return EXIT_FAILURE;
	}

	//read the headers
	readfits_header (Ifile, &Ihpar);
	readfits_header (Qfile, &Qhpar);
	readfits_header (Ufile, &Uhpar);
	readfits_header (Vfile, &Vhpar);
	
	if ((memcmp(Ihpar.naxis, Qhpar.naxis, 3) != 0) ||
		(memcmp(Ihpar.naxis, Uhpar.naxis, 3) != 0) ||
		(memcmp(Ihpar.naxis, Vhpar.naxis, 3) != 0) ) {
		printf("ERROR: naxis size mismatch between fits headers\n");
		return EXIT_FAILURE;
	}

	num_points = Ihpar.naxis[0] * Ihpar.naxis[1];
	num_planes = Ihpar.naxis[2];
	printf("num_points: %i\n", num_points);
	printf("num_planes: %i\n", num_planes);

	Iplane = malloc(sizeof(float) * num_points);
	Qplane = malloc(sizeof(float) * num_points);
	Uplane = malloc(sizeof(float) * num_points);
	Vplane = malloc(sizeof(float) * num_points);

	Iouthpar = Ihpar;
	Qouthpar = Qhpar;
	Uouthpar = Uhpar;
	Vouthpar = Vhpar;

	strcat(Iouthpar.object, " Empirical");
	strcat(Qouthpar.object, " Empirical");
	strcat(Uouthpar.object, " Empirical");
	strcat(Vouthpar.object, " Empirical");

	writefits_header(Ioutfile, Iouthpar);
	writefits_header(Qoutfile, Qouthpar);
	writefits_header(Uoutfile, Uouthpar);
	writefits_header(Voutfile, Vouthpar);

	for (i=0; i<num_planes; i++) 
	{
		int j, ix;
		float deltaInot, deltaQnot, deltaUnot, deltaVnot;
		float Inot, Qnot, Unot, Vnot;
		float mI, mQ, mU, mV;
		float Imax, Qmax, Umax, Vmax;

		//read a line from the empirical.dat file
		fscanf(datfile, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", 
			&ix, &mQ, &mU, &mV, 
			&Inot, &Qnot, &Unot, &Vnot, 
			&deltaInot, &deltaQnot, &deltaUnot, &deltaVnot, 
			&Imax, &Qmax, &Umax, &Vmax);
printf("mQ: %f  mU: %f  mV: %f\n", mQ, mU, mV);
		readfits_plane (Ifile, Iplane, &Ihpar);
		readfits_plane (Qfile, Qplane, &Qhpar);
		readfits_plane (Ufile, Uplane, &Uhpar);
		readfits_plane (Vfile, Vplane, &Vhpar);

		//do the actual corrections
		for (j=0; j<num_points; j++) {
			Iplane[j] = Iplane[j] - Inot;	
			Qplane[j] = Qplane[j] - Qnot - (mQ * Iplane[j]);	
			Uplane[j] = Uplane[j] - Unot - (mU * Iplane[j]);	
			Vplane[j] = Vplane[j] - Vnot - (mV * Iplane[j]);	
		}

		writefits_plane(Ioutfile, Iplane, &Iouthpar);
		writefits_plane(Qoutfile, Qplane, &Qouthpar);
		writefits_plane(Uoutfile, Uplane, &Uouthpar);
		writefits_plane(Voutfile, Vplane, &Vouthpar);
	}

	fclose(Ifile);
	fclose(Qfile);
	fclose(Ufile);
	fclose(Vfile);
	fclose(datfile);

	free(Iplane);
	free(Qplane);
	free(Uplane);
	free(Vplane);

	writefits_pad_end(Ioutfile, Iouthpar);
	writefits_pad_end(Qoutfile, Qouthpar);
	writefits_pad_end(Uoutfile, Uouthpar);
	writefits_pad_end(Voutfile, Vouthpar);

	fclose(Ioutfile);
	fclose(Qoutfile);
	fclose(Uoutfile);
	fclose(Voutfile);

}



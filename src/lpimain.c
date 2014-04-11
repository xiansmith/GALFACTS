#include <stdio.h>
#include <stdlib.h>
#include "programs/chardefs.h"
#include "programs/mathdefs.h"
#include "programs/misc.c"
#include "programs/misc_math.c"
#include "programs/coord_utils.c"
#include "programs/fitsio.c"

#include "common.h"


#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )



int main(int argc, char * argv[])
{
	header_param_list Qhpar, Uhpar, PIhpar;
	char * Qfilename, *Ufilename, *PIfilename;
	int num_points;
	float *Qdata, *Udata, *PIdata;
	int i;

	if (argc != 4) {
		printf("usage: %s <Qavg.fits> <Uavg.fits> <PIout.fits\n", argv[0]);
		return EXIT_FAILURE;
	}
	
	Qfilename = argv[1];
	Ufilename = argv[2];
	PIfilename = argv[3];

	//read average maps

	readfits_map (Qfilename, &Qdata, &Qhpar);
	readfits_map (Ufilename, &Udata, &Uhpar);
	if (memcmp(Qhpar.naxis, Uhpar.naxis, 2) != 0) {
		printf("ERROR: naxis size mismatch between fits headers\n");
		return EXIT_FAILURE;
	}
	num_points = Qhpar.naxis[0] * Uhpar.naxis[1];

	PIdata = (float*) malloc(sizeof(float) * num_points);

	for (i=0; i<num_points; i++) 
	{
		if (IS_BLANK_PIXEL(Qdata[i]) || IS_BLANK_PIXEL(Udata[i])) {
			PIdata[i] = BLANK_PIXEL;
		} else {
			PIdata[i] = sqrt(Qdata[i]*Qdata[i]+Udata[i]*Udata[i]);
		}
	}

	PIhpar = Qhpar;
	//TODO: change the name in the header
	writefits_map (PIfilename, PIdata, &PIhpar);

	free(Udata);
	free(Qdata);
	free(PIdata);

}



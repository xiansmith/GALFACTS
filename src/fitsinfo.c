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


int main(int argc, char *argv[]) 
{

	char * filename;
	FILE * file;
	header_param_list cube_hpar;


	if (argc < 3) {
		printf("Usage: %s <file.fits>\n", argv[0]);
		return EXIT_FAILURE;
	}
	filename = argv[1];


	file = fopen(filename, "r");
	readfits_header(cube_file, &cube_hpar);

	return EXIT_SUCCESS;
}

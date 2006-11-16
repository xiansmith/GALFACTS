#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "programs/chardefs.h"
#include "programs/mathdefs.h"
#include "programs/misc.c"
#include "programs/misc_math.c"
#include "programs/coord_utils.c"
#include "programs/fitsio.c"
#include <float.h>

#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )

#define MAX_SOURCES 128

typedef struct {
	int peakX;
	int peakY;
	int backX;
	int backY;
} SourceData;

int read_coordinate_file(const char * filename, SourceData data[], int maxsize)
{
	FILE * file;
	int i, num;

	file = fopen(filename, "r");
	if (file == NULL) {
		printf("ERROR: unable to open input file %s\n", filename);
		return -1;
	}

	i = 0;
	do {
		num = fscanf(file, "%i %i %i %i\n", &data[i].peakX, &data[i].peakY, &data[i].backX, &data[i].backY);
		if (num == 4) i++;
	} while (!feof(file) && i<maxsize);

	fclose(file);
	return i;
}

int main(int argc, char * argv[])
{
	int i, j, s;
	header_param_list hpar;
	float *cube;
	int num_x, num_y;
	int num_planes;
	FILE * outfile;
	SourceData sourcedata[MAX_SOURCES];
	int num_sources;

	if (argc != 4) {
		printf("usage: %s <in.dat> <in.fits> <out.dat>\n", argv[0]);
		return EXIT_FAILURE;
	}

	num_sources = read_coordinate_file(argv[1], sourcedata, MAX_SOURCES);
	readfits_cube (argv[2], &cube, &hpar);
	outfile = fopen(argv[3], "w");

	num_x = hpar.naxis[0];
	num_y = hpar.naxis[1];
	num_planes = hpar.naxis[2];
	printf("num_x: %i\n", num_x);
	printf("num_y: %i\n", num_y);
	printf("num_planes: %i\n", num_planes);
	
	for (s=0; s<num_sources; s++) 
	{
		SourceData source = sourcedata[s];
		int count = 0;
		float avg;
		float min = FLT_MAX;
		float max = FLT_MIN;

		for (j=0; j<num_planes; j++) 
		{
			int pospeak = j*num_x*num_y + source.peakY*num_x + source.peakX;
			int posback = j*num_x*num_y + source.backY*num_x + source.backX;
			float valback = cube[posback];
			float valpeak = cube[pospeak] - valback;

			if (j>44 && j<49) continue; //skip HI channels
			if (IS_BLANK_PIXEL(valpeak) || IS_BLANK_PIXEL(valback))  continue;

if (j == 0) {
printf("channel: %i valpeak: %f valback: %f\n", j, valpeak, valback);
}

			if (valpeak < min) {
				min = valpeak;
			} 
			if (valpeak > max) {
				max = valpeak;
			}
			
			avg += valpeak;
			count++;
		}
		avg /= count;
		
		fprintf(outfile, "%f %f %f %f\n", max, min, max-min, avg);
	}

	fclose(outfile);
	free(cube);
}



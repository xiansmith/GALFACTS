#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsio.h"


static void select_ranges(char * params, int frames[])
{
	int i;
	char *toka, *tokb;

printf("params: %s\n", params);
	for (i=0; i<256; i++) {
		frames[i] = 0;
	}

	//split the params based on comma
	//process each range
	toka = strtok(params, ",");
	do
	{
		int a, b;
		if ( (tokb = strchr(toka, '-')) ){
			*tokb = '\0';
			tokb++;
			a = atoi(toka);
			b = atoi(tokb);
		} else {
			a = b = atoi(toka);
		}

printf("a=%i b=%i\n", a, b);

		//set the flags in frames
		for (i=a; i<=b; i++) {
			frames[i] = 1;
		}
	} while ( (toka = strtok(NULL, ",")) );
}

int main(int argc, char *argv[]) 
{

	int i,j;
	char * cube_name;
	char * avg_name;
	FILE * cube_file;
	header_param_list cube_hpar;
	int data_len;
	float * plane_data;
	float * avg_data;
	float * count_data;
	int frames[256];


	if (argc < 3) {
		printf("Usage: %s <cube.fits> <avg.fits> <range1,range2,...>\n", argv[0]);
		return EXIT_FAILURE;
	}
	cube_name = argv[1];
	avg_name = argv[2];
	if (argc == 3) {
		for (i=0; i<256; i++) {
			frames[i] = 1;
		}
	} else {
		select_ranges(argv[3], frames);
	}


	cube_file = fopen(cube_name, "r");
	readfits_header(cube_file, &cube_hpar);

	//allocate memory for plane
	data_len = cube_hpar.naxis[0] * cube_hpar.naxis[1];
	plane_data = (float*) calloc(data_len, sizeof (float));
	avg_data = (float*) calloc(data_len, sizeof (float));
	count_data = (float*) calloc(data_len, sizeof (float));

	//iterate over the cube
	for (i=0; i<cube_hpar.naxis[2]; i++) 
	{
		//read a plane
		readfits_plane(cube_file, plane_data, &cube_hpar);
		//include in average only if stated in frames
		if (frames[i] == 1) {
			//printf("%i \n", i);
			for (j=0; j<data_len; j++) {
				//add to average and increment count if its not blank
				if (!IS_BLANK_PIXEL(plane_data[j]) && isfinite(plane_data[j])) {
					avg_data[j] += plane_data[j];
					count_data[j]++;
				}
			}
		}
	}
	free(plane_data);

	//calculate average
	for (j=0; j<data_len; j++) {
		if (count_data[j] == 0) {
			avg_data[j] = BLANK_PIXEL;
		} else {
			avg_data[j] /= count_data[j];
		}
	}
	free(count_data);
	
	//write the average fits
	cube_hpar.num_axes = 2;
	strcat(cube_hpar.object, " Average");
	writefits_map(avg_name, avg_data, &cube_hpar);
	free(avg_data);

	return EXIT_SUCCESS;
}

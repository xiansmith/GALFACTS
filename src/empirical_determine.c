#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "programs/fitsio.h"


#define MAX_NUM_PLANES 4096

typedef char byte;

/* Takes a comma seperated list of ranges, and sets flags in the 
 * planes array to specify if a plane is in one of the selected 
 * ranges.  Ranges can overlap and are inclusive.
 * All are selected if params is NULL;
 * eg: 1-3,5-6,8,9 would give [0,1,1,1,0,1,1,0,1,1,0,0,...]
 *
 * params- the ranges in string format
 * planes - set to 1 if plane of the cube is selected
 */
static void select_ranges(char * params, byte planes[])
{
	int i;
	char *toka, *tokb;

	printf("Ranges: %s\n", params);

	//handle special case where params are null
	if (params == NULL) {
		for (i=0; i<MAX_NUM_PLANES; i++) {
			planes[i] = 1;
		}
		return;
	}

	for (i=0; i<MAX_NUM_PLANES; i++) {
		planes[i] = 0;
	}

	//split the params based on comma
	//process each range
	toka = strtok(params, ",");
	do
	{
		int a, b;
		if (tokb = strchr(toka, '-')) {
			*tokb = '\0';
			tokb++;
			a = atoi(toka);
			b = atoi(tokb);
		} else {
			a = b = atoi(toka);
		}

		//printf("a=%i b=%i\n", a, b);

		//set the flags in planes
		for (i=a; i<=b; i++) {
			planes[i] = 1;
		}
	} while (toka = strtok(NULL, ","));
}

//average ignores blank pixels
double compute_average(float data[], int size)
{
	double sum;
	int count;
	int i;

	sum = 0.0;
	count = 0;
	for (i=0; i<size; i++) {
		if (!IS_BLANK_PIXEL(data[i])) {
			sum += data[i];
			count++;
		}
	}

	return sum/count;
}

//blank out pixels based on the refdata
void blank_corresponding_pixels(float refdata[], float data[], int size)
{
	int i;

	for (i=0; i<size; i++) {
		if (IS_BLANK_PIXEL(refdata[i])) {
			data[i] = BLANK_PIXEL;
		}
	}
}


//drop pixels if they don't meet the criteria
void blank_pixels(float data[], int size)
{
	double avg;
	int i;

	avg = compute_average(data, size);
	for (i=0; i<size; i++) {
		if (data[i] > (avg * 0.995)) {
			data[i] = BLANK_PIXEL;
		}
	}

}

#define SQR(X) ((X)*(X))
double compute_sigma(float data[], int size, double avg)
{
	double sum;
	int count;
	int i;

	sum = 0.0;
	count = 0;
	for (i=0; i<size; i++) {
		if (!IS_BLANK_PIXEL(data[i])) {
			sum += SQR(data[i] - avg);
			count++;	
		}
	}

	return sqrt(sum / (count-1) );
}


int main(int argc, char * argv[])
{

	header_param_list Ihpar, Qhpar, Uhpar, Vhpar, Iouthpar;
	header_param_list Iavghpar;
	FILE * Ifile, *Qfile, *Ufile, *Vfile, *datfile;
	char * Ifilename, *Qfilename, *Ufilename, *Vfilename;
	char * Iavgfilename;
	int num_planes, num_points;
	int i, j;
	char *datfilename = "empirical.dat";
	char *Ioutfilename = "empirical_Inoisemap.fits";
	float *Iplane, *Qplane, *Uplane, *Vplane;
	float *Iavgdata;
	int Imaxpos;
	float Imax;
	char * ranges;
	byte planes[MAX_NUM_PLANES];


	if (argc < 6) {
		printf("usage: %s <Iavg.fits> <I.fits> <Q.fits> <U.fits> <V.fits> [range1,range2,...]\n", argv[0]);
		return EXIT_FAILURE;
	}
	
	Iavgfilename = argv[1];
	Ifilename = argv[2];
	Qfilename = argv[3];
	Ufilename = argv[4];
	Vfilename = argv[5];
	ranges = (argc == 7) ? argv[6] : NULL;
	
	select_ranges(ranges, planes);
printf("ranges: %s\n", ranges);

	Ifile = fopen(Ifilename, "r");
	Qfile = fopen(Qfilename, "r");
	Ufile = fopen(Ufilename, "r");
	Vfile = fopen(Vfilename, "r");
	if (Ifile == NULL || Qfile == NULL || Ufile == NULL || Vfile == NULL) {
		printf("ERROR: unable to open input file for reading\n");
		return EXIT_FAILURE;
	}

	datfile = fopen(datfilename, "w");
	if (datfile == NULL) {
		printf("ERROR: unable to open output file for writing\n");
		return EXIT_FAILURE;
	}


	//average maps

	readfits_map (Iavgfilename, &Iavgdata, &Iavghpar);
	num_points = Iavghpar.naxis[0] * Iavghpar.naxis[1];

	//find max pixel in I
	Imaxpos = -1;
	Imax = 0;
	for (j=0; j<num_points; j++) {
		if (Imax < Iavgdata[j]) {
			Imax = Iavgdata[j];
			Imaxpos = j;
		}
	}

	blank_pixels(Iavgdata, num_points);
	strcat(Iavghpar.object, " Noisemap");
	writefits_map (Ioutfilename, Iavgdata, &Iavghpar);

	//cubes
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

	Iavgdata = malloc(sizeof(float) * num_points);
	Iplane = malloc(sizeof(float) * num_points);
	Qplane = malloc(sizeof(float) * num_points);
	Uplane = malloc(sizeof(float) * num_points);
	Vplane = malloc(sizeof(float) * num_points);


	for (i=0; i<num_planes; i++) 
	{
		float Qmax, Umax, Vmax;
		double mI, mQ, mU, mV;
		double Inot, Qnot, Unot, Vnot;
		double deltaInot, deltaQnot, deltaUnot, deltaVnot;
		int count = 1;

		readfits_plane (Ifile, Iplane, &Ihpar);
		readfits_plane (Qfile, Qplane, &Qhpar);
		readfits_plane (Ufile, Uplane, &Uhpar);
		readfits_plane (Vfile, Vplane, &Vhpar);

		//special case for handling exclusions
		if (i < num_planes-1 && planes[i+1] == 0) 
		{
			float *IplaneX, *QplaneX, *UplaneX, *VplaneX;
			IplaneX = malloc(sizeof(float) * num_points);
			QplaneX = malloc(sizeof(float) * num_points);
			UplaneX = malloc(sizeof(float) * num_points);
			VplaneX = malloc(sizeof(float) * num_points);
			do {
				readfits_plane (Ifile, IplaneX, &Ihpar);
				readfits_plane (Qfile, QplaneX, &Qhpar);
				readfits_plane (Ufile, UplaneX, &Uhpar);
				readfits_plane (Vfile, VplaneX, &Vhpar);
				count++;
				i++;
			} while (i<num_planes && planes[i] == 0);

			//simple average of the datapoints
			for (j=0; j<num_points; j++) {
				Iplane[j] = (Iplane[j] + IplaneX[j]) / 2.0;
				Qplane[j] = (Qplane[j] + QplaneX[j]) / 2.0;
				Uplane[j] = (Uplane[j] + UplaneX[j]) / 2.0;
				Vplane[j] = (Vplane[j] + VplaneX[j]) / 2.0;
			}
			free(IplaneX);
			free(QplaneX);
			free(UplaneX);
			free(VplaneX);
		}

		//max values from this plane/
		Imax = Iplane[Imaxpos];
		Qmax = Qplane[Imaxpos];
		Umax = Uplane[Imaxpos];
		Vmax = Vplane[Imaxpos];

		//Iavgdata is used as a mask to blank out non-background pixels
		blank_corresponding_pixels (Iavgdata, Iplane, num_points);
		blank_corresponding_pixels (Iavgdata, Qplane, num_points);
		blank_corresponding_pixels (Iavgdata, Uplane, num_points);
		blank_corresponding_pixels (Iavgdata, Vplane, num_points);

		//averages of background pixels
		Inot = compute_average (Iplane, num_points);
		Qnot = compute_average (Qplane, num_points);
		Unot = compute_average (Uplane, num_points);
		Vnot = compute_average (Vplane, num_points);

		//correction factors
		mQ = (Qmax - Qnot) / (Imax - Inot);
		mU = (Umax - Unot) / (Imax - Inot);
		mV = (Vmax - Vnot) / (Imax - Inot);

		//used for error analysis later
		deltaInot = compute_sigma (Iplane, num_points, Inot); 
		deltaQnot = compute_sigma (Qplane, num_points, Qnot);
		deltaUnot = compute_sigma (Uplane, num_points, Unot);
		deltaVnot = compute_sigma (Vplane, num_points, Vnot);

		//write out the results to the datafile using the following format:
		//plane mQ mU mV Inot Qnot Unot Vnot Isigma Qsigma Usigma Vsigma Imax Qmax Umax Vmax 
		for (j=i-count; j<i; j++) {
			fprintf(datfile, "%3i %.6f %.6f %.6f %.6f %.6f %.6f %.6f "
					"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
					j+1, mQ, mU, mV, 
					Inot, Qnot, Unot, Vnot, 
					deltaInot, deltaQnot, deltaUnot, deltaVnot, 
					Imax, Qmax, Umax, Vmax);
		}

	}

	free(Iplane);
	free(Qplane);
	free(Uplane);
	free(Vplane);
	free(Iavgdata);

	fclose(Ifile);
	fclose(Qfile);
	fclose(Ufile);
	fclose(Vfile);
	fclose(datfile);
}



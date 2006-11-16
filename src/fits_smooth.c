/* 
 * smooth.c
 *
 * Perform a smoothing operation through the cube
 *
 * Uses fitsio written by Steve Gibson.
 *
 * Jeff Dever - 2005 University of Calgary Radio Astronomy
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsio.h"

#define MAX_SMOOTH_WIDTH (4096)

static float H[MAX_SMOOTH_WIDTH];

enum SmoothType {HANNING, UNIFORM};

#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )


/* 
 * Hanning multipliers are stored in the first n positions of 
 * the static array H.  The sum of these must be 1.
 * n - the Hanning number and must be odd
 * 
 * This is quite a tricky function.  Relies on the fact that the 
 * Hanning function is a triangle, and that the numbers have a 
 * Triangular property.  I could probablly write a paper about this code  ...
 */
static void compute_hanning_function(int n)
{
	int k, i;
	float denom;

	k = n/2;
	denom = 2*((k*k-k)/2.0 + k) + (n+1)/2.0;
	for (i=1; i<=(n+1)/2; i++) {
		H[i-1] = H[n-i] = i/denom;
	}
}


static void print_usage(char * argv[])
{
	printf("\n");
	printf("Usage: %s <in.fits> <out.fits> [smooth_type] [smooth_width]\n", argv[0]);
	printf("smooth_type - one of {hanning, uniform}\n");
	printf("smooth_width - odd natural number less than %i\n", MAX_SMOOTH_WIDTH);
	printf("\n");
	printf("eg: %s cube.fits cube_smooth.fits hanning 3\n", argv[0]);
	printf("\n");
}


int main(int argc, char *argv[]) 
{

	int i, j, k, n;
	char *infilename, *outfilename;
	header_param_list inhpar, outhpar;
	float *indata, *outdata;
	int plane_size;
	int num_planes;
	const char * smooth_name;
	enum SmoothType smooth_type;
	int N;


	if (argc < 3) {
		print_usage(argv);
		return EXIT_FAILURE;
	}
	infilename = argv[1];
	outfilename = argv[2];
	smooth_name = (argc >= 4) ? argv[3] : "hanning";
	N = (argc >= 5) ? atoi(argv[4]) : 3;


	if (N/2 == (N+1)/2) {
		printf("ERROR: An even number (%i) for the smooth window is not valid\n", N);
		return EXIT_FAILURE;
	}

	if (strcmp(smooth_name, "hanning") == 0) 
	{
		smooth_type = HANNING;
		compute_hanning_function(N);
		printf("Using Hanning-%i smoothing\n", N);
	}
	else if (strcmp(smooth_name, "uniform") == 0) 
	{
		smooth_type = UNIFORM;
		printf("Using Uniform-%i smoothing\n", N);
	}
	else 
	{
		printf("ERROR: Unknown smoothing type: '%s'\n", smooth_name);
		print_usage(argv);
		return EXIT_FAILURE;
	}

	readfits_cube (infilename, &indata, &inhpar, 0);

	plane_size = inhpar.naxis[0] * inhpar.naxis[1];
	num_planes = inhpar.naxis[2];
	printf("num_planes: %i\n", num_planes);
	printf("plane_size: %i\n", plane_size);

	outdata = malloc (sizeof(float) * plane_size * num_planes);

	//iterate over the cube and smooth the planes
	for (i=0; i<plane_size; i++) 
	{
		for (j=0; j<num_planes; j++) 
		{
			double sum = 0.0;

			switch (smooth_type) {
				case HANNING:
					{
						for (k=j-N/2, n=0; n<N; k++, n++) 
						{
							int p;
							float datum;

							p = k;
							if (p < 0) p = 0;  
							if (p >= num_planes) p = num_planes-1;
							datum = indata[p*plane_size+i];

							if (IS_BLANK_PIXEL(datum)) {
								sum = indata[j*plane_size+i];
								break;
							} else {
								sum += datum * H[n];
							}
						}
						break;
					}
				case UNIFORM:
					{
						int count = 0;
						for (k=j-N/2, n=0; k<j+N/2 && n<=N; k++, n++)
						{
							float datum;

							if (k < 0) k = 0;  
							if (k >= num_planes) k = num_planes-1;
							datum = indata[k*plane_size+i];

							if (!IS_BLANK_PIXEL(datum)) {
								sum += datum;
							}
							count++;
						}
						sum = sum / count;
						break;
					}
			}
			outdata[j*plane_size+i] = (float) sum;
		}
	}

	//write the smoothed cube
	outhpar = inhpar;
	writefits_cube(outfilename, outdata, &outhpar, 0);

	//cleanup
	free(indata);
	free(outdata);

	return EXIT_SUCCESS;
}


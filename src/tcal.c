#include <stdio.h>
#include <stdlib.h>
#include "programs/chardefs.h"
#include "programs/mathdefs.h"
#include "programs/misc.c"
#include "programs/misc_math.c"
#include "programs/coord_utils.c"
#include "programs/fitsio.c"

#include "common.h"
#include "jsd/jsd_fit.h"


#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )


#define MAX_SMOOTH_WIDTH (4096)
static float H[MAX_SMOOTH_WIDTH];
int hs;



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
	hs = n;
}

void apply_hanning_function(double A[], double B[], int size)
{
	int n, i;
	for (i=hs/2; i<size-hs/2; i++) {
		double sum = 0;
		for (n=0; n<hs; n++) {
			sum += A[i-hs/2+n] * H[n];
		}
		B[i] = sum;
	}
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

//smooth the cal using polynomial fit
void polynomial_smooth(double Tcal[], double Tcal_smooth[], int size, int order)
{
	int i;
	double C[16];
	double chisq;
	double X[MAX_CHANNELS];
	const float nsigma = 3.0;

	for (i=0; i<size; i++) {
		X[i] = i;
	}
	jsd_poly_fit(X, Tcal, size, C, nsigma, order, &chisq);
	for (i=0; i<size; i++) {
		Tcal_smooth[i] = jsd_poly_eval(i, C, order);
	}
	jsd_print_poly(stdout, C, order);
}


int main(int argc, char * argv[])
{
	header_param_list XXhpar, XXavghpar, YYhpar, YYavghpar;
	char * XXavgfilename, *XXfilename, *YYavgfilename, *YYfilename;
	int num_planes, num_points;
	char *datfilename = "Tcal.dat";
	char *XXoutfilename = "XXnoisemap.fits";
	char *YYoutfilename = "YYnoisemap.fits";
	FILE * datfile;
	float *XXavgdata, *YYavgdata;
	float *XXdata, *YYdata;
	int XXmaxpos, YYmaxpos;
	float XXmax, YYmax;
	double xx_off[MAX_CHANNELS], yy_off[MAX_CHANNELS]; 
	double xx_on[MAX_CHANNELS], yy_on[MAX_CHANNELS]; 
	double xx[MAX_CHANNELS], yy[MAX_CHANNELS]; 
	double Tcalx[MAX_CHANNELS], Tcaly[MAX_CHANNELS];
	double Tcalx_smooth[MAX_CHANNELS], Tcaly_smooth[MAX_CHANNELS];
	double Tcalymean;
	int i;


	if (argc != 5) {
		printf("usage: %s <XXavg.fits> <XX.fits> <YYavg.fits> <YY.fits>\n", argv[0]);
		return EXIT_FAILURE;
	}
	
	XXavgfilename = argv[1];
	XXfilename = argv[2];
	YYavgfilename = argv[3];
	YYfilename = argv[4];

	//read average maps

	readfits_map (XXavgfilename, &XXavgdata, &XXavghpar);
	readfits_map (YYavgfilename, &YYavgdata, &YYavghpar);
	if (memcmp(XXavghpar.naxis, YYavghpar.naxis, 2) != 0) {
		printf("ERROR: naxis size mismatch between avg fits headers\n");
		return EXIT_FAILURE;
	}
	num_points = XXavghpar.naxis[0] * XXavghpar.naxis[1];

	//find max pixels
	XXmaxpos = -1;
	XXmax = 0;
	for (i=0; i<num_points; i++) {
		if (XXmax < XXavgdata[i]) {
			XXmax = XXavgdata[i];
			XXmaxpos = i;
		}
	}
	YYmaxpos = -1;
	YYmax = 0;
	for (i=0; i<num_points; i++) {
		if (YYmax < YYavgdata[i]) {
			YYmax = YYavgdata[i];
			YYmaxpos = i;
		}
	}

	if (YYmaxpos != XXmaxpos) {
		printf("WARN: XX and YY do not have the same peak pixel\n");
	}

	//blank and write noisemaps
	blank_pixels(XXavgdata, num_points);
	blank_pixels(YYavgdata, num_points);
	strcat(XXavghpar.object, " Noisemap");
	strcat(YYavghpar.object, " Noisemap");
	writefits_map (XXoutfilename, XXavgdata, &XXavghpar);
	writefits_map (YYoutfilename, YYavgdata, &YYavghpar);


	//read cubes

	readfits_cube (XXfilename, &XXdata, &XXhpar);
	readfits_cube (YYfilename, &YYdata, &YYhpar);
	
	if (memcmp(XXhpar.naxis, YYhpar.naxis, 3) != 0) {
		printf("ERROR: naxis size mismatch between fits headers\n");
		return EXIT_FAILURE;
	}

	if (XXhpar.naxis[0] * XXhpar.naxis[1] != num_points) {
		printf("ERROR: cube sizes do not match averages\n");
	}
	num_planes = XXhpar.naxis[2];
	printf("num_points: %i\n", num_points);
	printf("num_planes: %i\n", num_planes);

	//extract the on spectra
	for (i=0; i<num_planes; i++) {
		xx_on[i] = XXdata[i*num_points+XXmaxpos];
		yy_on[i] = YYdata[i*num_points+YYmaxpos];
	}

	//determine the off spectra
	for (i=0; i<num_planes; i++) {
		float *XXplane = &XXdata[i*num_points];
		float *YYplane = &YYdata[i*num_points];
		blank_corresponding_pixels (XXavgdata, XXplane, num_points);
		blank_corresponding_pixels (YYavgdata, YYplane, num_points);
		xx_off[i] = compute_average(XXplane, num_points);
		yy_off[i] = compute_average(YYplane, num_points);
	}

	//derive xx and yy spectra
	for (i=0; i<num_planes; i++) {
		xx[i] = xx_on[i] - xx_off[i];
		yy[i] = yy_on[i] - yy_off[i];
	}

	//compute initial tcal
	for (i=0; i<num_planes; i++) {
		Tcalx[i] = 1.0/xx[i];
		Tcaly[i] = 1.0/yy[i];
	}

	//compute average tcaly
	Tcalymean = 0.0;
	for (i=0; i<num_planes; i++) {
		Tcalymean += Tcaly[i];
	}
	Tcalymean /= num_planes;

	//compute final values of tcal
	for (i=0; i<num_planes; i++) {
		Tcalx[i] /= Tcalymean;
		Tcaly[i] /= Tcalymean;
	}
	
/*
	//smooth the tcal using Hanning
	{
		printf("performing Hanning smoothing\n");
		compute_hanning_function(11);
		for (i=0; i<num_planes; i++) {
			Tcalx_smooth[i] = Tcalx[i];
			Tcaly_smooth[i] = Tcaly[i];
		}
		apply_hanning_function(Tcalx, Tcalx_smooth, num_planes);
		apply_hanning_function(Tcaly, Tcaly_smooth, num_planes);
	}
*/

	//smooth using a polynomial
	printf("performing polynomial fit\n");
	polynomial_smooth(Tcalx, Tcalx_smooth, num_planes, 8);
	polynomial_smooth(Tcaly, Tcaly_smooth, num_planes, 8);


	//write out the results to the datafile
	datfile = fopen(datfilename, "w");
	if (datfile == NULL) {
		printf("ERROR: unable to open output file for writing\n");
		return EXIT_FAILURE;
	}

	fprintf(datfile, "#Tcalx_smooth Tcaly_smooth Tcalx Tcaly\n");
	for (i=0; i<num_planes; i++) {
		fprintf(datfile, "%.6f %.6f %.6f %.6f\n",
			Tcalx_smooth[i], Tcaly_smooth[i], Tcalx[i], Tcaly[i]);
	}

	fclose(datfile);

	free(XXavgdata);
	free(XXdata);
	free(YYavgdata);
	free(YYdata);

}



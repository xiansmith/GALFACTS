/*
Quick and dirty processing code to investigate cal variations.
Creates two seperate averages:
A - the low channels (below the median)
B - the high channels (above the median)

Mean and sigma are written to stdout for A and B, and then once more for A-B.
Data files are written for each of A and B, with columns at successively higher
smoothing levels.

Jeff Dever
2006/10/17

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsio.h"


//the low and high channel range to avoid the band rolloff channels.
#define LOWCHAN  25
#define HIGHCHAN  230


//average ignores blank pixels
static double compute_mean(double data[], int size)
{
	double sum;
	int count;
	int i;

	sum = 0.0;
	count = 0;
	for (i=0; i<size; i++) {
		if (isnormal(data[i])) {
			sum += data[i];
			count++;
		}
	}

	return sum/count;
}


#define SQR(X) ((X)*(X))
static double compute_sigma(double data[], int size, double avg)
{
	double sum;
	int count;
	int i;

	sum = 0.0;
	count = 0;
	for (i=0; i<size; i++) {
		if (isnormal(data[i])) {
			sum += SQR(data[i] - avg);
			count++;	
		}
	}

	return sqrt(sum / (count-1) );
}

/*
Performs uniform (boxcar) smoothing on a linear block of data.
@param A input data to be smoothed
@param B output data that is smoothed
@param N size of the arrays
@param width width of the smoothing window
*/
static void uniform_smooth(const double A[], double B[], int N, int width)
{
	int j, k;

	memset(B, 0, sizeof(double) * N);

	for (j=0; j<N; j++) 
	{
		int count = 0;
		for (k=j-width/2; k<j+width/2; k++)
		{
			if (k < 0) continue;
			if (k > N) break; 
			if (isnormal(A[k])) {
				B[j] += A[k];
				count++;
			}
		}
		B[j] /= count;
	}
}

/*
 * Produces a plain text file with the given filename that contains a row for each
 * time point and a column for various degrees of smoothing of the data.
 * The first column is unsmoothed.
 * This could have been much more adaptive, but whatever.
 */
static void produce_smoothies(const double data[], int N, const char *filename)
{
	int i;
	FILE * file;
	double *A5; 
	double *A50; 
	double *A300; 
	double *A3000;

	A5 = (double *) calloc(N, sizeof(double)); 
	A50 = (double *) calloc(N, sizeof(double)); 
	A300 = (double *) calloc(N, sizeof(double)); 
	A3000 = (double *) calloc(N, sizeof(double)); 

	uniform_smooth(data, A5, N, 5); //1sec
	uniform_smooth(data, A50, N, 50); //10sec
	uniform_smooth(data, A300, N, 300); //1min
	uniform_smooth(data, A3000, N, 3000); //10min

	file = fopen(filename, "w");
	if (file == NULL) {
		printf("ERROR: unable to open '%s'\n", filename);
		return;
	}

	fprintf(file, "#raw 1sec 10sec 1min 10min\n");

	for (i=0; i<N; i++) {
		fprintf(file, "%f5.8 %f5.8 %f5.8 %f5.8 %f5.8\n", data[i], A5[i], A50[i], A300[i], A3000[i]);
	}

	fclose(file);

	free(A5);
	free(A50);
	free(A300);
	free(A3000);
}


/*
 * Sets any datapoints that are nsigma away from the mean to be NAN.
 * Does this iteratively as required.
 */
static void remove_spikes(double data[], int N, float nsigma)
{
	double sigma,mean;
	int i;
	int repeat;

	do {
		repeat = 0;
		mean = compute_mean(data, N);
		sigma = compute_sigma(data, N, mean);
		for (i=0; i<N; i++) {
			if (isnormal(data[i])) {
				if (fabs(mean-data[i]) > (nsigma*sigma)) {
					data[i] = NAN;
					repeat = 1;
				}
			}
		}
	} while (repeat);
}



int main(int argc, char * argv[])
{

	header_param_list hpar;
	char * filename;
	int i, j;
	int N, C;
	float *data;
	double *Adata;
	double *Bdata;
	int lowchan = 25;
	int highchan = 230;

	if (argc < 2) {
		printf("usage: %s <spec.fits>]\n", argv[0]);
		return EXIT_FAILURE;
	}
	filename = argv[1];

	readfits_map (filename, &data, &hpar);

	//dimensions of the data
	C = hpar.naxis[0];
	N = hpar.naxis[1];

	//Average the data into the two sub bands A,B
	Adata = (double *) calloc(N, sizeof(double));
	Bdata = (double *) calloc(N, sizeof(double));
	for (i=0; i<N; i++) 
	{
		for (j=lowchan; j<C/2; j++) {
			Adata[i] += data[i*C+j];
		}
		Adata[i] /= C/2 - lowchan;
		for (j=C/2; j<highchan; j++) {
			Bdata[i] += data[i*C+j];
		}
		Bdata[i] /= highchan - C/2;
	}

	//remove the hical and doppler correction spikes a sigma multiplier
	remove_spikes(Adata, N, 5.0);
	remove_spikes(Bdata, N, 5.0);
	
	//now we have a sub band of data, do the smoothing on it
	produce_smoothies(Adata, N, "A_smooth.dat");
	produce_smoothies(Bdata, N, "B_smooth.dat");

	//do the subtraction, compute mean and sigma
	{
		double mean, sigma;
		double * Cdata = (double *) calloc(N, sizeof(double));

		for (i=0; i<N; i++) {
			Cdata[i] = Adata[i] - Bdata[i];
		}

		mean = compute_mean(Cdata, N);
		sigma = compute_sigma(Cdata, N, mean);
		printf("mean=%f sigma=%f\n", mean, sigma);
	}

	//done
	return EXIT_SUCCESS;
}



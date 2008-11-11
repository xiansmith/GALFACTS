#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "smoothing.h"



/* 
 * Uses Box-Muller transform to produce normally distributed
 * random numbers with the given sigma.
 *
 * You may wish to call srand48 to initialize the rand48 before
 * calling this function.
 *
 */
double rng_gaussian (const double sigma)
{
	double x, y, r2;

	do
	{
		x = 2 * drand48() - 1;
		y = 2 * drand48() - 1;
		r2 = x * x + y * y;
	}
	while (r2 > 1.0 || r2 == 0);

	return sigma * y * sqrt (-2.0 * log (r2) / r2);
}


void create_noisy_data(double A[], int array_size, double freq, double noise_sigma)
{
	int i;
	for (i=0; i<array_size; i++) {
		A[i] = sin(2*M_PI*i/(double)array_size * freq) + rng_gaussian(noise_sigma)/20;
	}
}


void test_smoothing(int array_size, int smooth_size)
{
	double *A, *B, *C, *H;
	int i;
	double sum;
	double noise_sigma = 5.0;
	double freq = 3;

	A = calloc(array_size, sizeof(double));
	B = calloc(array_size, sizeof(double));
	C = calloc(array_size, sizeof(double));
	H = calloc(smooth_size, sizeof(double));

	printf("#array size: %i\n", array_size);
	printf("#smooth size: %i\n", smooth_size);

	//create some noisy data
	create_noisy_data(A, array_size, freq, noise_sigma);

	//produce some bad data
	for (i=0; i<25; i++) {
		A[i] = NAN;
	}
	for (i=300; i<320; i++) {
		A[i] = NAN;
	}
	for (i=500; i<530; i++) {
		A[i] = NAN;
	}
	for (i=800; i<840; i++) {
		A[i] = NAN;
	}
	for (i=1200; i<1280; i++) {
		A[i] = NAN;
	}
	for (i=1600; i<1700; i++) {
		A[i] = NAN;
	}

	for (i=array_size-50; i<array_size; i++) {
		A[i] = NAN;
	}

	compute_boxcar_coefficients(H, smooth_size);
	apply_smoothing_function(A, B, array_size, H, smooth_size);
	sum = 0;
	for (i=0; i<smooth_size; i++) {
		sum += H[i];
	}
	printf("#uniform sum: %g\n", sum);


	compute_hanning_coefficients(H, smooth_size);
	apply_smoothing_function(A, C, array_size, H, smooth_size);
	sum = 0;
	for (i=0; i<smooth_size; i++) {
		sum += H[i];
	}
	printf("#hanning sum: %g\n", sum);


	for (i=0; i<array_size; i++) {
		printf("%i %g %g %g\n", i, A[i], B[i], C[i]);
	}

	free(H);
	free(C);
	free(B);
	free(A);
}


int main(int argc, char *argv[])
{
	test_smoothing(2000, 99);
	return EXIT_SUCCESS;
}


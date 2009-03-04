#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "smoothing.h"

/*
 * Smoothing functions must have the property that the sum of all the elements is 1
 *
 *
 */


static int interpolate_iteration(double A[], int size, int initial);


/* 
 * Hanning multipliers are stored in the first n positions of 
 * the array H.  The sum of these must be 1.
 * n - the Hanning number and must be odd
 * H - must be allocated to be of size n
 * 
 * This is quite a tricky function.  Relies on the fact that the 
 * Hanning function is a triangle, and that the numbers have a 
 * Triangular property.  I could probablly write a paper about this code  ...
 */
void compute_hanning_coefficients(double H[], int n)
{
	int k, i;
	double denom;

	k = n/2;
	denom = 2*((k*k-k)/2.0 + k) + (n+1)/2.0;
	for (i=1; i<=(n+1)/2; i++) {
		H[i-1] = H[n-i] = i/denom;
	}
}

void compute_boxcar_coefficients(double H[], int n)
{
	int i;
	double val = 1.0/n;
	for (i=0; i<n; i++) {
		H[i] = val;
	}
}

void interpolate_missing_data(double A[], int size)
{
	int pos;
	
	pos = 0;
	do {
		pos = interpolate_iteration(A, size, pos);
	} while (pos > 0);
}



/* returns the position of the largest element that is guranteed to 
 * be all finite.
 * On the first call the initial value should be zero, and subsequent calls should be
 * the return value intil -1 is returned.
 */
static int interpolate_iteration(double A[], int size, int initial)
{
	int i;
	int start, end;
	double val;
	
	start = end = -1;
	i = initial;
	
	//set start as the first non-finite element of the range
	//will be -1 if NAN extends off the low end
	for ( ; i<size; i++) {
		if (! finite(A[i])) {
			start = i-1;
			break;
		}
	}

	//set end as the last non-finite element of the range
	//will be -1 if the NAN extends off the high end
	for ( ; i<size; i++) {
		if (finite(A[i])) {
			end = i;
			break;
		}
	}

	//no need to continue since data is either all finite or all infinite
	if (start < 0 && end < 0) {
		return -1;
	}

	//do the interpolation in the case where NAN extends off to the left
	if (start < 0) {
		val = A[end];
		for (i=0; i<end; i++) {
			A[i] = val;
		}
		return end;  //start to check at the end position next time
	}


	//do the interpolation in the case where NAN extends off the the right
	if (end < 0) {
		val = A[start];
		for (i=start; i<size; i++) {
			A[i] = val;
		}
		return -1;  //there must be no more NAN
	}

	//since none of end or start are NAN, do the interpolation over a range
	val = (A[start] + A[end]) / 2.0;
	for (i=start+1; i<end; i++) {
		A[i] = val;
	}
	return end;   //start to check at the end position next time

}




/*
 * A - input data to smooth (must exist)
 * B - output array of smoothed data (must be memory allocated)
 * size - the size of array A and B
 * H - the array of smoothing co-efficients
 * n - the size of H
 *
 * 
 *
 */
void apply_smoothing_function(const double A[], double B[], int size, const double H[], int n)
{
	int i, j, indx;
	
	for (i=0; i<size; i++) 
	{
		double sum = 0.0;
		double scale = 0.0; //used to normalize when 
		for (j=0; j<n; j++) {
			indx = i + j-n/2;
			if (indx < 0) continue;
			if (indx >= size) break;
			if (!finite(A[indx])) continue;
			sum += A[indx] * H[j];
			scale += H[j];
		}
		B[i] = sum / scale;
	}
} 


void hanning_smooth_data(double A[], int size, int n)
{
	double *B, *H;
	int i;

	printf("Requesting malloc for %ld bytes of memory\n.",sizeof(double)*size);
	B = calloc(size, sizeof(double));
	if (B == NULL) {
		printf("ERROR: malloc failed in hanning_smooth_data() !\n");
	}
	printf("Requesting malloc for %ld bytes of memory\n.",sizeof(double)*n);
	H = calloc(n, sizeof(double));
	if (H == NULL) {
		printf("ERROR: malloc failed in hanning_smooth_data() !\n");
	}
		
	compute_hanning_coefficients(H, n);
	apply_smoothing_function(A, B, size, H, n);

	for (i=0; i<size; i++) {
		A[i] = B[i];
	}

	free(B);
	free(H);
}



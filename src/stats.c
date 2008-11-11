#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "stats.h"

/*
 * Computes the statistical mean of the given data array
 */
double compute_mean(double data[], int size)
{
	double sum;
	int count;
	int i;

	sum = 0.0;
	count = 0;
	for (i=0; i<size; i++) {
		if (isfinite(data[i])) {
			sum += data[i];
			count++;
		}
	}

	return (count>0) ? sum/count : NAN;
}

#define SQR(X) ((X)*(X))
double compute_sigma(double data[], int size, double mean)
{
	double sum;
	int count;
	int i;

	sum = 0.0;
	count = 0;
	for (i=0; i<size; i++) {
		if (isfinite(data[i])) {
			sum += SQR(data[i] - mean);
			count++;	
		}
	}
	return (double) sqrt(sum / (count-1) );
}

/* 
 * Assumes the data is gaussian, rejects outliers that are nsigma away from the mean.
 * Sets outlier values to be NAN in the input data array.
 * returns the number of rejected values.
 */
int reject_outliers(double data[], int size, float nsigma)
{
	double mean, sigma;
	int i, reject_count;

	mean = compute_mean(data, size);
	sigma = compute_sigma(data, size, mean);
	reject_count = 0;
	for (i=0; i<size; i++) {
		if (isnan(data[i])) continue;
		if (fabs(mean - data[i]) > nsigma*sigma) {
			data[i] = NAN;
			reject_count++;
		}
	}
	return reject_count;
}

double compute_clean_mean(double data[], int size, float nsigma)
{
	do {} while (reject_outliers(data, size, nsigma) > 0);
	return compute_mean(data, size);
}


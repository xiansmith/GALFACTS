//#include "decdependence.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "jsd/jsd_fit.h"
#include "jsd/jsd_futil.h"
#include "stats.h"

/*
 * use a single channel worth of data
 * iterate over each beam
 * iterate over every day
 * put points into DEC bins
 * reject outliers in the bins
 * fit a final curve to the dec bins
 * then subtract the amount of that curve from the data
 * 
 */

/*
static void compute_dirty_bins (double **binarrayX, int countX[], double binX[], int num_bins)
{
	int b;
	for (b=0; b<num_bins; b++)
	{
		binX[b] = compute_mean(binarrayX[b], 0,countX[b]);
	}
}
*/

//static void compute_clean_bins (double **binarrayX, int *countX, double *binX, int num_bins, float nsigma, int *outliercounts)
static void compute_clean_bins (double **binarrayX, int *countX, double *binX, double *sigmaX,int num_bins)
{
	int b, i;
	double mean, sigma;
	int outlier;
	double *binarray;
	int count;
	float nsigma = 2.5;
	//memset(outliercounts, 0, num_bins*sizeof(int));
	for (b=0; b<num_bins; b++)
	{
		//repeat until no more outliers
		do {
			outlier = 0;
			binarray = binarrayX[b];
			count = countX[b];
			mean = compute_mean(binarray, 0,count);
			sigma = compute_sigma(binarray, count, mean);

			//do the cleaning by setting outliers to NAN
			for (i=0; i<count; i++) {
				if (fabs(mean-binarray[i]) > nsigma*sigma) {
					//binarray[i] = NAN;
					binarray[i] = mean;
					outlier = 1;
	//				outliercounts[b]++;
				}
			}
		} while (outlier);
		binX[b] = mean;
		sigmaX[b] = sigma;
	}
}


static void create_dec_bins(double * I,double * V, int numvalues,double **binarrayV, int *counts, float decmin, float decgrain, int num_bins)
{
	int d, r;

	memset(counts, 0, num_bins*sizeof(int));

	for (d=0; d<numvalues; d++) 
	{
			int bin = (int) floor((I[d] - decmin) / decgrain);

			if (bin<0 || bin>=num_bins) continue;
//			if (isfinite(I)) {
				double *binarray;
				binarray = binarrayV[bin];
				binarray[counts[bin]] = V[d];
				counts[bin]++;
//			}
	}
}
#define MAX_BIN_SIZE 1000
static void beam_dec_dependence(double * I, double * V,  int numvalues, float decmin, float decmax, float decgrain)
{
	int i;
	int num_bins;
	double **binarrayV;
	double *cleanbinV;
	double *sigmabinV;
	int * counts, *outliercounts;
	FILE * decfile;
	char filename[32+1];
	float nsigma = 5.0; //for the per bin data point exclusion
	num_bins = (int) ceil((decmax-decmin)/decgrain);
	decgrain = (decmax-decmin)/num_bins; //actual grain (bin) size

	binarrayV = (double**) malloc(num_bins * sizeof(double*));
	for (i=0; i<num_bins; i++) {
		binarrayV[i] = (double*) malloc(MAX_BIN_SIZE * sizeof(double));
	}
	counts = malloc(num_bins * sizeof(int));
//	outliercounts = malloc(num_bins * sizeof(double));
	cleanbinV = malloc(num_bins * sizeof(double));
	sigmabinV = malloc(num_bins * sizeof(double));

	//bin up the values in declination
	create_dec_bins(I, V, numvalues,binarrayV, counts, decmin, decgrain, num_bins);
	//compute_clean_bins(binarrayV, counts, cleanbinV, num_bins, nsigma, outliercounts);
	compute_clean_bins(binarrayV, counts, cleanbinV, sigmabinV,num_bins);

	//print bin data
//	snprintf(filename, 32, "decbins_beam%i_chan%i.dat", beam, chan); 
	decfile = fopen("posterbins.txt", "w");
	fprintf(decfile, "#bin I counts cleanV sigmaV V%\n");
	for (i=0; i<num_bins; i++) {
		//fprintf(decfile, "%i %f %i %g%g\n", i, i*decgrain+decmin, counts[i], cleanbinI[i], cleanbinV[i]);
		printf("%i %f %i %g\n", i, i*decgrain+decmin, counts[i], cleanbinV[i]);
		fprintf(decfile,"%i %f %i %g %g %g\n", i, i*decgrain+decmin+decgrain/2, counts[i], cleanbinV[i],sigmabinV[i],50*cleanbinV[i]/(i*decgrain+decmin+decgrain/2));
	}
	fclose(decfile);


	//do it again just so we can print out the bin data and see the residuals
//	create_dec_bins(I, V,  binarrayV, counts,  decmin, decgrain, num_bins);
//	compute_clean_bins(binarrayV, counts, cleanbinV, num_bins, nsigma, outliercounts);

//	snprintf(filename, 32, "decbinsnew_beam%i_chan%i.dat", beam, chan); 
//	decfile = fopen(filename, "w");
//	fprintf(decfile, "#bin DEC counts cleanI cleanQ cleanU cleanV\n");
//	for (i=0; i<num_bins; i++) {
//		fprintf(decfile, "%i %f %i %g %g\n", i, i*decgrain+decmin, counts[i], cleanbinI[i], cleanbinV[i]);
//	}
//	fclose(decfile);


	//cleanup
	for (i=0; i<num_bins; i++) {
//		free(binarrayI[i]);
		free(binarrayV[i]);
	}
//	free(binarrayI);
	free(binarrayV);
//	free(cleanbinI);
	free(cleanbinV);
	free(sigmabinV);
	free(counts);
}

void main()
{
	FILE *posterfile;
	//posterfile =fopen("totaldata.txt","r");
	posterfile =fopen("VvsIdata.txt","r");
	int count = jsd_line_count(posterfile);
	double *I,*V;
	I = (double*)malloc(sizeof(double)*count);
	V = (double*)malloc(sizeof(double)*count);
	double temp;
	int i;
	for(i=0;i<count;i++)
	{
		fscanf(posterfile,"%lf %lf %lf %lf %lf %lf %lf %lf",&temp,&temp,&temp,&I[i],&temp,&V[i],&temp,&temp);
		//printf("%f %f\n",I[i],V[i]);
	}
	beam_dec_dependence(I, V, count, 0.0, 10.0, 0.5);
}

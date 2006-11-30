#include "decdependence.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "jsd/jsd_fit.h"

/*
 * use a single channel worth of data
 * iterate over each beam
 * iterate over every day
 * fit a curve to each scan, reject outliers (sources)
 * take remaining points and put them into DEC bins
 * fit a final curve to the dec bins
 * then subtract the amount of that curve from the data
 * 
 */

#define BIN_ORDER 4

/*
 * decgrain - granularity of the declination bins, ie the bin width, in degrees
 */
static void beam_dec_dependence(FluxWappData * wappdata, int beam, float decmin, float decmax, float decgrain)
{
	int i, d, r;
	int num_bins;
	double *binI, *binQ, *binU, *binV, *binDEC;
	int * bin_counts;
	FILE * decfile;
	char filename[16+1];
	
	snprintf(filename, 16, "dec_beam%i.dat", beam);
	decfile = fopen(filename, "w");
	

	num_bins = (int) ceil((decmax-decmin)/decgrain);
	decgrain = (decmax-decmin)/num_bins; //actual grain (bin) size

	binI = calloc(num_bins, sizeof(double));
	binQ = calloc(num_bins, sizeof(double));
	binU = calloc(num_bins, sizeof(double));
	binV = calloc(num_bins, sizeof(double));
	binDEC = calloc(num_bins, sizeof(double));
	bin_counts = calloc(num_bins, sizeof(int));

	//TODO: reject outliers as needed

	for (d=0; d<wappdata->numDays; d++) 
	{
		if (d%7!=beam) continue; //TODO: this is a workaround for non-beam aware datastructures
		FluxDayData * daydata;

		daydata = &wappdata->daydata[d];
		for (r=0; r<daydata->numRecords; r++) 
		{
			int bin;
			bin = (int) floor((daydata->records[r].DEC - decmin) / decgrain);
			if (bin<0 || bin>=num_bins) continue;
			binI[bin] += daydata->records[r].stokes.I;
			binQ[bin] += daydata->records[r].stokes.Q;
			binU[bin] += daydata->records[r].stokes.U;
			binV[bin] += daydata->records[r].stokes.V;
			bin_counts[bin]++;
		}
	}

	//normalize by the bin counts
	for (i=0; i<num_bins; i++) 
	{
		register int count = bin_counts[i];
		binI[i] /= count;
		binQ[i] /= count;
		binU[i] /= count;
		binV[i] /= count;
		binDEC[i] = i*decgrain+decmin; //TODO: compute this more intellegently
		fprintf(decfile, "%f %i %g %g %g %g\n", binDEC[i], count, binI[i], binQ[i], binU[i], binV[i]); 
	}

	fclose(decfile);
	

	//curve fit the binned data
	float nsigma = 5.0;
	double chisq;
	double cI[BIN_ORDER+1];
	double cQ[BIN_ORDER+1];
	double cU[BIN_ORDER+1];
	double cV[BIN_ORDER+1];
	double min, max;

	jsd_minmax(binDEC, num_bins, &min, &max);
	jsd_normalize(binDEC, num_bins, min, max);
//printf("decmin:%g min:%g decmax:%g, max: %g\n", decmin, min, decmax, max);

	jsd_poly_fit(binDEC, binI, num_bins, nsigma, cI, BIN_ORDER, &chisq);
	jsd_poly_fit(binDEC, binQ, num_bins, nsigma, cQ, BIN_ORDER, &chisq);
	jsd_poly_fit(binDEC, binU, num_bins, nsigma, cU, BIN_ORDER, &chisq);
	jsd_poly_fit(binDEC, binV, num_bins, nsigma, cV, BIN_ORDER, &chisq);

	free(binI);
	free(binQ);
	free(binU);
	free(binV);
	free(binDEC);
	free(bin_counts);

	for (d=0; d<wappdata->numDays; d++) 
	{
		if (d%7!=beam) continue; //TODO: this is a workaround for non-beam aware datastructures
		FluxDayData * daydata;

		daydata = &wappdata->daydata[d];
		for (r=0; r<daydata->numRecords; r++) 
		{
			double DEC = NORMALIZE(daydata->records[r].DEC, min, max);
			daydata->records[r].stokes.I -= jsd_poly_eval(DEC, cI, BIN_ORDER);
			daydata->records[r].stokes.Q -= jsd_poly_eval(DEC, cQ, BIN_ORDER);
			daydata->records[r].stokes.U -= jsd_poly_eval(DEC, cU, BIN_ORDER);
			daydata->records[r].stokes.V -= jsd_poly_eval(DEC, cV, BIN_ORDER);
		}
//if (d==0) jsd_print_poly(stdout, cI, BIN_ORDER);
	}
}


void remove_dec_dependence(FluxWappData * wappdata, float decmin, float decmax, float decgrain)
{
	int beam;

	for (beam=0; beam<7; beam++) {
		beam_dec_dependence(wappdata, beam, decmin, decmax, decgrain);
	}
}


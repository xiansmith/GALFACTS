#include "calibrate.h"
#include "rfi.h"
#include <stdlib.h>
#include "jsd/jsd_fit.h"
#include <math.h>
#include "smoothing.h"
#include "stats.h"
#include <string.h>

//pre-compute the raw cal for every channel
//void compute_raw_cal(SpecRecord dataset[], int size)
void compute_raw_cal(SpecRecord dataset[], int size, int lowchan, int highchan)
{
	int n, i;
	for (n=0; n<size; n++)
	{
		for (i=lowchan; i<highchan; i++)
		{
			dataset[n].cal.xx[i] = dataset[n].calon.xx[i] - dataset[n].caloff.xx[i];
			dataset[n].cal.yy[i] = dataset[n].calon.yy[i] - dataset[n].caloff.yy[i];
			dataset[n].cal.xy[i] = dataset[n].calon.xy[i] - dataset[n].caloff.xy[i];
			dataset[n].cal.yx[i] = dataset[n].calon.yx[i] - dataset[n].caloff.yx[i];
		}
	}
}



/* Averages up the cal in frequency, does a curve fit, and then applies it.
 * TODO: the curve fits fail on the typical data.  I hate NR!
 */
/*#define FIT_ORDER (8)
void curve_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI)
{
	int n, chan;
	double C[FIT_ORDER+1];
	float nsigma = 4.0;
	double sumsq;
	double *Xast, *Txx, *Tyy, *Txy, *Tyx;
	double *Sxx, *Syy, *Sxy, *Syx;
	double *Kxx, *Kyy, *Kxy, *Kyx;
	double Txxavg, Tyyavg, Txyavg, Tyxavg;
	FILE * eqfile;
	FILE * chifile;
	size_t count;

	eqfile = fopen("caleq.dat", "w");
	fprintf(eqfile, "#XXc0 XXc1 YYc0 YYc1 XYc0 XYc1 YXc0 YXc1\n");

	chifile = fopen("calchi.dat", "w");
	fprintf(chifile, "#XX YY XY YX\n");

	Xast = (double*) calloc(size, sizeof(double));
	Txx = (double*) calloc(size, sizeof(double));
	Tyy = (double*) calloc(size, sizeof(double));
	Txy = (double*) calloc(size, sizeof(double));
	Tyx = (double*) calloc(size, sizeof(double));
	Sxx = (double*) calloc(MAX_CHANNELS, sizeof(double));
	Syy = (double*) calloc(MAX_CHANNELS, sizeof(double));
	Sxy = (double*) calloc(MAX_CHANNELS, sizeof(double));
	Syx = (double*) calloc(MAX_CHANNELS, sizeof(double));
	Kxx = (double*) calloc(MAX_CHANNELS, sizeof(double));
	Kyy = (double*) calloc(MAX_CHANNELS, sizeof(double));
	Kxy = (double*) calloc(MAX_CHANNELS, sizeof(double));
	Kyx = (double*) calloc(MAX_CHANNELS, sizeof(double));


	//average up the cals in time (T__)
	for (n=0; n<size; n++) 
	{
		Xast[n] = dataset[n].AST;
		if (dataset[n].flagBAD) {
			Txx[n] = Tyy[n] = Tyx[n] = Txy[n] = NAN;
			continue;
		}
		count = 0;
		for (chan=lowchan; chan<=highchan; chan++) 
		{
			if (ignoreRFI || dataset[n].flagRFI[chan] == RFI_NONE) {
				Txx[n] += dataset[n].cal.xx[chan];
				Tyy[n] += dataset[n].cal.yy[chan];
				Txy[n] += dataset[n].cal.xy[chan];
				Tyx[n] += dataset[n].cal.yx[chan];
				count++;
			}
		}

		Txx[n] /= count;
		Tyy[n] /= count;
		Txy[n] /= count;
		Tyx[n] /= count;

	}

	//compute the scalar average tcals 
	Txxavg = compute_mean(Txx, 0,size);
	Tyyavg = compute_mean(Tyy, 0,size);
	Tyxavg = compute_mean(Tyx, 0,size);
	Txyavg = compute_mean(Txy, 0,size);
	printf("Txxavg:%g  Tyyavg:%g  Txyavg:%g  Tyxavg:%g \n", Txxavg, Tyyavg, Txyavg, Tyxavg);


	//average up the tcals in spectral channel (S__)
	for (chan=lowchan; chan<=highchan; chan++) 
	{
		count = 0;
		for (n=0; n<size; n++) 
		{
			if (dataset[n].flagBAD) continue;
			if (ignoreRFI || dataset[n].flagRFI[chan] == RFI_NONE) {
				Sxx[chan] += dataset[n].cal.xx[chan];
				Syy[chan] += dataset[n].cal.yy[chan];
				Sxy[chan] += dataset[n].cal.xy[chan];
				Syx[chan] += dataset[n].cal.yx[chan];
				count++;
			}
		}
		Sxx[chan] /= count;
		Syy[chan] /= count;
		Sxy[chan] /= count;
		Syx[chan] /= count;
	}

	//compute the correction array using the scalar average tcals and the average spectra
	for (chan=lowchan; chan<=highchan; chan++) {
		Kxx[chan] = Txxavg - Sxx[chan];
		Kyy[chan] = Tyyavg - Syy[chan];
		Kxy[chan] = Txyavg - Sxy[chan];
		Kyx[chan] = Tyxavg - Syx[chan];
	}

	{
		FILE * testfile;
		testfile = fopen("testfile.dat", "w");
		for (chan=lowchan; chan<=highchan; chan++) {
			fprintf(testfile, "%i %g %g\n", chan, Sxx[chan], Kxx[chan]);
		}
		//for (n=0; n<size; n++) {
		//	fprintf(testfile, "%g %g\n", Xast[n], Txx[n]);
		//}

		fclose(testfile);
	}


	// normalize the x values for the curve fit 
	double min, max;
	jsd_minmax(Xast, size, &min, &max);
	jsd_normalize(Xast, size, min, max);

	//XX curve fit and apply
	jsd_poly_fit(Xast, Txx, size, nsigma, C, FIT_ORDER, &sumsq);
	fprintf (eqfile, "XX: ");
	jsd_print_poly (eqfile, C, FIT_ORDER);
	fprintf(chifile, "%g ", sumsq);
	for (n=0; n<size; n++) {
		double val = jsd_poly_eval(NORMALIZE(Xast[n], min, max), C, FIT_ORDER);
		for (chan=lowchan; chan<=highchan; chan++) {
			dataset[n].cal.xx[chan] = val + Kxx[chan];
		}
	}

	//YY curve fit and apply
	jsd_poly_fit(Xast, Tyy, size, nsigma, C, FIT_ORDER, &sumsq);
	fprintf (eqfile, "YY: ");
	jsd_print_poly (eqfile, C, FIT_ORDER);
	fprintf(chifile, "%g ", sumsq);
	for (n=0; n<size; n++) {
		double val = jsd_poly_eval(NORMALIZE(Xast[n], min, max), C, FIT_ORDER);
		for (chan=lowchan; chan<=highchan; chan++) {
			dataset[n].cal.yy[chan] = val + Kyy[chan];
		}
	}

	//XY curve fit and apply
	jsd_poly_fit(Xast, Txy, size, nsigma, C, FIT_ORDER, &sumsq);
	fprintf (eqfile, "XY: ");
	jsd_print_poly (eqfile, C, FIT_ORDER);
	fprintf(chifile, "%g ", sumsq);
	for (n=0; n<size; n++) {
		double val = jsd_poly_eval(NORMALIZE(Xast[n], min, max), C, FIT_ORDER);
		for (chan=lowchan; chan<=highchan; chan++) {
			dataset[n].cal.xy[chan] = val + Kxy[chan];
		}
	}

	//YX curve fit and apply
	jsd_poly_fit(Xast, Tyx, size, nsigma, C, FIT_ORDER, &sumsq);
	fprintf (eqfile, "YX: ");
	jsd_print_poly (eqfile, C, FIT_ORDER);
	fprintf(chifile, "%g ", sumsq);
	for (n=0; n<size; n++) {
		double val = jsd_poly_eval(NORMALIZE(Xast[n], min, max), C, FIT_ORDER);
		for (chan=lowchan; chan<=highchan; chan++) {
			dataset[n].cal.yx[chan] = val + Kyx[chan];
		}
	}

	fprintf(chifile, "\n");


	fclose(eqfile);
	fclose(chifile);
	free(Xast);
	free(Txx);
	free(Tyy);
	free(Txy);
	free(Tyx);
	free(Sxx);
	free(Syy);
	free(Sxy);
	free(Syx);
	free(Kxx);
	free(Kyy);
	free(Kxy);
	free(Kyx);

}
*/
/*
 * Normalizes the data array by its average.
 */
static void normalize_data(double * data, int start, int end)
{
	double avg;
	int i;

	avg = compute_mean(data, start, end);
	for (i=start; i<end; i++) {
		data[i] /= avg;
	}
}



static void compute_avg_spectra(SpecRecord dataset[], double * xxf, double * yyf, double * xyf, double * yxf, int lowchan, int highchan, int tstart, int tend)
{
	int count;
	int i, chan;


	//average up the cal in time to produce an average cal spectra, 
	memset(xxf, 0, highchan*sizeof(double));
	memset(yyf, 0, highchan*sizeof(double));
	memset(xyf, 0, highchan*sizeof(double));
	memset(yxf, 0, highchan*sizeof(double));
//	for (chan=0; chan<MAX_CHANNELS; chan++)
	for (chan=lowchan; chan<highchan; chan++) 
	{
		count = 0;
		for (i=tstart; i<tend; i++) 
		{
			if ((dataset[i].flagBAD || dataset[i].flagRFI[chan] != RFI_NONE)) {
				continue;
			}
			count++;
			xxf[chan] += dataset[i].cal.xx[chan]; 
			yyf[chan] += dataset[i].cal.yy[chan]; 
			xyf[chan] += dataset[i].cal.xy[chan]; 
			yxf[chan] += dataset[i].cal.yx[chan]; 
		}
		xxf[chan] /= count;
		yyf[chan] /= count;
		xyf[chan] /= count;
		yxf[chan] /= count;
	}
}

/*
 * Averages up the cal across the band, for all of the auto and cross correlations,
 * performs hanning smoothing of the given smooth_width and then writes the
 * smoothed values back into dataset[i].cal.  Writes out a text file (calfile.dat)
 * that contains the smooth values for plotting purposes.
 *
 * - raw cals must be pre-computed in dataset[i].cal before this function is called.
 * - Only channels between lowchan and highchan are used in the average.  
 * - nsigma is used in the outlier rejection step before smoothing is performed.
 * - smooth_width must be given in bins (not seconds) and should be odd.
 * 
 */
void smooth_cal_bandaverage(SpecRecord dataset[], int size, int lowchan, int highchan, int t_smooth_width, float nsigma)
{
//	const int f_smooth_width = 3;
	double *xx, *yy, *xy, *yx, *ast;
	double xxf[MAX_CHANNELS];
	double yyf[MAX_CHANNELS];
	double xyf[MAX_CHANNELS];
	double yxf[MAX_CHANNELS];
	int count;
	int i, chan;
	const char *filename = "calfile.dat";
	FILE *calfile;
	printf("Requesting malloc for %u bytes of memory\n",sizeof(double)*size);
	xx = (double*) calloc(size, sizeof(double));
	if (xx == NULL) {
		printf("ERROR: malloc failed in smooth_cal_bandaverage() !\n");
	}
	printf("Requesting malloc for %u bytes of memory\n",sizeof(double)*size);
	yy = (double*) calloc(size, sizeof(double));
	if (yy == NULL) {
		printf("ERROR: malloc failed in smooth_cal_bandaverage() !\n");
	}
	printf("Requesting malloc for %u bytes of memory\n",sizeof(double)*size);
	xy = (double*) calloc(size, sizeof(double));
	if (xy == NULL) {
		printf("ERROR: malloc failed in smooth_cal_bandaverage() !\n");
	}
	printf("Requesting malloc for %u bytes of memory\n",sizeof(double)*size);
	yx = (double*) calloc(size, sizeof(double));
	if (yx == NULL) {
		printf("ERROR: malloc failed in smooth_cal_bandaverage() !\n");
	}
	printf("Requesting malloc for %u bytes of memory\n",sizeof(double)*size);
	ast = (double*) calloc(size, sizeof(double));
	if (ast == NULL) {
		printf("ERROR: malloc failed in smooth_cal_bandaverage() !\n");
	}
	
	//average up the cal in frequency to produce an average cal amplitude
	for (i=0; i<size; i++) 
	{
		if (dataset[i].flagBAD) {
			xx[i] = yy[i] = xy[i] = yx[i] = NAN;
			continue;
		}
		count = 0;
//		for (chan=0; chan<MAX_CHANNELS; chan++)
		for (chan=lowchan; chan<highchan; chan++) 
		{
			if ((dataset[i].flagRFI[chan] != RFI_NONE)) {
				continue;
			}
			count++;
			xx[i] += dataset[i].cal.xx[chan]; 
			yy[i] += dataset[i].cal.yy[chan]; 
			xy[i] += dataset[i].cal.xy[chan]; 
			yx[i] += dataset[i].cal.yx[chan]; 
		}
		xx[i] /= count; 
		yy[i] /= count; 
		xy[i] /= count; 
		yx[i] /= count; 
	}

	//do outlier rejection to kill the stupid highcals at endpoints
	reject_outliers(xx, size, nsigma);
	reject_outliers(xy, size, nsigma);
	reject_outliers(yx, size, nsigma);
	reject_outliers(yy, size, nsigma);

	//smooth the band average amplitudes
	hanning_smooth_data(xx, size, t_smooth_width);
	hanning_smooth_data(xy, size, t_smooth_width);
	hanning_smooth_data(yx, size, t_smooth_width);
	hanning_smooth_data(yy, size, t_smooth_width);


	//determine the average spectra
//	compute_avg_spectra(dataset, xxf, yyf, xyf, yxf, 0, MAX_CHANNELS, 0, size);
	compute_avg_spectra(dataset, xxf, yyf, xyf, yxf, lowchan, highchan, 0, size);
	//smooth the data values in frequency
	//hanning_smooth_data(xxf, MAX_CHANNELS, f_smooth_width);
	//hanning_smooth_data(xyf, MAX_CHANNELS, f_smooth_width);
	//hanning_smooth_data(yxf, MAX_CHANNELS, f_smooth_width);
	//hanning_smooth_data(yyf, MAX_CHANNELS, f_smooth_width);

	//normalize the data values to 1
//	normalize_data(xxf, MAX_CHANNELS);
//	normalize_data(yyf, MAX_CHANNELS);
//	normalize_data(xyf, MAX_CHANNELS);
//	normalize_data(yxf, MAX_CHANNELS);
	normalize_data(xxf, lowchan, highchan);
	normalize_data(yyf, lowchan, highchan);
	normalize_data(xyf, lowchan, highchan);
	normalize_data(yxf, lowchan, highchan);
	
	{
		double xxf1[MAX_CHANNELS], yyf1[MAX_CHANNELS], xyf1[MAX_CHANNELS], yxf1[MAX_CHANNELS];
		double xxf2[MAX_CHANNELS], yyf2[MAX_CHANNELS], xyf2[MAX_CHANNELS], yxf2[MAX_CHANNELS];
		double xxf3[MAX_CHANNELS], yyf3[MAX_CHANNELS], xyf3[MAX_CHANNELS], yxf3[MAX_CHANNELS];
		double xxf4[MAX_CHANNELS], yyf4[MAX_CHANNELS], xyf4[MAX_CHANNELS], yxf4[MAX_CHANNELS];
		FILE * caldat;

//		compute_avg_spectra(dataset, xxf1, yyf1, xyf1, yxf1, 0, MAX_CHANNELS, 0, size/4);
//		compute_avg_spectra(dataset, xxf2, yyf2, xyf2, yxf2, 0, MAX_CHANNELS, size/4, 2*size/4);
//		compute_avg_spectra(dataset, xxf3, yyf3, xyf3, yxf3, 0, MAX_CHANNELS, 2*size/4, 3*size/4);
//		compute_avg_spectra(dataset, xxf4, yyf4, xyf4, yxf4, 0, MAX_CHANNELS, 3*size/4, size);
		compute_avg_spectra(dataset, xxf1, yyf1, xyf1, yxf1, lowchan, highchan, 0, size/4);
		compute_avg_spectra(dataset, xxf2, yyf2, xyf2, yxf2, lowchan, highchan, size/4, 2*size/4);
		compute_avg_spectra(dataset, xxf3, yyf3, xyf3, yxf3, lowchan, highchan, 2*size/4, 3*size/4);
		compute_avg_spectra(dataset, xxf4, yyf4, xyf4, yxf4, lowchan, highchan, 3*size/4, size);

/*		normalize_data(xxf1, MAX_CHANNELS);
		normalize_data(xxf2, MAX_CHANNELS);
		normalize_data(xxf3, MAX_CHANNELS);
		normalize_data(xxf4, MAX_CHANNELS);

		normalize_data(yyf1, MAX_CHANNELS);
		normalize_data(yyf2, MAX_CHANNELS);
		normalize_data(yyf3, MAX_CHANNELS);
		normalize_data(yyf4, MAX_CHANNELS);
*/
		normalize_data(xxf1, lowchan, highchan);
		normalize_data(xxf2, lowchan, highchan);
		normalize_data(xxf3, lowchan, highchan);
		normalize_data(xxf4, lowchan, highchan);

		normalize_data(yyf1, lowchan, highchan);
		normalize_data(yyf2, lowchan, highchan);
		normalize_data(yyf3, lowchan, highchan);
		normalize_data(yyf4, lowchan, highchan);

		caldat = fopen("calspectra.dat", "w");
		fprintf(caldat, "#xxf xxf1 xxf2 xxf3 xxf4\n");
		for (i=lowchan; i<highchan; i++) 
		{
			fprintf(caldat, "%g %g %g %g %g %g %g %g %g %g\n", 
				xxf[i], xxf1[i], xxf2[i], xxf3[i], xxf4[i], yyf[i], yyf1[i], yyf2[i], yyf3[i], yyf4[i]);
		}
		fclose(caldat);

	}


	//Compute the cal as the scaled band average
	for (i=0; i<size; i++) 
	{
		//store these new smoothed values back in to the cal arrays of the dataset
		for (chan=lowchan; chan<highchan; chan++) 
		{
			dataset[i].cal.xx[chan] = xxf[chan] * xx[i];
			dataset[i].cal.yy[chan] = yyf[chan] * yy[i];
			dataset[i].cal.xy[chan] = xyf[chan] * xy[i];
			dataset[i].cal.yx[chan] = yxf[chan] * yx[i];
		}
	}

	// write out the cal file
	calfile = fopen(filename, "w");
	if (calfile == NULL) {
		printf("ERROR: unable to open file %s\n", filename);
		return;
	}
	fprintf(calfile, "# RA DEC AST calXX calYY calXY calYX\n");
	for (i=0; i<size; i++) 
	{
		fprintf(calfile,"%2.8f %2.8f %8.2f %4.8f %4.8f %4.8f %4.8f\n", 
				dataset[i].RA, dataset[i].DEC, dataset[i].AST, xx[i], yy[i], xy[i], yx[i]);
	}
	fclose(calfile);

	free(ast);
	free(xx);
	free(yy);
	free(xy);
	free(yx);

}

/*
void linear_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI)
{
	int n, chan;
	double C[2]; //assuming linear fit 
	float nsigma = 4.0;
	double sumsq;
	double *Xast, *Yxx, *Yyy, *Yxy, *Yyx;
	FILE * eqfile;
	FILE * chifile;
	size_t count;

	eqfile = fopen("caleq.dat", "w");
	fprintf(eqfile, "#chan XXc0 XXc1 YYc0 YYc1 XYc0 XYc1 YXc0 YXc1\n");

	chifile = fopen("calchi.dat", "w");
	fprintf(chifile, "#chan XX YY XY YX\n");

	//y is the cal value
	//x is the time steps
	Xast = (double*) malloc(sizeof(double) * size);
	Yxx = (double*) malloc(sizeof(double) * size);
	Yyy = (double*) malloc(sizeof(double) * size);
	Yxy = (double*) malloc(sizeof(double) * size);
	Yyx = (double*) malloc(sizeof(double) * size);

	for (chan=lowchan; chan<=highchan; chan++) 
	{
		fprintf(chifile, "%i ", chan);

		count = 0;
		for (n=0; n<size; n++) {
			if (dataset[n].flagBAD) continue;
			if (ignoreRFI || dataset[n].flagRFI[chan] == RFI_NONE) {
				Xast[count] = dataset[n].AST;
				Yxx[count] = dataset[n].cal.xx[chan];
				Yyy[count] = dataset[n].cal.yy[chan];
				Yxy[count] = dataset[n].cal.xy[chan];
				Yyx[count] = dataset[n].cal.yx[chan];
				count++;
			}
		}

		if (count == 0) {
			//NaN values will result in the dataset cal
		}

		// normalize the x values for the curve fit 
		{
			double min, max;
			jsd_minmax(Xast, count, &min, &max);
			jsd_normalize(Xast, count, min, max);
		}

		//XX
		jsd_linear_fit(Xast, Yxx, count, nsigma, C, &sumsq);
		fprintf (eqfile, "XX %i %g %g\n", chan, C[0], C[1]);
		fprintf(chifile, "%g ", sumsq);
		for (n=0; n<size; n++) //eval function
			dataset[n].cal.xx[chan] = jsd_linear_eval(Xast[n], C);

		//YY
		jsd_linear_fit(Xast, Yyy, count, nsigma, C, &sumsq);
		fprintf (eqfile, "YY %i %g %g\n", chan, C[0], C[1]);
		fprintf(chifile, "%g ", sumsq);
		for (n=0; n<size; n++) //eval function
			dataset[n].cal.yy[chan] = jsd_linear_eval(Xast[n], C);

		//XY
		jsd_linear_fit(Xast, Yxy, count, nsigma, C, &sumsq);
		fprintf (eqfile, "XY %i %g %g\n", chan, C[0], C[1]);
		fprintf(chifile, "%g ", sumsq);
		for (n=0; n<size; n++) //eval function
			dataset[n].cal.xy[chan] = jsd_linear_eval(Xast[n], C);

		//YX
		jsd_linear_fit(Xast, Yyx, count, nsigma, C, &sumsq);
		fprintf (eqfile, "YX %i %g %g\n", chan, C[0], C[1]);
		fprintf(chifile, "%g ", sumsq);
		for (n=0; n<size; n++) //eval function
			dataset[n].cal.yx[chan] = jsd_linear_eval(Xast[n], C);

		fprintf(chifile, "\n");

	}

	fclose(eqfile);
	fclose(chifile);
	free(Xast);
	free(Yxx);
	free(Yyy);
	free(Yxy);
	free(Yyx);
}

*/


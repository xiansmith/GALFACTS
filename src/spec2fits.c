#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>
#include "common.h"
#include "spec.h"
#include "smooth.h"
#include "rfi.h"
#include "programs/fitsLib.h"
#include "markdata.h"
#include <math.h>


#define FILENAME_SIZE 128
static float compute_mean(float data[], int size)
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

	return (float) sum/count;
}


#define SQR(X) ((X)*(X))
static float compute_sigma(float data[], int size, float avg)
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
	return (float) sqrt(sum / (count-1) );
}

/* Returns the error string for the given errcode and compiled regular expression.
 * Returned string is dynamicly allocated and must be freed by the caller.
 */
static char *get_regerror (int errcode, regex_t *preg)
{
	size_t length = regerror(errcode, preg, NULL, 0);
	char *buffer = malloc(length);
	regerror(errcode, preg, buffer, length);
	return buffer;
}


/* Gets the regex error and writes to the given out file.
 */
static void handle_regerr(int errcode, regex_t *preg, FILE *out)
{
	char * str = get_regerror(errcode, preg);
	fprintf(out, "ERROR: %s\n", str);
	free(str);
}


/* Copy the substring of src from start to end into dest.  Ensures that no more than
 * max characters (minus 1 for null) are copied.
 */
static void strnsub (char * dest, const char * src, int start, int end, int max)
{
	int num;
	if ((int)strlen(src) < start) return; //start is off the end of the string
	num = (end-start > max) ? max : end-start; //ensure we don't copy in too many
	strncpy(dest, &src[start], num); //do the copy
	dest[num] = '\0';
}


/* get the filename from the file path out of the filename
 */
static int get_name_from_filepath(char * name, char * path, int max, const char * filepath)
{
	regex_t reg;
	regex_t * preg = &reg;
	int errcode;
	const char * regex = "(.*)/([^/]*)$";
	regmatch_t pmatch[3];
	size_t nmatch = 3;
	errcode = regcomp(preg, regex, REG_EXTENDED);
	if (errcode != 0) {
		handle_regerr(errcode, preg, stderr);
		return -1;
	}

	errcode = regexec(preg, filepath, nmatch, pmatch, 0x0);
	if (errcode != 0) {
		handle_regerr(errcode, preg, stderr);
		return -2;
	}

	regfree(preg);

	//pmatch[0] is the entire match, pmatch[1] is the first subexpression match in '()'
	if (path != NULL) strnsub(path, filepath, pmatch[1].rm_so, pmatch[1].rm_eo, max);
	if (name != NULL) strnsub(name, filepath, pmatch[2].rm_so, pmatch[2].rm_eo, max);

	return 0;
}


/*
 * Writes pointing information for the observations.   Bad flags ARE observed.
 */
static void write_pointing(SpecRecord dataset[], int size, const char * fileroot)
{
	int n;
	FILE * file;
	char filename[FILENAME_SIZE+1];

	snprintf(filename, FILENAME_SIZE, "%s_pointing.dat", fileroot);
	file = fopen(filename, "w");
	if (file == NULL) {
		printf("ERROR: unable to open '%s' for writing\n", filename);
		return;
	}

	fprintf(file, "#RA DEC AST\n");
	for (n=0; n<size; n++)
	{
		if (dataset[n].flagBAD) continue;

		fprintf(file, "%7.6f %7.6f %7.2f\n",
			dataset[n].RA/15.0, dataset[n].DEC, dataset[n].AST);
	}

	fclose(file);
}


/*
 * Writes band averaged information for the observations.   Bad flags are NOT observed.
 */
static void write_band_average(SpecRecord dataset[], int size, const char * fileroot)
{
	int n, chan, count;
	FILE * file;
	char filename[FILENAME_SIZE+1];
	PolAvg avgcalon, avgcaloff;

	snprintf(filename, FILENAME_SIZE, "%s_bandavg.dat", fileroot);
	file = fopen(filename, "w");
	if (file == NULL) {
		printf("ERROR: unable to open '%s' for writing\n", filename);
		return;
	}

	fprintf(file, "#RA DEC AST CALON(XX XY YX YY) CALOFF(XX XY YX YY)\n");
	for (n=0; n<size; n++)
	{
		memset(&avgcalon, 0, sizeof(PolAvg));
		memset(&avgcaloff, 0, sizeof(PolAvg));
		count = 0;
		for (chan=0; chan<MAX_CHANNELS; chan++)
		{
			avgcalon.xx += dataset[n].calon.xx[chan];
			avgcalon.xy += dataset[n].calon.xy[chan];
			avgcalon.yx += dataset[n].calon.yx[chan];
			avgcalon.yy += dataset[n].calon.yy[chan];
			avgcaloff.xx += dataset[n].caloff.xx[chan];
			avgcaloff.xy += dataset[n].caloff.xy[chan];
			avgcaloff.yx += dataset[n].caloff.yx[chan];
			avgcaloff.yy += dataset[n].caloff.yy[chan];
			count++;
		}
		if (count <= 0) continue;
		fprintf(file, "%7.6f %7.6f %7.2f %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f\n",
			dataset[n].RA/15.0, dataset[n].DEC, dataset[n].AST,
			avgcalon.xx/chan, avgcalon.xy/chan, avgcalon.yx/chan, avgcalon.yy/chan,
			avgcaloff.xx/chan, avgcaloff.xy/chan, avgcaloff.yx/chan, avgcaloff.yy/chan );
	}

	fclose(file);
}

static void write_time_average(SpecRecord dataset[], int size, const char * fileroot)
{
	int n, chan;
	FILE * file;
	char filename[FILENAME_SIZE+1];
	PolAvg avgcalon, avgcaloff;

	snprintf(filename, FILENAME_SIZE, "%s_timeavg.dat", fileroot);
	file = fopen(filename, "w");
	if (file == NULL) {
		printf("ERROR: unable to open '%s' for writing\n", filename);
		return;
	}

	fprintf(file, "#CALON(XX XY YX YY) CALOFF(XX XY YX YY)\n");
	for (chan=0; chan<MAX_CHANNELS; chan++)
	{
		int count = 0;
		memset(&avgcalon, 0, sizeof(PolAvg));
		memset(&avgcaloff, 0, sizeof(PolAvg));
		for (n=0; n<size; n++)
		{
			avgcalon.xx += dataset[n].calon.xx[chan];
			avgcalon.xy += dataset[n].calon.xy[chan];
			avgcalon.yx += dataset[n].calon.yx[chan];
			avgcalon.yy += dataset[n].calon.yy[chan];
			avgcaloff.xx += dataset[n].caloff.xx[chan];
			avgcaloff.xy += dataset[n].caloff.xy[chan];
			avgcaloff.yx += dataset[n].caloff.yx[chan];
			avgcaloff.yy += dataset[n].caloff.yy[chan];
			count++;
		}
		fprintf(file, "%7.6f %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f\n",
			avgcalon.xx/count, avgcalon.xy/count, avgcalon.yx/count, avgcalon.yy/count,
			avgcaloff.xx/count, avgcaloff.xy/count, avgcaloff.yx/count, avgcaloff.yy/count);
	}

	fclose(file);
}

static void write_rfi_data(SpecRecord dataset[], int size, const char * fileroot, float fcen, float df)
{
	float numSigma = 4.5;
	float numSigmaThresh = 100.0;
	int ignoreA_low = 61;
	int ignoreA_high = 64;
	int lowchan = 25;
	int highchan = 230;
	int n, chan;
	float fstart;
	char filename[FILENAME_SIZE+1];
	FILE * file;

	snprintf(filename, FILENAME_SIZE, "%s_rfi.dat", fileroot);
	file = fopen(filename, "w");
	if (file == NULL) {
		printf("ERROR: unable to open '%s' for writing\n", filename);
		return;
	}
	fprintf(file, "#chan freq RA DEC AST\n");

	fstart = fcen - MAX_CHANNELS/2 * df;

	rfi_detection(dataset, size, lowchan, highchan, numSigma, numSigmaThresh, ignoreA_low, ignoreA_high, 0, 0);
	aerostat_rfi_blanking(dataset, size, lowchan, highchan);

	for (n=0; n<size; n++)
	{
		SpecRecord * pRec = &(dataset[n]);

		//now that the channels are marked as RFI, write them out
		for (chan=0; chan<MAX_CHANNELS; chan++) {
		 	if (pRec->flagRFI[chan] != RFI_NONE) {
				fprintf(file,
					"%3i %9.6f %9.6f %9.6f %8.2f\n",
					chan, fstart+df*chan, pRec->RA, pRec->DEC, pRec->AST);
			}
		}
	}

	fclose (file);
}

/* Takes a full days worth of data as dataset with size number of measurements
 * and then writes out raw fits files of all the data.  Two cubes are produced,
 * one for calon and one for cal off.  The four planes are the 4 auto/cross correlations.
 */
static void write_raw_pol_fits(SpecRecord dataset[], int size, float fcen, float df, const char * fileroot)
{
	int n, chan;
	header_param_list hpar;
	float * data;
	FILE * file;
	char filename[FILENAME_SIZE+1];

	//allocate memory for one plane only to minimize further memory requirements
	data = (float *) malloc (MAX_CHANNELS * size * sizeof(float));
	if (data == NULL) {
		printf("ERROR: out of memory!\n");
		return;
	}

	/* initilize the header parameter object for the cube */
	init_header_param_list (&hpar);  /* initialize parameter records */
	hpar.bitpix = -32;
	hpar.num_axes = 3;

	hpar.naxis[0] = MAX_CHANNELS;
	hpar.naxis[1] = size;
	hpar.naxis[2] = 4;

	sprintf (hpar.ctype[0], "Frequency");
	sprintf (hpar.ctype[1], "Time");
	sprintf (hpar.ctype[2], "Correlation");

	hpar.crval[0] = fcen;		/* ref channel value*/
	hpar.crval[1] = dataset[0].AST;	/* ref time value */
	hpar.crval[2] = 0;		/* ref correlation value */

	hpar.crpix[0] = MAX_CHANNELS / 2;	/* start channel position */
	hpar.crpix[1] = 0;			/* start time position */
	hpar.crpix[2] = 0;			/* correlation start position */

	hpar.cdelt[0] = df;						/* channel increment */
	hpar.cdelt[1] = (dataset[size-1].AST - dataset[0].AST) / size;	/* time per increment */
	hpar.cdelt[2] = 1;						/* correlation increment*/

	sprintf (hpar.bunit, "Power");
	sprintf (hpar.telescope, "Arecibo");


	/* CAL OFF */
	sprintf (hpar.object, "Cal Off XX, XY, YX, YY");
	snprintf(filename, FILENAME_SIZE, "%s_caloff.fits", fileroot);
	file = fopen(filename, "w");
	if (file == NULL) {
		printf("ERROR: unable to open file '%s' for writing\n", filename);
	}
	writefits_header(file, &hpar);

	for (n=0; n < size; n++)
		for (chan=0; chan<MAX_CHANNELS; chan++)
			data[n*MAX_CHANNELS+chan] = dataset[n].caloff.xx[chan];
	writefits_plane(file, data, &hpar);
	for (n=0; n < size; n++)
		for (chan=0; chan<MAX_CHANNELS; chan++)
			data[n*MAX_CHANNELS+chan] = dataset[n].caloff.xy[chan];
	writefits_plane(file, data, &hpar);
	for (n=0; n < size; n++)
		for (chan=0; chan<MAX_CHANNELS; chan++)
			data[n*MAX_CHANNELS+chan] = dataset[n].caloff.yx[chan];
	writefits_plane(file, data, &hpar);
	for (n=0; n < size; n++)
		for (chan=0; chan<MAX_CHANNELS; chan++)
			data[n*MAX_CHANNELS+chan] = dataset[n].caloff.yy[chan];
	writefits_plane(file, data, &hpar);

	writefits_pad_end(file, &hpar);
	fclose(file);


	/* CAL ON */
	sprintf (hpar.object, "Cal On XX, XY, YX, YY");
	snprintf(filename, FILENAME_SIZE, "%s_calon.fits", fileroot);
	file = fopen(filename, "w");
	if (file == NULL) {
		printf("ERROR: unable to open file '%s' for writing\n", filename);
	}
	writefits_header(file, &hpar);

	for (n=0; n < size; n++)
		for (chan=0; chan<MAX_CHANNELS; chan++)
			data[n*MAX_CHANNELS+chan] = dataset[n].calon.xx[chan];
	writefits_plane(file, data, &hpar);
	for (n=0; n < size; n++)
		for (chan=0; chan<MAX_CHANNELS; chan++)
			data[n*MAX_CHANNELS+chan] = dataset[n].calon.xy[chan];
	writefits_plane(file, data, &hpar);
	for (n=0; n < size; n++)
		for (chan=0; chan<MAX_CHANNELS; chan++)
			data[n*MAX_CHANNELS+chan] = dataset[n].calon.yx[chan];
	writefits_plane(file, data, &hpar);
	for (n=0; n < size; n++)
		for (chan=0; chan<MAX_CHANNELS; chan++)
			data[n*MAX_CHANNELS+chan] = dataset[n].calon.yy[chan];
	writefits_plane(file, data, &hpar);


	writefits_pad_end(file, &hpar);
	fclose(file);


	free(data);
}

/*
 * Sets any datapoints that are nsigma away from the mean to be NAN.
 * Does this iteratively as required.
 */
static void remove_spikes(float data[], int N, float nsigma)
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

static void write_noise_measurements(SpecRecord dataset[], int size, const char* fileroot, int chan)
{
	int n, i;
	char filename[FILENAME_SIZE+1];
	FILE * file;
	float *D[8];

	float *xxon;
	float *xxoff;
	float *xxdiffon;
	float *xxdiffoff;
	float *yyon;
	float *yyoff;
	float *yydiffon;
	float *yydiffoff;

	snprintf(filename, FILENAME_SIZE, "%s_noise.dat", fileroot);
	file = fopen(filename, "w");
	if (file == NULL) {
		printf("ERROR: unable to open '%s' for writing\n", filename);
		return;
	}

	for (i=0; i<8; i++) {
		D[i] = (float *) malloc(sizeof(float) * size);
	}

	//setup the alias names for clarity
	xxon = D[0];
	xxoff = D[1];
	xxdiffon = D[2];
	xxdiffoff = D[3];
	yyon = D[4];
	yyoff = D[5];
	yydiffon = D[6];
	yydiffoff = D[7];

	for (n=0; n<size; n++) {
		xxon[n] = dataset[n].calon.xx[chan];
		xxoff[n] = dataset[n].caloff.xx[chan];
		yyon[n] = dataset[n].calon.yy[chan];
		yyoff[n] = dataset[n].caloff.yy[chan];
	}

	for (n=0; n<size-1; n++)
	{
		xxdiffon[n] = xxon[n] - xxon[n+1];
		xxdiffoff[n] = xxoff[n] - xxoff[n+1];
		yydiffon[n] = yyon[n] - yyon[n+1];
		yydiffoff[n] = yyoff[n] - yyoff[n+1];
	}

	remove_spikes(xxdiffon, size-1, 5);
	remove_spikes(xxdiffoff, size-1, 5);
	remove_spikes(yydiffon, size-1, 5);
	remove_spikes(yydiffoff, size-1, 5);

	fprintf(file, "#xxon xxdiffon xxoff xxdiffoff yyon yydiffon yyoff yydiffoff\n");
	for (i=0; i<8; i++) {
		float sigma, mean;
		mean = compute_mean(D[i], size-1);
		sigma = compute_sigma(D[i], size-1, mean);
		printf("Mean:%6.8f Sigma:%6.8f\n", mean, sigma);
	}

	for (n=0; n<size-1; n++) {
		for (i=0; i<8; i++) {
			float *A = D[i];
			fprintf(file, "%6.8f ", A[n]);
		}
		fprintf(file, "\n");
	}

	for (i=0; i<8; i++) {
		free(D[i]);
	}

	return;
}


static void process_dataset(const char * filepath, int beam, int smooth)
{
	FILE * datafile, *cfgfile;
	int numRecords;
	SpecRecord * dataset;
	ConfigData cfgData;
	float df;
	float fcen;

	char datafilepath[256+1];
	char cfgfilepath[256+1];
	char badfilepath[256+1];
	char filename[80+1];
	char path[80+1];


	/* Open Datafile */
	snprintf(datafilepath, 256, "%s", filepath);
	if ( (datafile = fopen(datafilepath, "rb") ) == NULL )
	{
		printf("ERROR: can't open data file '%s'\n", datafilepath);
		return;
	}


	/* Open Configuration File */
/*	snprintf(cfgfilepath, 256, "%s_cfg", filepath);
	if ( (cfgfile = fopen(cfgfilepath, "r") ) == NULL )
	{
		printf("ERROR: can't open config file '%s'\n", cfgfilepath);
		return;
	}
*/

	/* read config file */
/*	printf("reading config file %s ... \n", cfgfilepath);
	read_cfgfile(cfgfile, &cfgData);
	fclose(cfgfile);
*/
	/* calculate the channel frequencies */
/*	fcen = cfgData.centerMHz;
	df = cfgData.bandwitdhkHz / MAX_CHANNELS / 1000;
	printf ("center frequency: %fMHz\n", fcen);
	printf("channel bandwidth: %fMHz\n", df);
*/
	/* Read Datafile */
	printf("Reading data file %s\n", datafilepath);
	numRecords = read_datafile(datafile, &dataset, beam);
	fclose(datafile);
	printf("Read %i records\n", numRecords);

	if (numRecords <= 0) {
		printf("ERROR: Skipping %s: there are no records\n", datafilepath);
		return;
	}

	if (smooth) {
		printf("Performing frequency smoothing\n");
		perform_freq_smoothing(dataset, numRecords, 0, MAX_CHANNELS);
	}



	get_name_from_filepath(filename, path, FILENAME_SIZE, filepath);
//	printf("Filename: %s\n",filename);

	snprintf(badfilepath, 256, "%s/bad_datapoints.dat", path);
	mark_bad_datapoints(badfilepath, dataset, numRecords);

//	printf("writing raw data to fits ...\n");
//	write_raw_pol_fits(dataset, numRecords, fcen, df, filename);

	printf("Writing rfi data\n");
	write_rfi_data(dataset, numRecords, filename, fcen, df);

	printf("Writing pointing data\n");
	write_pointing(dataset, numRecords, filename);

	printf("Writing band average data\n");
	write_band_average(dataset, numRecords, filename);

	printf("Writing time average data\n");
	write_time_average(dataset, numRecords, filename);

	printf("Writing noise measurements\n");
	write_noise_measurements(dataset, numRecords, filename, 90);

	printf("Done!\n");
	free(dataset);
	return;
}



/* get the beamid out of the filename
 */
static int get_beamid(const char * filename)
{
	regex_t reg;
	regex_t * preg = &reg;
	int errcode;
	const char * regex = "b(.?)"; //TODO: make this an input arg
	regmatch_t pmatch[2];
	size_t nmatch = 2;
	char name[8+1];

	errcode = regcomp(preg, regex, REG_EXTENDED);
	if (errcode != 0) {
		handle_regerr(errcode, preg, stderr);
		return -1;
	}

	errcode = regexec(preg, filename, nmatch, pmatch, 0x0);
	if (errcode != 0) {
		handle_regerr(errcode, preg, stderr);
		return -2;
	}

	regfree(preg);

	//pmatch[0] is the entire match, pmatch[1] is the first subexpression match in '()'
	strnsub(name, filename, pmatch[1].rm_so, pmatch[1].rm_eo, 8);
	return atoi(name);
}


int main(int argc, char *argv[])
{
	int beamid;
	int smooth;
	char * filename;

	if (argc != 3)
	{
		printf("Usage: %s <specfilename> <smooth>\n", argv[0]);
		printf("eg: %s A2174.perpuls_08_00_025346+250905.beam0.53989.spec 1 \n", argv[0]);
		return EXIT_FAILURE;
	}

	filename = argv[1];
	smooth = atoi(argv[2]);
//	beamid = atoi(argv[3]);

	beamid = get_beamid(filename);
	if (beamid < 0) {
		return EXIT_FAILURE;
	}
	printf("Beam:%d\n",beamid);
	process_dataset(filename, beamid, smooth);
	return EXIT_SUCCESS;
}


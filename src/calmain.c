#include "common.h"
#include "cal.h"
#include "rfi.h"
#include "stokes.h"
#include "calibrate.h"
#include "scan.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include "smooth.h"
#include <glob.h>
#include "markdata.h"
#include "programs/fitsio.h"

int multibeam; //SSG

static int read_tcal(const char * filename, float * Tcalx, float * Tcaly)
{
    FILE * file;
    int p;
    char line[80+1];
    char * tok;
const int START_CHAN = 30; //TODO: this is a hack for the start channel number
    
    file = fopen(filename, "r");
    if (file == NULL) {
        printf("ERROR: unable to open file %s\n", filename);
        return -1;
    }
    
    fgets(line, 80, file);
    p = START_CHAN;
    do {
        if (line[0] != '#') 
        {
            tok = strtok(line, " ");
            Tcalx[p] = atof(tok);
            tok = strtok(NULL, " ");
            Tcaly[p] = atof(tok);
            //printf("%f, %f\n", Tcalx[p], Tcaly[p]);
            p++;
        }
        
        fgets(line, 80, file);
    } while (!feof(file));
    
    fclose(file);
    return p-START_CHAN;
}

void write_cal_fits(SpecRecord dataset[], int size, float fcen, float df)
{
	int i, n;
	header_param_list hpar;

	const char * calxxfile = "calxx.fits";
	const char * calyyfile = "calyy.fits";
	const char * calxyfile = "calxy.fits";
	const char * calyxfile = "calyx.fits";

	float *caldataxx  = (float *) malloc (MAX_CHANNELS * size * sizeof(float));
	float *caldataxy  = (float *) malloc (MAX_CHANNELS * size * sizeof(float));
	float *caldatayx  = (float *) malloc (MAX_CHANNELS * size * sizeof(float));
	float *caldatayy  = (float *) malloc (MAX_CHANNELS * size * sizeof(float));

	for (n=0; n<size; n++) {
		for (i=0; i<MAX_CHANNELS; i++) {
			caldataxx[n*MAX_CHANNELS+i] = dataset[n].cal.xx[i];
			caldataxy[n*MAX_CHANNELS+i] = dataset[n].cal.xy[i];
			caldatayx[n*MAX_CHANNELS+i] = dataset[n].cal.yx[i];
			caldatayy[n*MAX_CHANNELS+i] = dataset[n].cal.yy[i];
		}
	}
	init_header_param_list (&hpar);  /* initialize parameter records */
	hpar.bitpix = -32;
	hpar.num_axes = 2;
	hpar.naxis[0] = MAX_CHANNELS;
	hpar.naxis[1] = n;
	sprintf (hpar.ctype[0], "Frequency");
	sprintf (hpar.ctype[1], "Time");
	hpar.crval[0] = fcen;		  /* channels */
	hpar.crval[1] = (dataset[size-1].AST + dataset[0].AST) / 2.0;               /* seconds */
	hpar.crpix[0] = MAX_CHANNELS / 2.0;	/* image center in pixels */
	hpar.crpix[1] = 0.5 + n / 2.0;
	hpar.cdelt[0] = df;                     /* channel */
	hpar.cdelt[1] = 0.2;                     /* seconds */
	sprintf (hpar.bunit, "Kelvin");
	sprintf (hpar.telescope, "Arecibo");

	sprintf (hpar.object, "Calbration values for XX");
	writefits_map (calxxfile, caldataxx, &hpar);
	sprintf (hpar.object, "Calbration values for XY");
	writefits_map (calxyfile, caldataxy, &hpar);
	sprintf (hpar.object, "Calbration values for YX");
	writefits_map (calyxfile, caldatayx, &hpar);
	sprintf (hpar.object, "Calbration values for YY");
	writefits_map (calyyfile, caldatayy, &hpar);
	
	free(caldataxx);
	free(caldataxy);
	free(caldatayx);
	free(caldatayy);
}




static void create_annotations(SpecRecord dataset[], int size)
{
	int n;
	FILE * annfile;
	const char * annfilename = "pointing.ann";
	double h, secderiv;
	int found;

	annfile = fopen(annfilename, "w");
	if (annfile == NULL) {
		printf("ERROR: unable to open '%s' for writing\n", annfilename);
		return;
	}

	fprintf(annfile, "#Annotations\n");

	found = 0;
	fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", dataset[0].RA, dataset[0].DEC, dataset[0].AST);
	for (n=1; n<size-1; n++)
	{
		if (dataset[n].flagBAD) {
			fprintf(annfile, "COLOUR %s\n", "RED"); 
		} else {
			fprintf(annfile, "COLOUR %s\n", "GREEN");
		}
		fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", 
			dataset[n].RA, dataset[n].DEC, dataset[n].AST);

		//write a red circle around pointing rages of change that are too high
		h = (dataset[n+1].AST - dataset[n-1].AST) / 2;
		secderiv = (dataset[n-1].DEC - 2*dataset[n].DEC + dataset[n+1].DEC) / (h*h);
		if (fabs(secderiv) > 0.08) {
			if (!found) {
				found = 1; 
				fprintf(annfile, "COLOUR %s\n", "RED"); 
				fprintf(annfile, "CIRCLE W %7.4f %7.4f %7.4f #%7.2f\n", 
					dataset[n].RA, dataset[n].DEC, 0.025, dataset[n].AST);
			}
		} else {
			found = 0;
		}
	}
	fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", dataset[n].RA, dataset[n].DEC, dataset[n].AST);
	
	fclose(annfile);
}



static void process_dataset(const char * datadirname, const char * datedir, const char * subdir, int beam, float numSigma, float numSigmaThresh, int lowchan, int highchan, int tsmooth, int ignoreRFI, int rfiTolerance, int rfiSpan, int ignoreA_low, int ignoreA_high, float decmin, float decmax, int scan_count_thresh, float Tcalx[], float Tcaly[])
{	
	
	FILE * datafile, *cfgfile;
	int numRecords;
	SpecRecord * dataset;
	ConfigData cfgData;	
	int i;
	float df;
	float freq[MAX_CHANNELS];
	float fcen;
	float maxRA, minRA, maxDEC, minDEC;

	glob_t globbuf;
	char datafilename[128+1];
	char cfgfilename[128+1];
	char badfilename[128+1];
	char globpattern[128+1];

	/* Open Datafile */
	globbuf.gl_offs = 1;
//	sprintf(globpattern, "%s/%s/*.za_scan*.%s.%s.spec", datadirname, datedir, subdir, datedir);
	sprintf(globpattern, "%s/%s/*.perpuls*.%s.%s.spec", datadirname, datedir, subdir, datedir); //ssg for a2174 data
	glob(globpattern, 0, NULL, &globbuf);
	if (globbuf.gl_pathc < 1) {
		printf("ERROR: unable to find file with pattern %s\n", globpattern);
	}
//	else //SSG
//		printf("DIAGNOSTIC: found file with pattern %s:: beam %d\n", globpattern,beam); //SSG

	strncpy(datafilename, globbuf.gl_pathv[0], 128);
	//globfree(&globbuf);	

	if ( (datafile = fopen(datafilename, "rb") ) == NULL )
	{ 
		printf("ERROR: can't open data file '%s'\n", datafilename); 
		return;
	}


	/* Open Configuration File */
	globbuf.gl_offs = 1;
//	sprintf(globpattern, "%s/%s/*.za_scan*.%s.%s.spec_cfg", datadirname, datedir, subdir, datedir);
	sprintf(globpattern, "%s/%s/*.perpuls*.%s.%s.spec_cfg", datadirname, datedir, subdir, datedir);//ssg for a2174 data
	glob(globpattern, 0, NULL, &globbuf);
	if (globbuf.gl_pathc < 1) {
		printf("ERROR: unable to find file with pattern %s\n", globpattern);
	}
	strncpy(cfgfilename, globbuf.gl_pathv[0], 128);
	//globfree(&globfree);	
	if ( (cfgfile = fopen(cfgfilename, "r") ) == NULL )
	{ 
		printf("ERROR: can't open config file '%s'\n", cfgfilename); 
		return;
	}


	/* read config file */
	printf("reading config file %s ... \n", cfgfilename);
	read_cfgfile(cfgfile, &cfgData);
	fclose(cfgfile);

	/* calculate the channel frequencies */
	fcen = cfgData.centerMHz; 
	df = cfgData.bandwitdhkHz / 256 / 1000;
	printf ("center frequency: %fMHz\n", fcen);
	printf("channel bandwidth: %fMHz\n", df);
	for (i=0; i<MAX_CHANNELS; i++) {
		freq[i] = fcen + ((float)(i - 127)) * df;
	}

	/* read data file */
	printf("reading data file %s ... \n", datafilename);
	numRecords = read_datafile(datafile, &dataset, beam);
	fclose(datafile);
	printf("Read %i records\n", numRecords);

	if (numRecords <= 0) {
		printf("ERROR: Skipping %s %s: there are no records!\n", datedir, subdir);
		return;
	}


	/* determine the field size */
	maxRA = 0.0;
	minRA = INFINITY;
	maxDEC = 0.0;
	minDEC = INFINITY;
	for (i=0; i<numRecords; i++) {
		SpecRecord *pRec = &(dataset[i]);
		minRA = (pRec->RA < minRA) ? pRec->RA : minRA;
		maxRA = (pRec->RA > maxRA) ? pRec->RA : maxRA;
		minDEC = (pRec->DEC < minDEC) ? pRec->DEC : minDEC;
		maxDEC = (pRec->DEC > maxDEC) ? pRec->DEC : maxDEC;
	}
	printf("data ranges from RA (%f, %f) DEC (%f, %f)\n", minRA, maxRA, minDEC, maxDEC);
	
	/* mark bad datapoints */
	sprintf(badfilename, "%s/%s/bad_datapoints.dat", datadirname, datedir);
	printf("marking bad data points ...\n");
	mark_bad_datapoints(badfilename, dataset, numRecords);

	/* count remaining good records */
//	for (i=0; i<numRecords; i++) {
//		SpecRecord *pRec = &(dataset[i]);
//	}


	printf("perform frequency smoothing ...\n");
	perform_freq_smoothing(dataset, numRecords);

	if (!ignoreRFI) 
	{
		/* Perform the RFI detection */
		printf("performing aerostat RFI blanking ...\n");
		aerostat_rfi_blanking(dataset, numRecords, lowchan, highchan);
		printf("performing RFI detection ...\n");
		rfi_detection(dataset, numRecords, lowchan, highchan, numSigma, numSigmaThresh, ignoreA_low, ignoreA_high, 0, 0);
		printf("performing RFI spanning ...\n");
		rfi_spanning(dataset, numRecords, lowchan, highchan, rfiSpan);
		printf("performing RFI blanking ...\n");
		rfi_blanking(dataset, numRecords, lowchan, highchan, rfiTolerance);
		printf("performing RFI write ...\n");
		rfi_write(dataset, numRecords, lowchan, highchan, freq);
	}

	/* Write out the annotations */
	printf("creating annotations ...\n");
	create_annotations(dataset, numRecords); 

	/* Compute the cal values */
	printf("compute the raw values of cal ...\n");
//ssg	compute_raw_cal(dataset, numRecords);
	compute_raw_cal(dataset, lowchan, highchan, numRecords);//ssg
	//print_cal(dataset, numRecords, lowchan, highchan, 0);

	/* Compute curve fit cal */
	printf("compute the curve fit cal ...\n");
	linear_fit_cal(dataset, numRecords, lowchan, highchan, ignoreRFI);
	//print_cal(dataset, numRecords, lowchan, highchan, ignoreRFI);

	printf("writing cal fits on curve fit cal ...\n");
	write_cal_fits(dataset, numRecords, fcen, df);

	/* Calculate Stokes parameters */
	printf("calculate stokes parameters ...\n");
	calculate_stokes(dataset, numRecords, lowchan, highchan, ignoreRFI, Tcalx, Tcaly);

	/* Write Channel Data */
	printf("writing channel data...\n");
//	write_channel_data(dataset, numRecords, lowchan, highchan);
//	write_channel_data(dataset, numRecords, 0, 1);
	write_channel_data(dataset, numRecords, lowchan, highchan, beam); //SSG
//	write_channel_data(dataset, numRecords, 0, 1,beam); //SSG

	/* Smooth */
	//printf("smooth the band averaged data in time ...\n");
	//smooth_stokes(dataset, numRecords, lowchan, highchan, tsmooth);

	/* done */
	printf("done!\n");
	free(dataset);
	return;
}

int main(int argc, char *argv[])
{
	int ignoreRFI;
	int tsmooth;
	int lowchan, highchan;
	int rfiTolerance;
	int rfiSpan;
	int beam;
	float numSigmaThresh;
	int ignoreA_low;
	int ignoreA_high;
	int numdirs;
	float decmin, decmax;
	int scan_count_thresh;
	char subdirname[5+1];//SSG
	int beamcounter;//SSG
	float numSigma;
	char * datadirname;
	char ** datedirs;

	float Tcalx[MAX_CHANNELS];
	float Tcaly[MAX_CHANNELS];
	int num_tcal;
	const char *tcalfilename = "Tcal.dat";

	int mjd, i, j;
	char *subdir;

	if (argc != 17) {
		printf("Usage: %s <datadir> <subdir> <beam> <lowchan> <highchan> <ignoreRFI> <RFITolerance> <ignoreA_low> <ignoreA_high> <numsigma> <numsigmathresh> <tsmooth> <rfispan> <decmin> <decmax> <scan_count_thresh>\n", argv[0]);
		printf("eg: %s /n/swift2/galfacts/data/A1863/SPEC beam1 1 3.5 25 230 0 25 0 0 3.5 0.03 5 10\n", argv[0]);
		return EXIT_FAILURE;
	} 
	else 
	{
		datadirname = argv[1];
		subdir = argv[2];
        beam = atoi(argv[3]);
		lowchan = atoi(argv[4]);
		highchan = atoi(argv[5]);
		ignoreRFI = atoi(argv[6]);
		rfiTolerance = atoi(argv[7]);
		ignoreA_low = atoi(argv[8]);
		ignoreA_high = atoi(argv[9]);
		numSigma = atof(argv[10]);
		numSigmaThresh = atof(argv[11]);
		tsmooth = atoi(argv[12]);
		rfiSpan = atoi(argv[13]);
		decmin = atof(argv[14]);
		decmax = atof(argv[15]);
		scan_count_thresh = atoi(argv[16]);
		if(beam == 8) //SSG
			multibeam = 1; //SSG
		else
			multibeam = 0; //SSG
	}

	//Handle the Tcal correction file
	num_tcal = read_tcal(tcalfilename, Tcalx, Tcaly);
	if (num_tcal < 0) {
		printf("ERROR: unable to read %s\n", tcalfilename);
		return EXIT_FAILURE;
	} else if (num_tcal < highchan-lowchan) {
		printf("ERROR: Only read %i factors from %s, but require %i\n", 
			num_tcal, tcalfilename, highchan-lowchan);
		return EXIT_FAILURE;
	} else {
		printf("Read %i factors from %s\n", num_tcal, tcalfilename);
	}

	printf("Channels (%i, %i]\n", lowchan, highchan);
	if (ignoreRFI) {
		printf("Ignoring RFI\n");
	} else {
		printf("RFI Tolerance: %i%%\n", rfiTolerance);
		printf("RFI Spanning: %i\n", rfiSpan);
		printf("Ignore Range: (%i,%i)\n", ignoreA_low, ignoreA_high);
		printf("Num Sigma: %g\n", numSigma);
		printf("Num Sigma Threshold: %g\n", numSigmaThresh);
	}
	printf("Time Smooth: %i\n", tsmooth);

	numdirs = get_date_dirs(datadirname, &datedirs);
	printf("Found %i data directories in %s\n", numdirs, datadirname);
	beamcounter = -1; //SSG
	for (mjd=0; mjd<numdirs; mjd++) 
	{
		const char * datedir = datedirs[mjd];
		mode_t mode = S_IRWXU|S_IRWXG|S_IRWXO;
		mkdir(datedir, mode);
		chdir(datedir);
		//SSG
		printf("Now processing day %s...\n",datedirs[mjd]); //SSG
		if(beam == 8)
		{
			for(j = 0;j < 7;j++) //for each beam
			{
				if(beamcounter == 6)
					beamcounter = 0;
				else
					beamcounter++;
				sprintf(subdirname,"beam%d",beamcounter);
				subdir = subdirname;
				mkdir(subdir, mode);
				chdir(subdir);
				printf("Processing %s\n",subdirname); //SSG
				process_dataset(datadirname, datedir, subdir, beamcounter, numSigma, numSigmaThresh, lowchan, highchan, tsmooth, ignoreRFI, rfiTolerance, rfiSpan, ignoreA_low, ignoreA_high, decmin, decmax, scan_count_thresh, Tcalx, Tcaly); //SSG
				chdir("..");
			}
		}
		else
		{
			beamcounter = beam;
			mkdir(subdir, mode);
			chdir(subdir);
			process_dataset(datadirname, datedir, subdir, beamcounter, numSigma, numSigmaThresh, lowchan, highchan, tsmooth, ignoreRFI, rfiTolerance, rfiSpan, ignoreA_low, ignoreA_high, decmin, decmax, scan_count_thresh, Tcalx, Tcaly); //SSG
			chdir("..");
		}
		//SSG
//		mkdir("beam0", mode);//SSG
//		chdir("beam0");//SSG

//SSG		process_dataset(datadirname, datedir, subdir, beam, numSigma, numSigmaThresh, lowchan, highchan, tsmooth, ignoreRFI, rfiTolerance, rfiSpan, ignoreA_low, ignoreA_high, decmin, decmax, scan_count_thresh, Tcalx, Tcaly);

		chdir("..");
	}

	for (i=0; i<numdirs; i++) {
		free(datedirs[i]);
	}
	free(datedirs);
	return EXIT_SUCCESS;
}


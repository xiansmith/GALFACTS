#include "common.h"
#include "/u/sguram/cfitsio/fitsio.h"
#include "cal.h"
#include "rfi.h"
#include "stokes.h"
#include "calibrate.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include "smooth.h"
#include <glob.h>
#include "markdata.h"
#include "programs/fitsLib.h"


int multibeam; //SSG

static int read_tcal(const char * filename, float * Tcalx, float * Tcaly)
{
    FILE * file;
    int p;
    char line[80+1];
    char * tok;
//	const int START_CHAN = 30; //TODO: this is a hack for the start channel number

    file = fopen(filename, "r");
    if (file == NULL) {
        printf("ERROR: unable to open file %s\n", filename);
        return -1;
    }

    fgets(line, 80, file);
    p = 0;
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
    return p;
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

	//plot the starting dot as yellow
	//TODO: perhaps this should be a circle
	fprintf(annfile, "COLOUR %s\n", "YELLOW");
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

		//write a red circle around pointing ranges of change that are too high
		h = (dataset[n+1].AST - dataset[n-1].AST) / 2;
		secderiv = (dataset[n-1].DEC - 2*dataset[n].DEC + dataset[n+1].DEC) / (h*h);
		if (fabs(secderiv) > 0.08)
		{
			if (!found) { //so we only plot one circle per error
				found = 1;
				fprintf(annfile, "COLOUR %s\n", "RED");
				fprintf(annfile, "CIRCLE W %7.4f %7.4f %7.4f #%7.2f\n",
					dataset[n].RA, dataset[n].DEC, 0.025, dataset[n].AST);
				printf("WARN: DEC pointing error detected\n");
			}

		} else {
			found = 0;
		}

	}
	fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", dataset[n].RA, dataset[n].DEC, dataset[n].AST);

	fclose(annfile);
}

void write_sdfits(SpecRecord dataset[], int size, ConfigData *pCfg,int beam)
{
	fitsfile *fptr;
	int status = 0,true = TRUE,i,j;
	char filename[FLEN_FILENAME];
	double RA,DEC,AST;
	double I,Q,U,V;
	char projid[] = "A2186";
	sprintf(filename,"!%s.%04i%02i%02i.b%d.sdfits",projid,pCfg->year,pCfg->month,pCfg->day,beam);
	fits_create_file(&fptr,filename,&status);
//	fits_write_key(fptr,TLOGICAL,"SIMPLE",&true,"file does conform to FITS standard",&status);
//	temp = -32;
//	fits_write_key(fptr,TINT,"BITPIX",&temp,"data is in single precision floats",&status);
//	temp = 0;
//	fits_write_key(fptr,TINT,"NAXIS",&temp,"number of data axes",&status);
//	fits_write_key(fptr,TLOGICAL,"EXTEND",&true,"FITS dataset may contain extensions",&status);
//	fits_write_comment(fptr,"File written using GALFACTS code by sguram@phas.ucalgary.ca",&status);
	if(status)
	{
		fits_read_errmsg(filename);
		printf("%s\n",filename);
	}

	char *ttype[] = {"TIME","DATA","CRVAL3","CRVAL4","RESTFRQ"};
	char *tform[] = {"1D","2048D","1D","1D","1D"};
	char *tunit[] = {"s","K","deg","deg","Hz"};
	double stokespix =1,stokesval= 1,stokesdelt = 1;
	double freqpix = (MAX_CHANNELS+1)/2.0;
	double freqval= pCfg->centerMHz*1000000,freqdelt = pCfg->bandwitdhkHz*1000/MAX_CHANNELS;
	float bandwidth =pCfg->bandwitdhkHz*1000;
	float exposure = pCfg->integrationTime/1000.0*5,tsys = 40.1;
	char date[11];
	sprintf(date,"%04i-%02i-%02i",pCfg->year,pCfg->month,pCfg->day);
	int tfields = 5;
	long int naxes[] = {2048,1,1,1};
	double restfrq = 1420405752.00;
	fits_create_tbl(fptr,BINARY_TBL,0,tfields,ttype,tform,tunit,"SINGLE DISH",&status);
//	fits_write_key(fptr,TSTRING,"EXTNAME","SINGLE DISH","name of this binary table extension",&status);
	fits_write_key(fptr,TSTRING,"PROJID",projid,"project id",&status);
	fits_write_key(fptr,TSTRING,"TELESCOP","Arecibo 305m","telescope name",&status);
	fits_write_key(fptr,TINT,"BEAM",&beam,"ALFA beam",&status);
	fits_write_key(fptr,TSTRING,"DATE-OBS",date,"observation date",&status);
//	fits_write_key(fptr,TSTRING,"TTYPE1","TIME","AST",&status);
//	fits_write_key(fptr,TINT,"TFORM1","1D","field name",&status);
//	fits_write_key(fptr,TINT,"TUNIT1","s","field name",&status);
	fits_write_key(fptr,TFLOAT,"EXPOSURE",&exposure,"integration time per point in seconds",&status);
	printf("Bandwidth: %f\n",bandwidth);
	fits_write_key(fptr,TFLOAT,"BANDWID",&bandwidth,"spectral bandwidth",&status);
	fits_write_key(fptr,TFLOAT,"TSYS",&tsys,"nominal tsys in Kelvin",&status);
	fits_write_key(fptr,TSTRING,"SPECSYS","TOPOCENT","Doppler reference frame of observation",&status);
//	fits_write_key(fptr,TSTRING,"TTYPE2","DATA","field name",&status);
//	fits_write_key(fptr,TINT,"TFORM2","8192D","field name",&status);
//	fits_write_key(fptr,TINT,"TUNIT2","K","field name",&status);
//	fits_write_key(fptr,TSTRING,"TDIM4","(4,2048)","field name",&status);
	fits_write_tdim(fptr,2,4,naxes,&status);
	fits_write_key(fptr,TSTRING,"CTYPE1","FREQ","field name",&status);
	fits_write_key(fptr,TDOUBLE,"CRPIX1",&freqpix,"field name",&status);
	fits_write_key(fptr,TDOUBLE,"CRVAL1",&freqval,"field name",&status);
	fits_write_key(fptr,TDOUBLE,"CDELT1",&freqdelt,"field name",&status);
	fits_write_key(fptr,TSTRING,"CUNIT1","Hz","field name",&status);
	fits_write_key(fptr,TSTRING,"CTYPE2","STOKES","field name",&status);
	fits_write_key(fptr,TDOUBLE,"CRPIX2",&stokespix,"field name",&status);
	fits_write_key(fptr,TDOUBLE,"CRVAL2",&stokesval,"field name",&status);
	fits_write_key(fptr,TDOUBLE,"CDELT2",&stokesdelt,"field name",&status);
	fits_write_key(fptr,TSTRING,"CTYPE3","RA","field name",&status);
	fits_write_key(fptr,TDOUBLE,"CRPIX3",&stokespix,"field name",&status);
//	fits_write_key(fptr,TINT,"TTYPE3","CRVAL3","field name",&status);
//	fits_write_key(fptr,TINT,"TFORM3","1D","field name",&status);
//	fits_write_key(fptr,TINT,"TUNIT3","deg","field name",&status);
	fits_write_key(fptr,TDOUBLE,"CDELT3",&stokesdelt,"field name",&status);
	fits_write_key(fptr,TSTRING,"CTYPE4","DEC","field name",&status);
	fits_write_key(fptr,TDOUBLE,"CRPIX4",&stokespix,"field name",&status);
//	fits_write_key(fptr,TINT,"TTYPE4","CRVAL4","field name",&status);
//	fits_write_key(fptr,TINT,"TFORM4","1D","field name",&status);
//	fits_write_key(fptr,TINT,"TUNIT4","deg","field name",&status);
	fits_write_key(fptr,TDOUBLE,"CDELT4",&stokesdelt,"field name",&status);

	fits_write_key(fptr,TSTRING,"OBJECT","Field D","field name",&status);

	size = size - size%5;

	for (i=0; i<size;)
	{
		int i1 = i;
		SpecRecord * pRec1 = &(dataset[i]);
		i++;
      	switch(beam)
            	{
				case 0:
					break;
				case 1:
		        	pRec1->RA -= 0.05422437;
				    pRec1->DEC -= 0.01237836;
					break;
				case 2:
	        		pRec1->RA -= 0.01685216;
			        pRec1->DEC -= 0.05022825;
			        break;
				case 3:
	        		pRec1->RA += 0.03753338;
			        pRec1->DEC -= 0.03779388;
					break;
				case 4:
	        		pRec1->RA += 0.05388680;
			        pRec1->DEC += 0.01303928;
					break;
				case 5:
	        		pRec1->RA += 0.01623990;
			        pRec1->DEC += 0.05000566;
					break;
				case 6:
	        		pRec1->RA -= 0.03743130;
			        pRec1->DEC += 0.03702619;
					break;
				default:
					printf("WARN: invalid beam (%i) specified.\n", beam);
	    }

		double nulval;
		int anynul;
		//average to 1 second
		SpecRecord * pRec2 = &(dataset[i]);
		i++;
		SpecRecord * pRec3 = &(dataset[i]);
		i++;
		SpecRecord * pRec4 = &(dataset[i]);
		i++;
		SpecRecord * pRec5 = &(dataset[i]);
		i++;

		RA = pRec1->RA;
		DEC = pRec1->DEC;
		AST = pRec1->AST;
		fits_write_col(fptr,TDOUBLE,1,i1+1,1,1,&AST,&status);

		for(j = 0;j < MAX_CHANNELS;j++)
		{
			if(isnan(pRec1->stokes.I[j]))
				pRec1->stokes.I[j] = 0.0;
//				printf("NAN in Record: %d Channel: %d\n",i1,j);
			if(isnan(pRec2->stokes.I[j]))
				pRec2->stokes.I[j] = 0.0;
//				printf("NAN in Record: %d Channel: %d\n",i1+1,j);
			if(isnan(pRec3->stokes.I[j]))
				pRec3->stokes.I[j] = 0.0;
//				printf("NAN in Record: %d Channel: %d\n",i1+2,j);
			if(isnan(pRec4->stokes.I[j]))
				pRec4->stokes.I[j] = 0.0;
//				printf("NAN in Record: %d Channel: %d\n",i1+3,j);
			if(isnan(pRec5->stokes.I[j]))
				pRec5->stokes.I[j] = 0.0;
//				printf("NAN in Record: %d Channel: %d\n",i1+4,j);


			I = (double)(pRec1->stokes.I[j])+(double)(pRec2->stokes.I[j])+(double)(pRec3->stokes.I[j])\
					+(double)(pRec4->stokes.I[j])+(double)(pRec5->stokes.I[j]);
			I = I/5.0;
//			if(j == 100)
	//			printf("%lf %lf %lf %lf %lf\n",pRec1->stokes.I[j],pRec2->stokes.I[j],pRec3->stokes.I[j],pRec4->stokes.I[j],pRec5->stokes.I[j]);
			fits_write_col(fptr,TDOUBLE,2,i1+1,j+1,1,&I,&status);
		//	fits_write_col(fptr,TDOUBLE,2,i+1,j+1+2048,1,&Q,&status);
		//	fits_write_col(fptr,TDOUBLE,2,i+1,j+1+4096,1,&U,&status);
		//	fits_write_col(fptr,TDOUBLE,2,i+1,j+1+6144,1,&V,&status);
	//		if(j == 100)
	//		printf("Writing %lf %lf %lf %lf\n",RA,DEC,restfrq,I);
//			if(j == 100)
//			{
//				fits_read_col(fptr,TDOUBLE,2,i1+1,j+1,1,&nulval,&I,&anynul,&status);
		//		fits_read_col(fptr,TDOUBLE,2,i+1,j+1+2048,1,&nulval,&Q,&anynul,&status);
		//		fits_read_col(fptr,TDOUBLE,2,i+1,j+1+4096,1,&nulval,&U,&anynul,&status);
		//		fits_read_col(fptr,TDOUBLE,2,i+1,j+1+6144,1,&nulval,&V,&anynul,&status);
		//		printf("Reading %lf\n",I);
//			}
		}

		fits_write_col(fptr,TDOUBLE,3,i1+1,1,1,&RA,&status);
		fits_write_col(fptr,TDOUBLE,4,i1+1,1,1,&DEC,&status);
		fits_write_col(fptr,TDOUBLE,5,i1+1,1,1,&restfrq,&status);
	//	fits_read_col(fptr,TDOUBLE,2,i1+1,100+1,1,&nulval,&I,&anynul,&status);
	//	fits_read_col(fptr,TDOUBLE,3,i1+1,1,1,&nulval,&RA,&anynul,&status);
	//	fits_read_col(fptr,TDOUBLE,4,i1+1,1,1,&nulval,&DEC,&anynul,&status);
	//	fits_read_col(fptr,TDOUBLE,5,i1+1,1,1,&nulval,&restfrq,&anynul,&status);
	//	printf("Reading %lf %lf %lf %lf\n",RA,DEC,restfrq,I);
	} //end for timestep
	if(status)
	{
		fits_read_errmsg(filename);
		printf("%s\n",filename);
	}
	fits_close_file(fptr,&status);
}

static void process_dataset(const char * datadirname, const char * datedir, const char * subdir, int beam, int band, float numSigma, float numSigmaThresh, int lowchan, int highchan, int tsmooth, int ignoreRFI, int rfiTolerance, int rfiSpan, int ignoreA_low, int ignoreA_high, float Tcalx[], float Tcaly[])
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

	//fix to allow for both precursor and main run files to be run
	if(band == -1)
		sprintf(globpattern, "%s/%s/*.*.beam%i.*.spec", datadirname, datedir, beam);
	else
		sprintf(globpattern, "%s/%s/*.*.b*%is*%i.*.spec", datadirname, datedir, beam, band);

	glob(globpattern, 0, NULL, &globbuf);
	if (globbuf.gl_pathc < 1) {
		printf("ERROR: unable to find file with pattern %s\n", globpattern);
	}
	else //SSG
		printf("DIAGNOSTIC: found file %s \n", globbuf.gl_pathv[0]); //SSG

	strncpy(datafilename, globbuf.gl_pathv[0], 128);
	//globfree(&globbuf);

	if ( (datafile = fopen(datafilename, "rb") ) == NULL )
	{
		printf("ERROR: can't open data file '%s'\n", datafilename);
		return;
	}


	/* Open Configuration File */
	globbuf.gl_offs = 1;

	//fix band for ZOA
	band = -1;
	//fix to allow for both precursor and main run files to be run
	if(band == -1)
		sprintf(globpattern, "%s/%s/*.*.beam%i.*.spec_cfg", datadirname, datedir, beam);
	else
		sprintf(globpattern, "%s/%s/*.*.b*%is*%i.*.spec_cfg", datadirname, datedir, beam, band);

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
	printf("Reading config file %s \n", cfgfilename);
	read_cfgfile(cfgfile, &cfgData);
	fclose(cfgfile);

	/* calculate the channel frequencies */
	fcen = cfgData.centerMHz;
	df = cfgData.bandwitdhkHz / MAX_CHANNELS / 1000;
	printf("Center frequency: %fMHz\n", fcen);
	printf("Channel bandwidth: %fMHz\n", df);
	for (i=0; i<MAX_CHANNELS; i++) {
		freq[i] = fcen + ((float)(i - MAX_CHANNELS/2 - 1)) * df;
	}

	/* read data file */
	printf("Reading data file %s. \n", datafilename);
	numRecords = read_datafile(datafile, &dataset, beam);
	fclose(datafile);
	printf("Read %i records\n", numRecords);

	if (numRecords <= 0) {
		printf("ERROR: Skipping %s %s: there are no records!\n", datedir, subdir);
		return;
	}

	numRecords = 100;
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
	printf("Data ranges from RA (%f, %f) DEC (%f, %f)\n", minRA, maxRA, minDEC, maxDEC);

//	printf("Performing frequency smoothing\n");
//	perform_freq_smoothing(dataset, numRecords, lowchan, highchan);

/*	if (!ignoreRFI)
	{
*/		/* Perform the RFI detection */
/*		printf("Performing aerostat RFI blanking\n");
		aerostat_rfi_blanking(dataset, numRecords, lowchan, highchan);
		printf("Performing RFI detection\n");
		rfi_detection(dataset, numRecords, lowchan, highchan, numSigma, numSigmaThresh, ignoreA_low, ignoreA_high, 0, 0);
		printf("Performing RFI spanning\n");
		rfi_spanning(dataset, numRecords, lowchan, highchan, rfiSpan);
		printf("Performing RFI blanking\n");
		rfi_blanking(dataset, numRecords, lowchan, highchan, rfiTolerance);
		printf("Performing RFI write\n");
		rfi_write(dataset, numRecords, lowchan, highchan, freq);
	}
*/
	/* Compute the cal values */
	printf("Compute the raw values of cal\n");
	compute_raw_cal(dataset, numRecords, lowchan, highchan);

	/* Compute curve fit cal */
//	printf("Compute smoothed cal\n");
//	smooth_cal_bandaverage(dataset, numRecords, lowchan, highchan, 120 * 5 - 1, 2.5);
        printf("Compute linear cal\n");
        linear_fit_cal(dataset, numRecords, lowchan, highchan, ignoreRFI);

	/* Calculate Stokes parameters */
	printf("Calculating stokes parameters\n");
	calculate_stokes(dataset, numRecords, lowchan, highchan, ignoreRFI, Tcalx, Tcaly);

	printf("Writing SDFITS\n");
	write_sdfits(dataset,numRecords,&cfgData,beam);

	printf("Done!\n");
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
	int band;
	float numSigmaThresh;
	int ignoreA_low;
	int ignoreA_high;
	int numdirs;
//	float decmin, decmax;
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

	if (argc != 15) {
		printf("Usage: %s <datadir> <subdir> <beam> <band> <lowchan> <highchan> <ignoreRFI> \n\
		 <RFITolerance> <ignoreA_low> <ignoreA_high> <numsigma> <numsigmathresh> <tsmooth>  \n\
		<rfispan> \n", argv[0]);
//		printf("eg: %s /n/swift2/galfacts/data/A1863/SPEC beam1 1 3.5 25 230 0 25 0 0 3.5 0.03 5 10\n", argv[0]);
		return EXIT_FAILURE;
	}
	else
	{
		datadirname = argv[1];
		subdir = argv[2];
	        beam = atoi(argv[3]);
	        band = atoi(argv[4]);
		lowchan = atoi(argv[5]);
		highchan = atoi(argv[6]);
		ignoreRFI = atoi(argv[7]);
		rfiTolerance = atoi(argv[8]);
		ignoreA_low = atoi(argv[9]);
		ignoreA_high = atoi(argv[10]);
		numSigma = atof(argv[11]);
		numSigmaThresh = atof(argv[12]);
		tsmooth = atoi(argv[13]);
		rfiSpan = atoi(argv[14]);
//		decmin = atof(argv[15]);
//		decmax = atof(argv[16]);
//		scan_count_thresh = atoi(argv[17]);
		if(beam == MULTIBEAM) //SSG
			multibeam = 1; //SSG
		else
			multibeam = 0; //SSG
	}
	//fix ignore RFI for ZOA
	ignoreRFI = 1;
	//Handle the Tcal correction file
	num_tcal = read_tcal(tcalfilename, Tcalx, Tcaly);
	if (num_tcal < 0)
	{
		printf("ERROR: unable to read %s\n", tcalfilename);
		printf("Setting all Tcal parameters to 1.00\n");
		for(i = 0;i < MAX_CHANNELS;i++)
		{
			Tcalx[i] = 1.0;
			Tcaly[i] = 1.0;
		}
//		return EXIT_FAILURE;
	}
	else if (num_tcal < MAX_CHANNELS)
	{
		printf("ERROR: Only read %i factors from %s, but require %i\n",
		num_tcal, tcalfilename, MAX_CHANNELS);
		printf("Setting all Tcal parameters to 1.00\n");
		for(i = 0;i < MAX_CHANNELS;i++)
		{
			Tcalx[i] = 1.0;
			Tcaly[i] = 1.0;
		}
//		return EXIT_FAILURE;
	}
	else
	{
		printf("Read %i factors from %s\n", num_tcal, tcalfilename);
	}

	printf("Channels (%i, %i]\n", lowchan, highchan);

	if (ignoreRFI)
	{
		printf("Ignoring RFI\n");
	}
	else
	{
		printf("RFI Tolerance: %i%%\n", rfiTolerance);
		printf("RFI Spanning: %i\n", rfiSpan);
		printf("Ignore Range: (%i,%i)\n", ignoreA_low, ignoreA_high);
		printf("Num Sigma: %g\n", numSigma);
		printf("Num Sigma Threshold: %g\n", numSigmaThresh);
	}
	printf("Time Smooth: %i\n", tsmooth);

	char filename[40+1];
	sprintf(filename,"Days.list");
	FILE * dayfile;
	dayfile = fopen(filename,"r");
	if(dayfile != NULL)
	{
		printf("Found file: \"Days.list\"\n");
		numdirs = jsd_line_count(dayfile);
		datedirs = (char **) malloc(sizeof(char*) * numdirs);
		char tempstr[15];
		for(i=0;i<numdirs;i++)
		{
			fscanf(dayfile,"%s",tempstr);
			datedirs[i] = (char *) malloc(sizeof(char) * strlen(tempstr));
			strcpy(datedirs[i],tempstr);
		}
		fclose(dayfile);
	}
	else
	{
		numdirs = get_date_dirs(datadirname, &datedirs);
	}

	printf("Found %i data directories in %s\n", numdirs, datadirname);
	printf("The following days will be processed:\n");
	for(i=0;i<numdirs;i++)
		printf("%s,",datedirs[i]);
	printf("\n");
	beamcounter = -1; //SSG
	for (mjd=0; mjd<numdirs; mjd++)
	{
		const char * datedir = datedirs[mjd];
		mode_t mode = S_IRWXU|S_IRWXG|S_IRWXO;
		mkdir(datedir, mode);
		chdir(datedir);
		//SSG
		printf("Now processing day %s...\n",datedirs[mjd]); //SSG
		if(beam == MULTIBEAM)
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
				process_dataset(datadirname, datedir, subdir, beamcounter, band, numSigma, numSigmaThresh, lowchan, highchan, tsmooth, ignoreRFI, rfiTolerance, rfiSpan, ignoreA_low, ignoreA_high, Tcalx, Tcaly); //SSG
				chdir("..");
			}
		}
		else
		{
			beamcounter = beam;
			mkdir(subdir, mode);
			chdir(subdir);
			process_dataset(datadirname, datedir, subdir, beamcounter, band, numSigma, numSigmaThresh, lowchan, highchan, tsmooth, ignoreRFI, rfiTolerance, rfiSpan, ignoreA_low, ignoreA_high, Tcalx, Tcaly); //SSG
			chdir("..");
		}

		chdir("..");
	}

	for (i=0; i<numdirs; i++) {
		free(datedirs[i]);
	}
	free(datedirs);
	return EXIT_SUCCESS;
}


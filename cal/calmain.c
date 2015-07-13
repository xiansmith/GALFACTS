#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <glob.h>

#include "common.h"
#include "cal.h"
#include "rfi.h"
#include "stokes.h"
#include "calibrate.h"
#include "smooth.h"
#include "markdata.h"
#include "jsd_futil.h"
#include "chebyshev.h"
#include "tcal_on_fly.h"
//--------------------------------------------------------------------------------------------------------
int multibeam;

clock_t sclock;
double sec_start;
//--------------------------------------------------------------------------------------------------------
static void start_clock(void)
{
	struct timeval time;
	gettimeofday(&time, NULL);
	sec_start = (double)time.tv_sec + (double)time.tv_usec * .000001;
	sclock = clock();
}
//--------------------------------------------------------------------------------------------------------
static void read_clock(void)
{
	printf("Computation time = %f\n", ((double)clock() - sclock)/CLOCKS_PER_SEC);

	struct timeval time;
	gettimeofday(&time, NULL);
	double sec_end = (double)time.tv_sec + (double)time.tv_usec * .000001;
	printf("Wall time = %f\n", sec_end - sec_start);
}
//--------------------------------------------------------------------------------------------------------
static void read_tcal(const char *filename, float *Tcalx, float *Tcaly) 
{
    FILE * file;
    int p, i;
    char line[80+1];
    char * tok;

    file = fopen(filename, "r");
    if(file == NULL) 
    	{
        printf("ERROR: unable to open file %s\n", filename);
		printf("Setting all Tcal parameters to 1.00\n");
 		for(i = 0; i < MAX_CHANNELS; i++)
			{
			Tcalx[i] = 1.0;
			Tcaly[i] = 1.0;
			}
    	}
	else
		{
   		fgets(line, 80, file);
   		p = 0;
   		do 
   			{
       		if (line[0] != '#')
       			{
           		tok = strtok(line, " ");
           		Tcalx[p] = atof(tok);
           		tok = strtok(NULL, " ");
           		Tcaly[p] = atof(tok);
           		p++;
       			}
       		fgets(line, 80, file);
   			}
   			while (!feof(file));
   		fclose(file);
		if (p < MAX_CHANNELS)
			{
			printf("ERROR: Only read %i factors from %s, but require %i\n", p, filename, MAX_CHANNELS);
			printf("Setting all Tcal parameters to 1.00\n");
			for(i = 0; i < MAX_CHANNELS; i++)
				{
				Tcalx[i] = 1.0;
				Tcaly[i] = 1.0;
				}
			}
		else
			{
			printf("Read %i factors from %s\n", p, filename);
			}
		} 		 
}
//--------------------------------------------------------------------------------------------------------
static void create_annotations(SpecRecord dataset[], int size)
{
	int n;
	FILE * annfile;

	const char * annfilename = "pointing.ann";
	double h, secderiv;
	int found;

	annfile = fopen(annfilename, "w");
	if(annfile == NULL) 
		{
		printf("ERROR: unable to open '%s' for writing\n", annfilename);
		return;
		}

	fprintf(annfile, "#Annotations\n");

	fprintf(annfile, "COLOUR %s\n", "YELLOW");
	fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", dataset[0].RA, dataset[0].DEC, dataset[0].AST);
	found = 0;
	for(n=1; n < size-1; n++)
		{
		if (dataset[n].flagBAD) fprintf(annfile, "COLOUR %s\n", "RED"); else fprintf(annfile, "COLOUR %s\n", "GREEN");
		fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", dataset[n].RA, dataset[n].DEC, dataset[n].AST);
		h = (dataset[n+1].AST - dataset[n-1].AST) / 2;
		secderiv = (dataset[n-1].DEC - 2*dataset[n].DEC + dataset[n+1].DEC) / (h*h);
		if(fabs(secderiv) > 0.08)
			{
			if(!found) 
				{
				found = 1;
				fprintf(annfile, "COLOUR %s\n", "RED");
				fprintf(annfile, "CIRCLE W %7.4f %7.4f %7.4f #%7.2f\n", dataset[n].RA, dataset[n].DEC, 0.025, dataset[n].AST);
				}
			} 
			else found = 0;
		}

	fclose(annfile);
}
//--------------------------------------------------------------------------------------------------------
static void process_dataset(const char *field, const char *datadirname, const char *datedir, const char *subdir, int beam, int band, int lowchan, int highchan, \
int RFIF, int RFIT, float numSigmaF, float numSigmaT, int freqSmoothing, \
int uvDenoising, float uvDenoisingTau, float uvDenoisingLambda, float hidrogenfreq, float hidrogenband, int calskyfiles, int annfiles, int fit_smooth, int twindow, int fwindow, int *badchannels, float RAmin,float RAmax,float DECmin,float DECmax)
{
	FILE * datafile, *cfgfile;
	glob_t globbuf;
	char globpattern[128+1];
	ConfigData cfgData;
	float fcen, df, freq[MAX_CHANNELS];
	int numRecords, k;
	SpecRecord * dataset;



	start_clock();
	printf("Locating/opening files\n");
	
	globbuf.gl_offs = 1;
	if(band == -1) sprintf(globpattern, "%s/%s/*.*.beam%i.*.spec", datadirname, datedir, beam); 
	else sprintf(globpattern, "%s/%s/*.*.b*%is*%i.*.spec", datadirname, datedir, beam, band); 
	glob(globpattern, 0, NULL, &globbuf);
	if(globbuf.gl_pathc < 1) printf("ERROR: unable to find file with pattern %s\n", globpattern); 
	else printf("Found file: %s \n", globbuf.gl_pathv[0]); 
	char datafilename[128+1]; strncpy(datafilename, globbuf.gl_pathv[0], 128);
	if((datafile = fopen(datafilename, "rb")) == NULL){printf("ERROR: can't open data file '%s'\n", datafilename); return;}

	globbuf.gl_offs = 1;
	if(band == -1) sprintf(globpattern, "%s/%s/*.*.beam%i.*.spec_cfg", datadirname, datedir, beam); 
	else sprintf(globpattern, "%s/%s/*.*.b*%is*%i.*.spec_cfg", datadirname, datedir, beam, band);
	glob(globpattern, 0, NULL, &globbuf);
	if (globbuf.gl_pathc < 1) printf("ERROR: unable to find file with pattern %s\n", globpattern); 
	else printf("Found file: %s \n", globbuf.gl_pathv[0]); 
	char cfgfilename[128+1]; strncpy(cfgfilename, globbuf.gl_pathv[0], 128);
	if((cfgfile = fopen(cfgfilename, "r")) == NULL){printf("ERROR: can't open config file '%s'\n", cfgfilename); return;}
	
	read_clock(); start_clock();
	printf("Reading config file: %s \n", cfgfilename);
	read_cfgfile(cfgfile, &cfgData);
	fclose(cfgfile);
	fcen = cfgData.centerMHz;
	df = cfgData.bandwitdhkHz / MAX_CHANNELS / 1000;
	printf("Center frequency: %fMHz\n", fcen);
	printf("Channel bandwidth: %fMHz\n", df);	
	printf("Hidrogen frequency: %fMHz\n", hidrogenfreq);
	printf("Hidrogen bandwidth: %fMHz, %fMHz\n", hidrogenfreq-hidrogenband, hidrogenfreq+hidrogenband);
	
	for(k=0; k<MAX_CHANNELS; k++) freq[k] = fcen - ((float)(k + 1 - MAX_CHANNELS/2)) * df; // !!!!!!!!!

	read_clock(); start_clock();
	printf("Writing freq.dat\n");
	FILE *freq_file;
	freq_file=fopen("freq.dat", "w");
	for(k=0; k<MAX_CHANNELS; k++) fprintf(freq_file, "%d\t%f\n", k, freq[k]);
	fclose(freq_file);
		
	printf("Total bandwidth: %fMHz, %fMHz\n", freq[0], freq[MAX_CHANNELS-1]);
	
	read_clock(); start_clock();
	printf("Reading data file: %s. \n", datafilename);
	numRecords = read_datafile(datafile, &dataset, beam, RAmin, RAmax, DECmin, DECmax);
	fclose(datafile);
	printf("Read %i records\n", numRecords);
	if (numRecords <= 0){printf("ERROR: Skipping %s %s: there are no records!\n", datedir, subdir); return;}
	read_clock(); start_clock();

	printf("Marking bad data points\n");
	char badfilename[128+1]; 
	sprintf(badfilename, "%s/%s/bad_datapoints.dat", datadirname, datedir);
	mark_bad_datapoints(badfilename, dataset, numRecords);
	read_clock(); start_clock();
	
	if(freqSmoothing)
		{
		printf("Performing frequency smoothing\n");  
		perform_freq_smoothing(dataset, numRecords, lowchan, highchan);
		read_clock(); start_clock();
		}
	
	if(RFIF)
		{
		printf("Performing RFI detection in frequency domain\n");
		rfi_detection_frequency_domain(dataset, numRecords, lowchan, highchan, numSigmaF, hidrogenfreq, hidrogenband, freq);
		}
	
	read_clock();start_clock();
	
	switch( RFIT )
		{
		case 1:	
					printf("Performing RFI detection in time domain 1\n");
					rfi_detection_time_domain1(field, dataset, numRecords, lowchan, highchan, numSigmaT, hidrogenfreq, hidrogenband, freq); 
					break;
		case 2 :
					printf("Performing RFI detection in time domain 2\n");
					rfi_detection_time_domain2(field, dataset, numRecords, lowchan, highchan, numSigmaT, hidrogenfreq, hidrogenband, freq); 
					break;
		case 3 : 	
					printf("Performing RFI detection in time domain 3\n");
					rfi_detection_time_domain3(field, dataset, numRecords, lowchan, highchan, numSigmaT, hidrogenfreq, hidrogenband, freq); 			
					break;
		default  : 	printf( "No RFI detection in time domain!\n" );
					break;
		}		

	read_clock(); start_clock();
	printf("Marking bad channels\n");
	mark_bad_channels(dataset, numRecords, lowchan, highchan, numSigmaT, hidrogenfreq, hidrogenband, freq, badchannels);
	read_clock(); start_clock();
	if(annfiles)
	{
		if( RFIF || RFIT )
		{
			printf("Create out of band RFI annotation files\n");
			outofbandrfi_ann(dataset, numRecords, lowchan, highchan);
			rfi_ann(dataset, numRecords, lowchan, highchan, freq);
			read_clock(); start_clock();
		}
		printf("Creating pointing annotation file\n"); 
		create_annotations(dataset, numRecords);
		read_clock(); start_clock();
	}
	
	printf("Writing RFI plot data\n");
	write_rfi_data( dataset, numRecords, lowchan, highchan );
	read_clock(); start_clock();
	
	printf("Compute the raw values of cal\n");
	compute_raw_cal(dataset, numRecords, lowchan, highchan);
	read_clock(); start_clock();

	printf("Initialize Tcal arrays\n");
	float Tcalx[MAX_CHANNELS];
	float Tcaly[MAX_CHANNELS];

	float Tcalx_s[MAX_CHANNELS];
	float Tcaly_s[MAX_CHANNELS];

        int i;
        for(i = 0;i < MAX_CHANNELS;i++)
        {
                Tcalx[i] = 1.0;
                Tcaly[i] = 1.0;
                Tcalx_s[i] = 1.0;
                Tcaly_s[i] = 1.0;
        }
	read_clock(); start_clock();

	if(fit_smooth == 1)
	{	
		printf("Compute linear cal\n"); 
		linear_fit_cal(dataset, numRecords, 0, MAX_CHANNELS, RFIF, fwindow);	
	}
	else if(fit_smooth == 2)
	{	
		printf("Compute simple smooth cal\n");

		//Change according to some days to deal with  jumps
		simple_smooth_cal(dataset, 0, numRecords ,numRecords, lowchan, highchan, twindow,fwindow);
	}
	else if(fit_smooth == 0)
	{
		printf("Compute rolling smooth cal\n"); 
		rolling_smooth_cal(dataset, numRecords, lowchan, highchan, twindow,fwindow);
	}
	read_clock(); start_clock();
	printf("Calculating Stokes parameters\n");
	int r;
	
	if(fit_smooth)
	{	
		//compute_tcal(dataset, numRecords,  0, highchan, hidrogenfreq, hidrogenband,freq, badchannels, Tcalx, Tcaly, 0, numRecords);
		//norm_one_tcal(0,MAX_CHANNELS,badchannels,Tcalx,Tcaly);
		//read_clock(); start_clock();
		calculate_stokes(dataset, numRecords, 0, MAX_CHANNELS, RFIF, calskyfiles, Tcalx, Tcaly, uvDenoising, uvDenoisingTau, uvDenoisingLambda,0,numRecords);
		//write_tcal(Tcalx,Tcaly,numRecords,0,MAX_CHANNELS);
	}
	else
	{
/*		for(r = 0;r < numRecords;r++ )
		{

		        compute_tcal(dataset, numRecords,  lowchan, highchan, hidrogenfreq, hidrogenband,freq, badchannels, Tcalx, Tcaly, r, twindow);
		        for(k=0;k<MAX_CHANNELS;k++)
		        {
		        	Tcalx_s[k] = Tcalx[k];
		        	Tcalx_s[k] = Tcalx[k];
	    		}

			norm_one_tcal(0,MAX_CHANNELS,badchannels,Tcalx_s,Tcaly_s);
			diffusion_filter(Tcalx_s, MAX_CHANNELS, fwindow);
			diffusion_filter(Tcaly_s, MAX_CHANNELS, fwindow);

		        if(!(r%10000))
        		        write_tcal(Tcalx_s,Tcaly_s,r,0,MAX_CHANNELS);

			calculate_stokes(dataset, numRecords, lowchan, highchan, RFIF, calskyfiles, Tcalx_s, Tcaly_s, uvDenoising, uvDenoisingTau, uvDenoisingLambda,r,r+1);

		}
*/
		calculate_stokes(dataset, numRecords, lowchan, highchan, RFIF, calskyfiles, Tcalx_s, Tcaly_s, uvDenoising, uvDenoisingTau, uvDenoisingLambda,0,numRecords);
	}

	read_clock(); start_clock();
	printf("Writing channel data to single file\n");
	write_binary_channel_data_single_file(dataset, numRecords, lowchan, highchan);
	//write_channel_data(dataset, numRecords, lowchan, highchan);
	// old method: write_binary_channel_data(dataset, numRecords, lowchan, highchan);

	read_clock(); start_clock();
	printf("Averaging data\n");
	average_stokes(dataset, numRecords, lowchan, highchan, hidrogenfreq, hidrogenband, freq);
	read_clock(); start_clock();
	
	printf("Writing binary channel data\n");
	write_binary_channel_data(dataset, numRecords, 0, 1);
	read_clock(); start_clock();

	// new method, new filename but old format
	// average data is still in channel 0, but no longer overwrite the real channel 0
	printf("Writing channel data\n");
	write_channel_data(dataset, numRecords, 0, 1);

	read_clock();
	free(dataset);
	return;
}
//--------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	int i, j;
	clock_t time0 = clock();

	struct timeval time;
	gettimeofday(&time, NULL);
	double sec0 = (double)time.tv_sec + (double)time.tv_usec * .000001;


	if ((argc < 26) || (argc > 27))
		{
		printf("Usage: %s <parameters_list>\n", argv[0]);
		printf("Got %d arguments\n", argc );
		return EXIT_FAILURE;
		}
	printf("***************************************************\n"); 
	printf("***************************************************\n");
	printf("**                                               **\n");
	printf("**   CALMAIN, GALFACTS calibration program       **\n");
	printf("**                                               **\n");
	printf("***************************************************\n");
	printf("***************************************************\n");
	char *field; field = argv[1]; printf("field = %s\n", field);
	char *datadirname; datadirname = argv[2]; printf("datadirname = %s\n", datadirname);
	char *beam; beam = argv[3]; printf("beam = %s\n", beam);
	int band = atoi(argv[4]); printf("band = %d\n", band);
	int lowchan = atoi(argv[5]); printf("lowchan = %d\n", lowchan);
	int highchan = atoi(argv[6]); printf("highchan = %d\n", highchan);
	int RFIF = atoi(argv[7]); printf("RFI-F method = %d\n", RFIF);
	int RFIT = atoi(argv[8]); printf("RFI-T method = %d\n", RFIT);
	float numSigmaF = atof(argv[9]); printf("numSigma-F = %f\n", numSigmaF);
	float numSigmaT = atof(argv[10]); printf("numSigma-T = %f\n", numSigmaT);
	int freqSmoothing = atoi(argv[11]); printf("freqSmoothing = %d\n", freqSmoothing);
	int uvDenoising = atoi(argv[12]); printf("uvDenoising = %d\n", uvDenoising);
	float uvDenoisingTau = atof(argv[13]); printf("uvDenoisingTau = %f\n", uvDenoisingTau);
	float uvDenoisingLambda = atof(argv[14]); printf("uvDenoisingLambda = %f\n", uvDenoisingLambda);
	float hidrogenfreq = atof(argv[15]); printf("hidrogenfreq = %f\n", hidrogenfreq);
	float hidrogenband = atof(argv[16]); printf("hidrogenband = %f\n", hidrogenband);
	int calskyfiles = atoi(argv[17]); printf("calskyfiles = %d\n", calskyfiles);
	int annfiles = atoi(argv[18]); printf("annfiles = %d\n", annfiles);
	int fit_smooth = atoi(argv[19]); printf("fit/smooth selection= %d\n", fit_smooth);
	int twindow = atoi(argv[20]); printf("smoothing window in time= %d\n", twindow);
        int fwindow = atoi(argv[21]); printf("smoothing window in frequency or max iterations for the rolling smooth function= %d\n", fwindow);

	// The imaging window. This helps reject highcal datapoints which mes up fitting
	float RAmin = atof(argv[22]); printf("RA min = %f\n",RAmin);
	float RAmax = atof(argv[23]); printf("RA max = %f\n",RAmax);
	float DECmin = atof(argv[24]); printf("DEC min = %f\n",DECmin);
	float DECmax = atof(argv[25]); printf("DEC max = %f\n",DECmax);

	char *day = argv[26];

	FILE * BadChannelsFile; BadChannelsFile = fopen("BadChannels.list", "r");

	int *badchannels = (int*)calloc(MAX_CHANNELS, sizeof(int));
	
	if(BadChannelsFile != NULL)
		{
		int chan1, chan2;
		do 
			{
			fscanf(BadChannelsFile, "%d %d\n", &chan1, &chan2);
			printf("Excluded: chan1 = %d chan2 = %d\n", chan1, chan2);
			for(i=chan1; i<=chan2; i++)
				{
				badchannels[i] = 1;
				}
			}
		while(!feof(BadChannelsFile));	
		fclose(BadChannelsFile);
		}


	// read and process bad data file
	FILE *baddatafile = fopen("baddata","rt");
	char header[80+1];
	float lowRA = 0.0, highRA = 0.0;
	int badlowchan = 0, badhighchan = 0, bad = 0;
	char badmjd[6];
	char badbeam[8]={0};
		
	// Handle the Days.list file
	FILE * dayfile; dayfile = fopen("Days.list", "r");
	int numdirs;
	char ** datedirs; 
	if(dayfile != NULL)
		{
		printf("Found file: \"Days.list\"\n");
		numdirs = jsd_line_count(dayfile);	
		datedirs = (char **) malloc(sizeof(char*) * numdirs);
		char tempstr[15];
		for(i = 0; i < numdirs; i++)
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

	// if we give a day via command line, just run that one day
	if( day != NULL ) {
		strncpy(datedirs[0], day, sizeof(day));
		numdirs = 1;
	}

	printf("Found %i data directories in %s\n", numdirs, datadirname);
	printf("The following days will be processed:\n");
	for(i=0; i < numdirs; i++) printf("%s,", datedirs[i]); printf("\n");

	char *subdir;
	int beamcounter; 
	char subdirname[5+1]; 
	mode_t mode = S_IRWXU|S_IRWXG|S_IRWXO;
	for (j = 0; j < numdirs; j++)
		{
		const char *datedir = datedirs[j];
		mkdir(datedir, mode);
		chdir(datedir);
		printf("***************************************************\n"); 
		printf("Now processing day: %s...\n", datedirs[j]); 
		printf("***************************************************\n");
		
		for(beamcounter = 0; beamcounter < 7; beamcounter++) 
			{
			if(beam[beamcounter]=='1')
			{

	                        // Read the Tcal correction file for the right beam
//        	                float Tcalx[MAX_CHANNELS];
//                	        float Tcaly[MAX_CHANNELS];
//	                        char tcalfile[32];
//        	                sprintf(tcalfile,"../Tcal%d.dat",beamcounter);
//                	        read_tcal(tcalfile, Tcalx, Tcaly);

				sprintf(subdirname,"beam%d", beamcounter);
				subdir = subdirname;
				mkdir(subdir, mode);
				chdir(subdir);
				bad = 0;

				if(baddatafile != NULL)
				{

				rewind( baddatafile );
				//read out the #header
				fgets(header,80,baddatafile);
				while(!feof(baddatafile))
					{
					// read out a line of bad data file

					//throws a warning is this correct !!
					fscanf(baddatafile,"%s %s", &badmjd, badbeam);
					//printf("bad beam pattern: %s %d %d %f %f\n", badbeam, badlowchan, badhighchan, lowRA, highRA);
					// check to see if mjd beam and channel match the current ones, if so stop reading further
					// we now have the RA range for the bad data for this day

					//printf( "\n ** %s %s", badmjd, datedir);
					//printf( "beam compare = %s", badbeam );

					if( (atoi(badmjd) == atoi(datedir)) && (badbeam[beamcounter] == '1') )
						{
						bad = 1;
						break;
						//for implementation read code at line no 361
						}
					}

				}

				if ( bad ) {
					printf( "\nSkipping");
					//printf("Skipping processing of day %s beam %s due to entry in bad data file", subdir, beamcounter);
				}
				else {
					printf("\nProcessing: %s\n",subdirname);
					process_dataset(field, datadirname, datedir, subdir, beamcounter, band, lowchan, highchan, \
									RFIF, RFIT, numSigmaF, numSigmaT, freqSmoothing, \
									uvDenoising, uvDenoisingTau, uvDenoisingLambda, \
									hidrogenfreq, hidrogenband, calskyfiles, annfiles, fit_smooth, \
									twindow, fwindow, badchannels, RAmin,RAmax,DECmin,DECmax);
				}

				printf("---------------------------------------------------\n");
				chdir("..");


				}
			}
		chdir("..");
		}
	if(baddatafile != NULL) fclose(baddatafile);
	for (i=0; i<numdirs; i++) free(datedirs[i]);
	free(datedirs); //free(subdir);
	free(badchannels);
	printf("Total computation time = %f\n", ((double)clock() - time0)/CLOCKS_PER_SEC);

	gettimeofday(&time, NULL);
	double sec1 = (double)time.tv_sec + (double)time.tv_usec * .000001;
	printf("Total wall time = %f\n", sec1-sec0);
	return EXIT_SUCCESS;
}



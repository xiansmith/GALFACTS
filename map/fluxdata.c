#include "map.h"
#include "fluxdata.h"
#include "jsd/jsd_futil.h"
#include <string.h>
#include <values.h>
#include <math.h>
#include <fcntl.h>

extern int multibeam;
/**********************************************************************/
int fluxrecord_read_binary(FluxRecord * pRec, FILE * file)
{
int k = 0;

k += fread(&pRec->RA, sizeof(float), 1, file);
k += fread(&pRec->DEC, sizeof(float), 1, file);
k += fread(&pRec->AST, sizeof(float), 1, file);
k += fread(&pRec->stokes.I, sizeof(float), 1, file);
k += fread(&pRec->stokes.Q, sizeof(float), 1, file);
k += fread(&pRec->stokes.U, sizeof(float), 1, file);
k += fread(&pRec->stokes.V, sizeof(float), 1, file);
//printf("RA %f\n", pRec->RA);

return k;
}
/**********************************************************************/

/**********************************************************************/
FluxWappData * fluxwappdata_alloc(const char *wapp, char **days, int numDays)
{
int i,j;
FluxWappData * wappdata;

wappdata = (FluxWappData*) malloc(sizeof(FluxWappData));
strncpy(wappdata->wapp, wapp, WAPP_LEN);
wappdata->numDays = numDays;
wappdata->daydata = (FluxDayData*) malloc(sizeof(FluxDayData) * numDays);
wappdata->scanDayData = (ScanDayData*) malloc(sizeof(ScanDayData) * numDays);
for(i=0; i<numDays; i++) 
	{
	if(!strcmp(wapp,"multibeam"))
		{
		j = i/7;
	    strncpy(wappdata->daydata[i].mjd, days[j], MJD_LEN);
		}
	else strncpy(wappdata->daydata[i].mjd, days[i], MJD_LEN);
    wappdata->daydata[i].numRecords = 0;
    wappdata->daydata[i].records = NULL;
    wappdata->scanDayData[i].numScans = 0;
    wappdata->scanDayData[i].scans = NULL;
    }
return wappdata;

}
/**********************************************************************/
void fluxwappdata_free(FluxWappData * wappdata)
{
int i;
if(wappdata == NULL) return;
if(wappdata->daydata != NULL) 
	{
    for(i=0; i<wappdata->numDays; i++) 
		{
        if(wappdata->daydata[i].records != NULL) 
			{
            free(wappdata->daydata[i].records);
            wappdata->daydata[i].records = NULL;
            }
        }
	free(wappdata->daydata);
	wappdata->daydata = NULL;
    }
free(wappdata);
}
/**********************************************************************/
int fluxdaydata_read_binary(const char *field, FluxDayData *daydata, FILE *infile, int beam, int chan, int day)
{
int numRecords;
char header[80+1];
float RAmax = FLT_MIN;
float RAmin = FLT_MAX;

fread(&numRecords, sizeof(int), 1, infile); 
if(daydata->records != NULL) free(daydata->records);
daydata->records = (FluxRecord*) malloc(numRecords * sizeof(FluxRecord));
if(daydata->records == NULL)
	{
	printf("ERROR: malloc failed !\n");
	exit(0);
	}
static float RAp[200000],DECp[200000],ASTp[200000];


// read and process bad data file

FILE *baddatafile = fopen("Baddata.list","rt");
if(baddatafile != NULL)
printf(":: Opened with fd %d\n",*baddatafile);
float lowRA = 0.0, highRA = 0.0;

int badmjd = 0, badlowchan = 0, badhighchan = 0, bad = 0;
char badbeam[8]={0};

if(baddatafile != NULL)
	{
	//read out the #header
	//printf(":: BAD file found.\n");
	fgets(header,80,baddatafile);
	while(!feof(baddatafile))
		{
		// read out a line of bad data file
		fscanf(baddatafile,"%d %s %d %d %f %f", &badmjd, badbeam, &badlowchan, &badhighchan, &lowRA, &highRA);
		//printf("bad beam pattern: %s %d %d %f %f\n", badbeam, badlowchan, badhighchan, lowRA, highRA);
		// check to see if mjd beam and channel match the current ones, if so stop reading further
		// we now have the RA range for the bad data for this day
		if((badmjd == day) && (badbeam[beam] == '1') && (chan >= badlowchan) && (chan <= badhighchan))
			{
			//fclose(baddatafile);
			bad = 1;
			break;
			//for implementation read code at line no 361
			}
		}
	}
else
{
//	printf(":: BAD file not found.\n");
	bad = 0;
}


////////////

int k = 0;
int flag = 0;
int m = 0;
while(!feof(infile) && k<numRecords)
    {
	FluxRecord *pRec = &daydata->records[k];
  //  int num = fluxrecord_read_binary(pRec, infile);
    int num = fluxrecord_read_binary(pRec, infile);
    if(num==0) continue;
	//fix for certain missing data in beam 0 (a real pain!)

// Just for N1
// ----------------------
if(field[0] == 'N' && field[1] == '1')
	{
	if(ASTp[m] != daydata->records[k].AST && beam != 0) continue;
	//pointing fix once and for all dec > 18
	if(beam == 0)
		{
		RAp[m] = pRec->RA;
		DECp[m] = pRec->DEC;
		ASTp[m] = pRec->AST;
		m++;
		if(m > 200000) printf("Pointing fix array overrun !\n");
		}
	}
// ---------------------- N1 ends

	if(num == 7)
		{
        float RA = daydata->records[k].RA;
        float DECR = daydata->records[k].DEC*M_PI/180.0;
        if(RA > RAmax) RAmax = RA;
        if (RA < RAmin) RAmin = RA;

// Just for N1
// ----------------------
	if(field[0] == 'N' && field[1] == '1')
		{
		switch(beam)
			{
		    case 0: break;
		    case 1:
                    daydata->records[k].RA = RAp[m] + 2.7417/(60*cos(DECR));
					m++;
					break;
		    case 2:
					daydata->records[k].RA =RAp[m]+ 5.4833/(60*cos(DECR));
					m++;
					break;
		    case 3:
                    daydata->records[k].RA =RAp[m]+ 2.7417/(60*cos(DECR));
					m++;
					break;
		    case 4:
                    daydata->records[k].RA =RAp[m] -2.7417/(60*cos(DECR));
					m++;
					break;
		    case 5:
                    daydata->records[k].RA =RAp[m] -5.4833/(60*cos(DECR));
					m++;
					break;
		    case 6:
                    daydata->records[k].RA =RAp[m] -2.7417/(60*cos(DECR));
					m++;
					break;
		    default:
					printf("WARN: invalid beam (%i) specified.\n", beam);
			}
	}
	// ---------------------- N1 ends
	if(bad && daydata->records[k].RA > highRA )
	{
		if( baddatafile == NULL ) break;

		if(!feof(baddatafile))
		{
			// read out a line of bad data file
			fscanf(baddatafile,"%d %s %d %d %f %f", &badmjd, badbeam, &badlowchan, &badhighchan, &lowRA, &highRA);
			// check to see if mjd beam and channel match the current ones, if so stop reading further
			// we now have the RA range for the bad data for this day
			if((badmjd == day) && (badbeam[beam] == '1') && (chan >= badlowchan) && (chan <= badhighchan) )
			{
				bad = 1;
			} else {
				bad = 0;
				//printf("Closing file\n");
				fclose( baddatafile );
			}

		} else
		{
			bad = 0;
			//printf("Closing file\n");
			fclose (baddatafile );
		}
	}

        if(bad && daydata->records[k].RA>=lowRA && daydata->records[k].RA<=highRA)
            {
            //if( daydata->records[k].RA > 117 ) printf("bad data at RA %f\n", daydata->records[k].RA );
//		printf(":BAD: %d %s %d %f\n", badmjd, badbeam, chan, daydata->records[k].RA);
            daydata->records[k].stokes.I = NAN;
            daydata->records[k].stokes.Q = NAN;
            daydata->records[k].stokes.U = NAN;
            daydata->records[k].stokes.V = NAN;
            }

		k++;
		}
	else if(num <= 0) break;
	else
		{
		printf("ERROR: flux file record only read %i fields\n", num);
		exit(1);
		}
    }
    daydata->numRecords = k;
    daydata->RAmin = RAmin;
    daydata->RAmax = RAmax;

	// check if not closed to avoid leaking 
/*	if ( fcntl(fileno(baddatafile), F_GETFD) != -1 )  {
		fclose(baddatafile);
		printf(":BAD: Closed file\n");
	}
*/
    return k;
}

/**********************************************************************/
int fluxdaydata_read_binary_single_file(const char *field, FluxDayData *daydata, FILE *infile, FILE *configfile, int beam, int chan, int day)
{
	int numRecords;
	char header[80+1];
	float RAmax = FLT_MIN;
	float RAmin = FLT_MAX;
	int k = 0;
	int channelcount;
	int configlines;
	int *records;
	int numRead;
	int cfgChan, cfgStart, cfgNumRecords;
	char filename[64];

	// for single file processing, we need the config file or fail and exit
	// even missing channels will have NANs to read
	if( configfile != NULL ) {
		do
		{
			numRead = fscanf(configfile, "%d %d %d\n", &cfgChan, &cfgStart, &cfgNumRecords);

			if( numRead == 3 ) {
				if( cfgChan == chan ) {
					break;
					//found the channel we are looking for
				}
			} else {
				printf("Error in the fluxtime config file. Read %d records instead of 3. Exiting.\n", numRead);
				exit(1);
			}
		}
		while(!feof(configfile) && (cfgChan < chan ));

		if( cfgChan != chan ) {
			printf("Could not find the channel data in the fluxdata file. Can not recover\n");
			printf("Looking for chan = %d, cfgChan %d\n", chan, cfgChan);
			exit(1);
		}

	} else {
		// average image is old method, so don't exit
		if( chan != 0 ) {
			printf("No fluxtime config found. Can not recover, exiting\n");
			printf("couldn't find config file for %d %d %d \n", day, beam, chan);
			exit(1);
		}
		//exit(1);
	}

	//printf( "Get data from single binary file. Fetching channel %d at position %d for %d numRecords\n", cfgChan, cfgStart, cfgNumRecords);

	// we could assume numRecords never changes, but just in case
	numRecords = cfgNumRecords;
	//printf("numRecords after is %d\n", numRecords );

	if(daydata->records != NULL) free(daydata->records);
	daydata->records = (FluxRecord*) malloc(numRecords * sizeof(FluxRecord));
	if(daydata->records == NULL)
		{
		printf("ERROR: malloc failed !\n");
		exit(0);
	}
	static float RAp[200000],DECp[200000],ASTp[200000];


// read and process bad data file

FILE *baddatafile = fopen("Baddata.list","rt");
if(baddatafile != NULL)
printf(":Day %d: Opened with fd %d\n",day,*baddatafile);

float lowRA = 0.0, highRA = 0.0;

int badmjd = 0, badlowchan = 0, badhighchan = 0, bad = 0;
char badbeam[8]={0};

if(baddatafile != NULL)
	{
	//printf(":: BAD file found.\n");
	//read out the #header
	fgets(header,80,baddatafile);
	while(!feof(baddatafile))
		{
		// read out a line of bad data file
		fscanf(baddatafile,"%d %s %d %d %f %f", &badmjd, badbeam, &badlowchan, &badhighchan, &lowRA, &highRA);
		//printf("bad beam pattern: %s %d %d %f %f\n", badbeam, badlowchan, badhighchan, lowRA, highRA);
		// check to see if mjd beam and channel match the current ones, if so stop reading further
		// we now have the RA range for the bad data for this day
		if((badmjd == day) && (badbeam[beam] == '1') && (chan >= badlowchan) && (chan <= badhighchan))
			{
			//fclose(baddatafile);
			bad = 1;
			break;
			//for implementation read code at line no 361
			}
		}
	}
else
{
	//printf(":Day %d: BAD file not found.\n",day);
	bad = 0;
}

////////////

		// seek to right record, without overflowing the long
		long recordSize = sizeof(float) * 7;
		int ret = fseek(infile, recordSize * cfgStart, SEEK_SET);

		if ( ret != 0 ) {
			printf("ERROR: Error seeking in binary file with error %d\n", ret );
		}


k = 0;
int flag = 0;
int m = 0;
while( !feof(infile) &&  k<numRecords)
    {
	FluxRecord *pRec = &daydata->records[k];

	// read the floats from the file after seeking to channel
	int num = fluxrecord_read_binary(pRec, infile);

	if(num==0) continue;
	//fix for certain missing data in beam 0 (a real pain!)

// Just for N1
// ----------------------
if(field[0] == 'N' && field[1] == '1')
	{
	if(ASTp[m] != daydata->records[k].AST && beam != 0) continue;
	//pointing fix once and for all dec > 18
	if(beam == 0)
		{
		RAp[m] = pRec->RA;
		DECp[m] = pRec->DEC;
		ASTp[m] = pRec->AST;
		m++;
		if(m > 200000) printf("Pointing fix array overrun !\n");
		}
	}
// ---------------------- N1 ends
 if(field[0] == 'N' && field[1] == '4') {
                                RAmax = 390.0;
                        }

                        if(field[0] == 'S' && field[1] == '4') {
                                RAmax = 415.0;
                        }

	
	if(num == 7)
		{
        float RA = daydata->records[k].RA;
        float DECR = daydata->records[k].DEC*M_PI/180.0;
        if(RA > RAmax) RAmax = RA;
        if (RA < RAmin) RAmin = RA;

// Just for N1
// ----------------------
	if(field[0] == 'N' && field[1] == '1')
		{
		switch(beam)
			{
		    case 0: break;
		    case 1:
                    daydata->records[k].RA = RAp[m] + 2.7417/(60*cos(DECR));
					m++;
					break;
		    case 2:
					daydata->records[k].RA =RAp[m]+ 5.4833/(60*cos(DECR));
					m++;
					break;
		    case 3:
                    daydata->records[k].RA =RAp[m]+ 2.7417/(60*cos(DECR));
					m++;
					break;
		    case 4:
                    daydata->records[k].RA =RAp[m] -2.7417/(60*cos(DECR));
					m++;
					break;
		    case 5:
                    daydata->records[k].RA =RAp[m] -5.4833/(60*cos(DECR));
					m++;
					break;
		    case 6:
                    daydata->records[k].RA =RAp[m] -2.7417/(60*cos(DECR));
					m++;
					break;
		    default:
					printf("WARN: invalid beam (%i) specified.\n", beam);
			}
	}
// ---------------------- N1 ends

	if(bad && daydata->records[k].RA > highRA )
	{
		if(!feof(baddatafile))
		{
			// read out a line of bad data file
			fscanf(baddatafile,"%d %s %d %d %f %f", &badmjd, badbeam, &badlowchan, &badhighchan, &lowRA, &highRA);
			//printf("AFTER bad beam pattern: %s %d %d %f %f\n", badbeam, badlowchan, badhighchan, lowRA, highRA);
			// check to see if mjd beam and channel match the current ones, if so stop reading further
			// we now have the RA range for the bad data for this day
			if((badmjd == day) && (badbeam[beam] == '1') && (chan >= badlowchan) && (chan <= badhighchan) )
			{
				bad = 1;
			} else {
				bad = 0;
			    printf("Closing file\n");
				fclose( baddatafile );
			}

		} else
		{
			bad = 0;
			    printf("Closing file\n");
			fclose (baddatafile );
		}
	}

        if(bad && daydata->records[k].RA>=lowRA && daydata->records[k].RA<=highRA)
            {
		//printf(":BAD: %d %s %d %f\n", badmjd, badbeam, chan, daydata->records[k].RA);
		//exit(1);
            //if( daydata->records[k].RA > 117 ) printf("bad data at RA %f\n", daydata->records[k].RA );
            daydata->records[k].stokes.I = NAN;
            daydata->records[k].stokes.Q = NAN;
            daydata->records[k].stokes.U = NAN;
            daydata->records[k].stokes.V = NAN;
            }
	
	
		 if(daydata->records[k].RA < RAmax)
			k++;
		}
	else if(num <= 0) break;
	else
		{
		printf("ERROR: flux file record only read %i fields\n", num);
		exit(1);
		}
    }
	daydata->numRecords = k;
    daydata->RAmin = RAmin;
    daydata->RAmax = RAmax;
  
	// check if not closed to avoid leaking 
/*	if ( fcntl(fileno(baddatafile), F_GETFD) != -1 )  {
		fclose(baddatafile);
		//printf(":BAD: Closed file\n");
	}
*/	
	return k;
}
/**********************************************************************/

/**********************************************************************/
void fluxwappdata_readchan_binary(const char *field, int band, FluxWappData * wappdata, int chan, int id, int avg, float decmin, float decmax)
{
int  m,j,n, i, numread, flag = 0, navg;
FILE *infile, *configfile;
char beamno[6];
float weight = 0;

float interp[10] = {1.23376620, 1.18181813, 1.12987018, 1.07792211, 1.02597404, 0.97402596, 0.92207789, 0.87012988, 0.81818181, 0.76623374};

if(avg == 0) navg = 1; else navg = avg;

for(m=0; m<wappdata->numDays; m++)
	{
	printf("Reading Day %d\n",m);
	fflush(stdout);
	FluxDayData * daydata = &wappdata->daydata[m];	
	FluxDayData * tempdata_avg = NULL;
	if(tempdata_avg == NULL) 
		{
		tempdata_avg = malloc(sizeof(FluxDayData));
		tempdata_avg->records = NULL;
		}
	flag = 0;
	for(n=0; n<navg; n++)
		{
		FluxDayData * tempdata = NULL;
		if(tempdata == NULL) 
			{
			tempdata = malloc(sizeof(FluxDayData));
			tempdata->records = NULL;
			}
		char filename[64+1];
		char configfilename[64+1];
		if(!strcmp(wappdata->wapp,"multibeam"))
		{
			j = m%7;
			sprintf(beamno,"beam%d",j);
			if(id == CLEAN)
			{
				sprintf(filename, "%s/%s/balance.dat", daydata->mjd, beamno);
				sprintf(configfilename, "%s/%s/balance.dat_cfg", daydata->mjd, beamno);

			}
			if(id == BASKETWEAVE)
			{
				sprintf(filename, "%s/%s/fluxtime.dat", daydata->mjd, beamno);
				sprintf(configfilename, "%s/%s/fluxtime.dat_cfg", daydata->mjd, beamno);
			}
		}
		else
		{
			j = atoi(&wappdata->wapp[4]);
			sprintf(beamno,"beam%d",j);
			if(id == CLEAN)
			{
				sprintf(filename, "%s/%s/balance.dat", daydata->mjd, wappdata->wapp);
				sprintf(configfilename, "%s/%s/balance.dat_cfg", daydata->mjd, wappdata->wapp);

			}
			if(id == BASKETWEAVE)
			{
				sprintf(filename, "%s/%s/fluxtime.dat", daydata->mjd, wappdata->wapp);
				sprintf(configfilename, "%s/%s/fluxtime.dat_cfg", daydata->mjd, wappdata->wapp);
			}
		}
		infile = fopen(filename, "rb");
		configfile = fopen(configfilename, "r");
		if( configfile == NULL ) printf("ERROR: configfile is NULL\n");
		if(infile != NULL) 
			{
			// we are looking for average image, read average.dat instead
			if( chan == 0 && avg == 0 ) {
				fclose( infile );
				sprintf(filename, "%s/%s/average.dat", daydata->mjd, beamno);
				FILE * avgfile = fopen(filename, "rb");

				if( avgfile != NULL )
				{
					numread = fluxdaydata_read_binary(field, tempdata, infile, j, 0, atoi(daydata->mjd));
				}
				else
				{
					printf("Error. Looking for %s for average image, failed. Exitting.\n",filename);
					exit(1);

				}
			}
			else
			{
				//numread = fluxdaydata_read_binary(field, tempdata, infile, j, chan+n, atoi(daydata->mjd));
				numread = fluxdaydata_read_binary_single_file(field, tempdata, infile, configfile, j, chan+n, atoi(daydata->mjd));
			}
			printf("Opened file %s\n",filename);
			printf("read %d records.\n",tempdata->numRecords );
			fclose(infile);
			fclose(configfile);
			if(n==0)
				{			
				if(tempdata_avg->records != NULL) free(tempdata_avg->records);	
				tempdata_avg->records = (FluxRecord*) malloc(tempdata->numRecords * sizeof(FluxRecord));
				tempdata_avg->numRecords = tempdata->numRecords;
				tempdata_avg->RAmin = tempdata->RAmin;
				tempdata_avg->RAmax = tempdata->RAmax;	
				for(i=0; i<numread; i++)
					{
					tempdata_avg->records[i].RA = tempdata->records[i].RA;
					tempdata_avg->records[i].DEC = tempdata->records[i].DEC;
					tempdata_avg->records[i].AST = tempdata->records[i].AST;
					tempdata_avg->records[i].stokes.I = 0.0;
					tempdata_avg->records[i].stokes.Q = 0.0;
					tempdata_avg->records[i].stokes.U = 0.0;
					tempdata_avg->records[i].stokes.V = 0.0;
					tempdata_avg->records[i].count = 0;
					tempdata_avg->records[i].weight = 0;
					}
				}	
			for(i=0; i<numread; i++)
			{
				if(isfinite(tempdata->records[i].stokes.I))
				{
					if( band == 1 )
					{
						tempdata_avg->records[i].stokes.I += tempdata->records[i].stokes.I*interp[n];
						tempdata_avg->records[i].stokes.Q += tempdata->records[i].stokes.Q*interp[n];
						tempdata_avg->records[i].stokes.U += tempdata->records[i].stokes.U*interp[n];
						tempdata_avg->records[i].stokes.V += tempdata->records[i].stokes.V*interp[n];
						tempdata_avg->records[i].weight += interp[n];
						tempdata_avg->records[i].count++;
					}
					else
					{
						tempdata_avg->records[i].stokes.I += tempdata->records[i].stokes.I;
						tempdata_avg->records[i].stokes.Q += tempdata->records[i].stokes.Q;
						tempdata_avg->records[i].stokes.U += tempdata->records[i].stokes.U;
						tempdata_avg->records[i].stokes.V += tempdata->records[i].stokes.V;
						tempdata_avg->records[i].count++;
					}
				}
			}
			free(tempdata->records);
			}
		else
			{
			flag++;
			printf("Error could not open file %s\n",filename);
			}
		free(tempdata);
		}
	if(flag < navg)
		{
		for(i=0; i<numread; i++)
			{
			if(tempdata_avg->records[i].count > 0)
				{
					if( band == 1 )
					{
						tempdata_avg->records[i].stokes.I /= tempdata_avg->records[i].weight;
						tempdata_avg->records[i].stokes.Q /= tempdata_avg->records[i].weight;
						tempdata_avg->records[i].stokes.U /= tempdata_avg->records[i].weight;
						tempdata_avg->records[i].stokes.V /= tempdata_avg->records[i].weight;
					}
					if ( band == 0 )
					{
						tempdata_avg->records[i].stokes.I /= tempdata_avg->records[i].count;
						tempdata_avg->records[i].stokes.Q /= tempdata_avg->records[i].count;
						tempdata_avg->records[i].stokes.U /= tempdata_avg->records[i].count;
						tempdata_avg->records[i].stokes.V /= tempdata_avg->records[i].count;
					}
				}
				else
				{
					tempdata_avg->records[i].stokes.I = NAN;
					tempdata_avg->records[i].stokes.Q = NAN;
					tempdata_avg->records[i].stokes.U = NAN;
					tempdata_avg->records[i].stokes.V = NAN;
				}
			}
		daydata->records = (FluxRecord*) malloc(tempdata_avg->numRecords * sizeof(FluxRecord));
		daydata->RAmin = tempdata_avg->RAmin;
		daydata->RAmax = tempdata_avg->RAmax;
		int jj = 0; 
		for(i=0; i < numread; i++)
			{
			int moonflag = 0;
			if(isfinite(tempdata_avg->records[i].stokes.I))
				{
				daydata->records[jj].stokes.I = tempdata_avg->records[i].stokes.I;
				daydata->records[jj].stokes.Q = tempdata_avg->records[i].stokes.Q;
				daydata->records[jj].stokes.U = tempdata_avg->records[i].stokes.U;
				daydata->records[jj].stokes.V = tempdata_avg->records[i].stokes.V;
				daydata->records[jj].RA = tempdata_avg->records[i].RA;
				daydata->records[jj].DEC = tempdata_avg->records[i].DEC;
				daydata->records[jj].AST = tempdata_avg->records[i].AST;				

// Just for N1
// ----------------------
		if(field[0] == 'N' && field[1] == '1')
			{
				if(!strcmp(daydata->mjd,"54811") && (tempdata_avg->records[i].AST > 9505)  && (tempdata_avg->records[i].AST < 9580)) moonflag = 1;
				if(!strcmp(daydata->mjd,"54812") && (tempdata_avg->records[i].AST > 13475)  && (tempdata_avg->records[i].AST < 13555)) moonflag = 1;
				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 42.0 && tempdata_avg->records[i].RA > 40.5))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 45.5 && tempdata_avg->records[i].RA > 44.5))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 48.25 && tempdata_avg->records[i].RA > 47.0))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 77.0 && tempdata_avg->records[i].RA > 75.0))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 79.25 && tempdata_avg->records[i].RA > 78.0))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54820") && ((tempdata_avg->records[i].RA < 66.25 && tempdata_avg->records[i].RA > 65.75))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54820") && ((tempdata_avg->records[i].RA < 79.75 && tempdata_avg->records[i].RA > 78.75))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54820") && ((tempdata_avg->records[i].RA < 88.25 && tempdata_avg->records[i].RA > 87.0))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54813") && ((tempdata_avg->records[i].RA < 88.00 && tempdata_avg->records[i].RA > 87.5))) moonflag = 1;
				if(!strcmp(daydata->mjd,"54785") && ((tempdata_avg->records[i].RA < 79.75 && tempdata_avg->records[i].RA > 79.0))) moonflag = 1;
			}
// ---------------------- N1 ends

//				Remove alien spacecraft !!
				if(i>0)
					{
					if((tempdata_avg->records[i].DEC < decmax && tempdata_avg->records[i].DEC > decmin  && fabs(tempdata_avg->records[i].DEC - tempdata_avg->records[i-1].DEC) > 0.0003 && (tempdata_avg->records[i].RA - tempdata_avg->records[i-1].RA) > 0 && !moonflag)) 
					jj++;
					}
				if(i==0) jj++;				
				}
			}		
		daydata->numRecords = jj;
		free(tempdata_avg->records);
		}
	free(tempdata_avg);
	} 
}
/**********************************************************************/
//void fluxwappdata_readchan_binary1(const char *field, FluxWappData * wappdata, int chan, int id, int avg, float decmin, float decmax) //band 1
//{
//int  m,j,n, i, numread, flag = 0, navg;
//FILE *infile;
//char beamno[6];
//float weight = 0;
//
//
//float interp[10] = {1.23376620, 1.18181813, 1.12987018, 1.07792211, 1.02597404, 0.97402596, 0.92207789, 0.87012988, 0.81818181, 0.76623374};
//
//if(avg == 0) navg = 1; else navg = avg;
//
//for(m=0; m<wappdata->numDays; m++)
//	{
//	FluxDayData * daydata = &wappdata->daydata[m];
//	FluxDayData * tempdata_avg = NULL;
//	if(tempdata_avg == NULL)
//		{
//		tempdata_avg = malloc(sizeof(FluxDayData));
//		tempdata_avg->records = NULL;
//		}
//	flag = 0;
//	for(n=0; n<navg; n++)
//		{
//		FluxDayData * tempdata = NULL;
//		if(tempdata == NULL)
//			{
//			tempdata = malloc(sizeof(FluxDayData));
//			tempdata->records = NULL;
//			}
//		char filename[64+1];
//		if(!strcmp(wappdata->wapp,"multibeam"))
//			{
//			j = m%7;
//			sprintf(beamno,"beam%d",j);
//			if(id == CLEAN) sprintf(filename, "%s/%s/balanceB%04i.dat", daydata->mjd, beamno, chan+n);
//			if(id == BASKETWEAVE) sprintf(filename, "%s/%s/fluxtime%04i.dat", daydata->mjd, beamno, chan+n);
//			}
//		else
//			{
//			j = atoi(&wappdata->wapp[4]);
//			if(id == CLEAN) sprintf(filename, "%s/%s/balanceB%04i.dat", daydata->mjd, wappdata->wapp, chan+n);
//			if(id == BASKETWEAVE) sprintf(filename, "%s/%s/fluxtime%04i.dat", daydata->mjd, wappdata->wapp, chan+n);
//			}
//		infile = fopen(filename, "rb");
//		if(infile != NULL)
//			{
//			printf("Opened file %s\n",filename);
//			numread = fluxdaydata_read_binary(field, tempdata, infile, j, chan+n, atoi(daydata->mjd));
//			printf("read %d records.\n",tempdata->numRecords );
//			fclose(infile);
//			if(n==0)
//				{
//				if(tempdata_avg->records != NULL) free(tempdata_avg->records);
//				tempdata_avg->records = (FluxRecord*) malloc(tempdata->numRecords * sizeof(FluxRecord));
//				tempdata_avg->numRecords = tempdata->numRecords;
//				tempdata_avg->RAmin = tempdata->RAmin;
//				tempdata_avg->RAmax = tempdata->RAmax;
//				for(i=0; i<numread; i++)
//					{
//					tempdata_avg->records[i].RA = tempdata->records[i].RA;
//					tempdata_avg->records[i].DEC = tempdata->records[i].DEC;
//					tempdata_avg->records[i].AST = tempdata->records[i].AST;
//					tempdata_avg->records[i].stokes.I = 0.0;
//					tempdata_avg->records[i].stokes.Q = 0.0;
//					tempdata_avg->records[i].stokes.U = 0.0;
//					tempdata_avg->records[i].stokes.V = 0.0;
//					tempdata_avg->records[i].count = 0;
//					tempdata_avg->records[i].weight = 0;
//					}
//				}
//			for(i=0; i<numread; i++)
//				{
//				if(isfinite(tempdata->records[i].stokes.I))
//					{
//					tempdata_avg->records[i].stokes.I += tempdata->records[i].stokes.I*interp[n];
//					tempdata_avg->records[i].stokes.Q += tempdata->records[i].stokes.Q*interp[n];
//					tempdata_avg->records[i].stokes.U += tempdata->records[i].stokes.U*interp[n];
//					tempdata_avg->records[i].stokes.V += tempdata->records[i].stokes.V*interp[n];
//					tempdata_avg->records[i].weight += interp[n];
//					tempdata_avg->records[i].count++;
//					}
//				}
//			free(tempdata->records);
//			}
//		else
//			{
//			flag++;
//			printf("Error could not open file %s\n",filename);
//			}
//		free(tempdata);
//		}
//	if(flag < navg)
//		{
//		for(i=0; i<numread; i++)
//			{
//			if(tempdata_avg->records[i].count > 0)
//				{
//				tempdata_avg->records[i].stokes.I /= tempdata_avg->records[i].weight;
//				tempdata_avg->records[i].stokes.Q /= tempdata_avg->records[i].weight;
//				tempdata_avg->records[i].stokes.U /= tempdata_avg->records[i].weight;
//				tempdata_avg->records[i].stokes.V /= tempdata_avg->records[i].weight;
//				}
//			else
//				{
//				tempdata_avg->records[i].stokes.I = NAN;
//				tempdata_avg->records[i].stokes.Q = NAN;
//				tempdata_avg->records[i].stokes.U = NAN;
//				tempdata_avg->records[i].stokes.V = NAN;
//				}
//			}
//		daydata->records = (FluxRecord*) malloc(tempdata_avg->numRecords * sizeof(FluxRecord));
//		daydata->RAmin = tempdata_avg->RAmin;
//		daydata->RAmax = tempdata_avg->RAmax;
//		int jj = 0;
//		for(i=0; i < numread; i++)
//			{
//			int moonflag = 0;
//			if(isfinite(tempdata_avg->records[i].stokes.I))
//				{
//				daydata->records[jj].stokes.I = tempdata_avg->records[i].stokes.I;
//				daydata->records[jj].stokes.Q = tempdata_avg->records[i].stokes.Q;
//				daydata->records[jj].stokes.U = tempdata_avg->records[i].stokes.U;
//				daydata->records[jj].stokes.V = tempdata_avg->records[i].stokes.V;
//				daydata->records[jj].RA = tempdata_avg->records[i].RA;
//				daydata->records[jj].DEC = tempdata_avg->records[i].DEC;
//				daydata->records[jj].AST = tempdata_avg->records[i].AST;
//
//// Just for N1
//// ----------------------
//			if(field[0] == 'N' && field[1] == '1')
//				{
//				if(!strcmp(daydata->mjd,"54811") && (tempdata_avg->records[i].AST > 9505)  && (tempdata_avg->records[i].AST < 9580)) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54812") && (tempdata_avg->records[i].AST > 13475)  && (tempdata_avg->records[i].AST < 13555)) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 42.0 && tempdata_avg->records[i].RA > 40.5))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 45.5 && tempdata_avg->records[i].RA > 44.5))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 48.25 && tempdata_avg->records[i].RA > 47.0))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 77.0 && tempdata_avg->records[i].RA > 75.0))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54793") && ((tempdata_avg->records[i].RA < 79.25 && tempdata_avg->records[i].RA > 78.0))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54820") && ((tempdata_avg->records[i].RA < 66.25 && tempdata_avg->records[i].RA > 65.75))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54820") && ((tempdata_avg->records[i].RA < 79.75 && tempdata_avg->records[i].RA > 78.75))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54820") && ((tempdata_avg->records[i].RA < 88.25 && tempdata_avg->records[i].RA > 87.0))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54813") && ((tempdata_avg->records[i].RA < 88.00 && tempdata_avg->records[i].RA > 87.5))) moonflag = 1;
//				if(!strcmp(daydata->mjd,"54785") && ((tempdata_avg->records[i].RA < 79.75 && tempdata_avg->records[i].RA > 79.0))) moonflag = 1;
//				}
//// ---------------------- N1 ends
//
////				Remove alien spacecraft !!
//				if(i>0)
//					{
//					if((tempdata_avg->records[i].DEC < decmax && tempdata_avg->records[i].DEC > decmin  && fabs(tempdata_avg->records[i].DEC - tempdata_avg->records[i-1].DEC) > 0.0003 && (tempdata_avg->records[i].RA - tempdata_avg->records[i-1].RA) > 0 && !moonflag))
//					jj++;
//					}
//				if(i==0) jj++;
//				}
//			}
//		daydata->numRecords = jj;
//		free(tempdata_avg->records);
//		}
//	free(tempdata_avg);
//	}
//}
/**********************************************************************/
int fluxwappdata_writechan(FluxWappData * wappdata, int chan)
{
    int m;
    int count;

    count = 0;
    for(m=0; m<wappdata->numDays; m++)
		{
        int k;
        FILE *file;
        int numRecords;
        char filename[64+1];
		char tempstring[6];
        FluxDayData * daydata = &wappdata->daydata[m];
		//SSG
		if(multibeam)
			{
			sprintf(tempstring,"beam%d",m%7);
	        sprintf(filename, "%s/%s/balance%04i.dat", daydata->mjd, tempstring, chan);
			}
		else
	    sprintf(filename, "%s/%s/balance%04i.dat", daydata->mjd, wappdata->wapp, chan);
		//SSG
		file = fopen(filename, "w");
		if(file == NULL) 
			{
			printf("ERROR: can't open output file %s\n", filename);
			continue;
			}
		else //SSG
		fprintf(file, "#RA DEC AST I Q U V\n"); //SSG
        numRecords = daydata->numRecords;
        for (k=0; k<numRecords; k++) 
			{
            fluxrecord_write(&daydata->records[k], file);
			}
        fclose(file);
        count++;
		}
return count;
}
/**********************************************************************/
int fluxwappdata_writechan_binary(FluxWappData * wappdata, int chan)
{
    int m;
    int count;

    count = 0;
    for(m=0; m<wappdata->numDays; m++)
		{
        int k;
        FILE *file;
        int numRecords;
        char filename[64+1];
		char tempstring[6];
        FluxDayData * daydata = &wappdata->daydata[m];
		//SSG
		if(multibeam)
			{
			sprintf(tempstring,"beam%d",m%7);
	        sprintf(filename, "%s/%s/balance%04i.dat", daydata->mjd, tempstring, chan);
			}
		else
	    sprintf(filename, "%s/%s/balance%04i.dat", daydata->mjd, wappdata->wapp, chan);
		//SSG
		file = fopen(filename, "wb");
		if(file == NULL) 
			{
			printf("ERROR: can't open output file %s\n", filename);
			continue;
			}
		else //SSG
		//fprintf(file, "#RA DEC AST I Q U V\n"); //SSG
        numRecords = daydata->numRecords;
        for (k=0; k<numRecords; k++) 
			{
			fluxrecord_write_binary(&daydata->records[k], file);
			}
        fclose(file);
        count++;
		}
return count;
}
/**********************************************************************/
int fluxwappdata_writechan_binary_single(FluxWappData * wappdata, int chan)
{
    int m;
    int count;

    count = 0;
    for(m=0; m<wappdata->numDays; m++)
		{
        int k;
        FILE *file;
        FILE *configfile;
        int numRecords;
        char filename[64+1];
		char tempstring[6];
        char configfilename[64+1];
		FluxDayData * daydata = &wappdata->daydata[m];

        //SSG
		if(multibeam)
		{
			sprintf(tempstring,"beam%d",m%7);
	        sprintf(filename, "%s/%s/balance.dat", daydata->mjd, tempstring);
	        sprintf(configfilename, "%s/%s/balance.dat_cfg", daydata->mjd, tempstring);
		}
		else
		{
			sprintf(filename, "%s/%s/balance.dat", daydata->mjd, wappdata->wapp);
			sprintf(configfilename, "%s/%s/balance.dat_cfg", daydata->mjd, tempstring);
		}

		configfile = fopen(configfilename, "w");
		if(configfile == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
		} else {
			// write config file entry
			fprintf(configfile, "%d\n", numRecords);
		}
		fclose( configfile );

		//SSG
		file = fopen(filename, "wb");
		if(file == NULL)
			{
			printf("ERROR: can't open output file %s\n", filename);
			continue;
			}
		else //SSG
		//fprintf(file, "#RA DEC AST I Q U V\n"); //SSG

		numRecords = daydata->numRecords;
        for (k=0; k<numRecords; k++)
			{
			fluxrecord_write_binary(&daydata->records[k], file);
			}
        fclose(file);
        count++;
		}
return count;
}
/**********************************************************************/
int fluxrecord_write(FluxRecord * pRec, FILE * file)
{
return fprintf(file,"%2.8f %2.8f %8.2f %4.6f %4.6f %4.6f %4.6f\n",
            pRec->RA, pRec->DEC, pRec->AST,
            pRec->stokes.I, pRec->stokes.Q, pRec->stokes.U, pRec->stokes.V);
}
/**********************************************************************/
int fluxrecord_write_binary(FluxRecord * pRec, FILE * file)
{
int s = 0;
s += fwrite(&pRec->RA, sizeof(float), 1, file);
s += fwrite(&pRec->DEC, sizeof(float), 1, file);
s += fwrite(&pRec->AST, sizeof(float), 1, file);
s += fwrite(&pRec->stokes.I, sizeof(float), 1, file);
s += fwrite(&pRec->stokes.Q, sizeof(float), 1, file);
s += fwrite(&pRec->stokes.U, sizeof(float), 1, file);
s += fwrite(&pRec->stokes.V, sizeof(float), 1, file);

return s;
}
/**********************************************************************/

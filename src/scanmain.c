#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "common.h"
#include "fluxdata.h"

enum ScanTag {UPSCAN, DOWNSCAN, UPENDPOINT, DOWNENDPOINT, UPSCAN_START, UPSCAN_END, DOWNSCAN_START, DOWNSCAN_END};

static void print_usage(const char * prog)
{
	printf(
	"Usage: %s wapp<n> <lowchan> <highchan> <decmin> <decmax>\n"
	"\n", prog);
}

static void carve_scanlines(FluxWappData *wappdata, int lowchan, int highchan, float decmin, float decmax)
{
	int chan,i,j;
	FILE *scanfile;
	char filename[64+1];


	//for each channel
	for (chan=lowchan; chan<highchan; chan++)
	{
		//ChannelData chan_data;
		//chan_data.num_days = wappdata->numDays;
		//chan_data.chan = chan;

		printf("Channel: %i\n", chan);

		//read channel data files for all the days
		printf("reading channel data ...\n");
		fluxwappdata_readchan(wappdata, chan);
		
		printf ("seperating out scanlines for each day ...\n");

		//iterate over the days
		for (i=0; i<wappdata->numDays; i++)
		{
			FluxDayData * daydata = &wappdata->daydata[i];
			const int numRecords = daydata->numRecords;
			enum ScanTag *tags = malloc(sizeof(enum ScanTag) * numRecords);
			int upcount, downcount;
			int upcheck, downcheck;

			if (numRecords <= 0) {
				continue;
			}

			printf("day %s\n", daydata->mjd);

			//mark the data points with the direction tags
			for (j=0; j<numRecords-1; j++) 
			{
				if (daydata->records[j].DEC > daydata->records[j+1].DEC) {
					tags[j] = DOWNSCAN;
				} else if (daydata->records[j].DEC < daydata->records[j+1].DEC) {
					tags[j] = UPSCAN;
				}
			}

			//mark the endpoints at the limits
			for (j=0; j<numRecords; j++) 
			{
				if (daydata->records[j].DEC > decmax) {
					tags[j] = UPENDPOINT;
				} else if (daydata->records[j].DEC < decmin) {
					tags[j] = DOWNENDPOINT;
				}
			}

			//set a starting condition
			if (tags[0] == UPSCAN) {
			   tags[0] = UPSCAN_START;
			} else if (tags[0] == DOWNSCAN) {
				tags[0] = DOWNSCAN_START;
			}

			//delimit the scan tags 
			for (j=1; j<numRecords; j++) 
			{
				switch (tags[j]) {
					case UPSCAN:
						if (tags[j-1] == DOWNENDPOINT) {
							tags[j] = UPSCAN_START;
						} 
						break;
					case DOWNSCAN:
						if (tags[j-1] == UPENDPOINT) {
							tags[j] = DOWNSCAN_START;
						} 
						break;
					case UPENDPOINT:
						if (tags[j-1] == UPSCAN) {
							tags[j-1] = UPSCAN_END;
						} 
						break;
					case DOWNENDPOINT:
						if (tags[j-1] == DOWNSCAN) {
							tags[j-1] = DOWNSCAN_END;
						} 
						break;
					default:
						break;
				}
			}

			//set an ending condition
			if (tags[numRecords-1] == UPSCAN) {
				tags[numRecords-1] = UPSCAN_END;
			} else if (tags[numRecords-1] == DOWNSCAN) {
				tags[numRecords-1] = DOWNSCAN_END;
			}

			// sanity check the scan terminations
			upcount = downcount = 0;
			upcheck = downcheck = 0;
			for (j=0; j<numRecords; j++) 
			{
				switch (tags[j]) {
					case UPSCAN_START:
						upcount++;
						upcheck++;
						break;
					case DOWNSCAN_START:
						downcount++;
						downcheck++;
						break;
					case UPSCAN_END:
						upcheck--;
						break;
					case DOWNSCAN_END:
						downcheck--;
						break;
					default:
						break;
				}
				if (upcheck<0 || upcheck>1 || downcheck<0 || downcheck>1) {
					printf("WARN: Scan delimiters are incorrect (up:%i down:%i)\n", upcheck, downcheck);
				}
			}



			sprintf(filename, "%s/%s/chan_%03i", 
							daydata->mjd, wappdata->wapp, chan);
			mkdir(filename, S_IRWXU|S_IRWXG|S_IRWXO);

			
			scanfile = NULL;
			for (j=0; j<numRecords; j++) 
			{
				//skip over endpoints
				if (tags[j] == UPENDPOINT || tags[j] == DOWNENDPOINT) {
					continue;
				}

				//create a new file if needed
				if (tags[j] == UPSCAN_START) {
					sprintf(filename, "%s/%s/chan_%03i/upscan_%02i.dat", 
							daydata->mjd, wappdata->wapp, chan, upcount);
					scanfile = fopen(filename, "w");
				} else if (tags[j] == DOWNSCAN_START) {
					sprintf(filename, "%s/%s/chan_%03i/downscan_%02i.dat", 
							daydata->mjd, wappdata->wapp, chan, downcount);
					scanfile = fopen(filename, "w");
				}

				//write the data to file
				if (scanfile != NULL) {
					fluxrecord_write(&daydata->records[j], scanfile);
				} else {
					printf("ERROR: %s is not open\n", filename);
				}
				

				//close the file if needed
				if (tags[j] == UPSCAN_END) {
					fclose(scanfile);
					scanfile = NULL;
					upcount++;
				} else if (tags[j] == DOWNSCAN_END) {
					fclose(scanfile);
					scanfile = NULL;
					downcount++;
				}

			}

			free(tags);
		}

	}
}


int main (int argc, char* argv[])
{
	int numDays;
	char ** files;
	FluxWappData * wappdata;

	char * wapp;
	int lowchan, highchan;
	float decmin, decmax;

	if (argc != 6) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	} else {
		wapp = argv[1];
		lowchan = atoi(argv[2]);
		highchan = atoi(argv[3]);
		decmin = (float) atof(argv[4]);
		decmax = (float) atof(argv[5]);
	}
		
	numDays = get_date_dirs("./", &files);
	if (numDays <= 0) {
		printf("ERROR: could not find any date dirs\n");
		return EXIT_FAILURE;
	}

	// allocate and initialize the wapp data days
	wappdata = fluxwappdata_alloc(wapp, files, numDays);

	// do the seperation of the data into scanlines
	carve_scanlines(wappdata, lowchan, highchan, decmin, decmax);

	//free the data
	fluxwappdata_free(wappdata);

	return EXIT_SUCCESS;
}


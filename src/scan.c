#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "common.h"
#include "fluxdata.h"
#include "scan.h"

enum ScanTag {UPSCAN, DOWNSCAN, UPENDPOINT, DOWNENDPOINT};

/*
   For each element in dataset, tags the element with one of the values
   of the ScanTag enumeration.  The tag determination is made by checking the
   DEC field of the datapoint to determine if it is strictly increasing (UPSCAN),
   decreasing (DOWNSCAN) or beyond decmin and decmax (DOWNENDPOINT, UPENDPOINT).
   
   Its too bad this has to use decmin,decmax to determine what is in an endpoint.
   The problem is that the endpoints have changes in them that are not strictly
   decreasing/increasing.
 */
static void tag_scanlines(enum ScanTag *tags, FluxRecord dataset[], int numRecords, float decmin, float decmax)
{
	int j;

	//mark the data points with the direction tags
	for (j=0; j<numRecords-1; j++) 
	{
		if (dataset[j].DEC > (dataset[j+1].DEC-0.001)) { //decreasing
			tags[j] = DOWNSCAN;
		} else if (dataset[j].DEC < (dataset[j+1].DEC+0.001)) { //increasing
			tags[j] = UPSCAN;
		}
	}

	//mark the endpoints at the limits
	for (j=0; j<numRecords; j++) 
	{
		if (dataset[j].DEC > decmax) {
			tags[j] = UPENDPOINT;
		} else if (dataset[j].DEC < decmin) {
			tags[j] = DOWNENDPOINT;
		}
	}

	//ensure the last scan is correctly terminated with an endpoint
//	if (tags[numRecords-1] == DOWNSCAN) 
//		tags[numRecords-1] = DOWNENDPOINT;
//	if (tags[numRecords-1] == UPSCAN) 
//		tags[numRecords-1] = UPENDPOINT;
}

static void tag_scanlines2(enum ScanTag *tags, FluxRecord dataset[], int numRecords, float decmin, float decmax)
{
	int j;

	if (dataset[0].DEC > dataset[1].DEC) {
		tags[j] = DOWNSCAN;
	} else {
		tags[j] = UPSCAN;
	}

	//mark the data points with the direction tags
	for (j=1; j<numRecords-1; j++) 
	{
		float thisDEC = dataset[j].DEC;
		float nextDEC = dataset[j+1].DEC;
		if (tags[j-1] == DOWNSCAN) 
		{
			if (thisDEC > (nextDEC-0.001)) { //decreasing
				tags[j] = DOWNSCAN;
			} else if (thisDEC < (nextDEC+0.001) && thisDEC < decmin) { //increasing
				tags[j] = UPSCAN;
			}
		}
		if (tags[j-1] == UPSCAN) 
		{
			if (thisDEC > (nextDEC-0.001) && thisDEC > decmax) { //decreasing
				tags[j] = DOWNSCAN;
			} else if (thisDEC < (nextDEC+0.001)) { //increasing
				tags[j] = UPSCAN;
			}
		}
	}

	//mark the endpoints at the limits
	for (j=1; j<numRecords; j++) 
	{
		if (tags[j-1] == UPSCAN && tags[j] == DOWNSCAN) {
			tags[j] = UPENDPOINT;
		}
		if (tags[j-1] == DOWNSCAN && tags[j] == UPSCAN) {
			tags[j] = DOWNENDPOINT;
		}
	}

}



static int count_scanlines(enum ScanTag tags[], int size)
{
	int i;
	int upcount = 0;
	int downcount = 0;

	for (i=0; i<size-1; i++) 
	{
		if (tags[i] == DOWNSCAN && tags[i+1] != DOWNSCAN) {
			downcount++;
		}
		if (tags[i] == UPSCAN && tags[i+1] != UPSCAN) {
			upcount++;
		}
	}
	
	return upcount + downcount + 1;
}

static void output_tags(enum ScanTag tags[], FluxRecord dataset[], int size, char* day)
{
	int i;
	char filename[32];
	FILE * outfile;

	sprintf(filename, "tags_%s.ann", day);
	outfile = fopen(filename, "w");
	for (i=0; i<size; i++) {
		switch (tags[i]) {
			case UPSCAN:
				fprintf(outfile, "COLOUR BLUE\n");
				fprintf(outfile, "DOT %f %f\n", dataset[i].RA, dataset[i].DEC);
				break;
			case DOWNSCAN:
				fprintf(outfile, "COLOUR YELLOW\n");
				fprintf(outfile, "DOT %f %f\n", dataset[i].RA, dataset[i].DEC);
				break;
			default:
				fprintf(outfile, "COLOUR RED\n");
				fprintf(outfile, "DOT %f %f\n", dataset[i].RA, dataset[i].DEC);
				fprintf(outfile, "CIRCLE %f %f 0.005\n", dataset[i].RA, dataset[i].DEC);
				break;
		}
	}
}



void determine_scan_lines(FluxWappData * wappdata, float decmin, float decmax)
{
	int d;

	for (d=0; d<wappdata->numDays; d++) 
	{
		FluxDayData * daydata;
		ScanDayData * scanDayData;
		enum ScanTag *tags;
		enum ScanTag direction;
		int i;
		int numRecords;
		int numScans;
		int scanCount;

		daydata = &wappdata->daydata[d];
		scanDayData = &wappdata->scanDayData[d];
//		printf("day %s\n", daydata->mjd);
		numRecords = daydata->numRecords;

		tags = malloc(sizeof(enum ScanTag) * numRecords);
		tag_scanlines(tags, daydata->records, numRecords, decmin, decmax);

		numScans = count_scanlines(tags, numRecords);
		scanDayData->numScans = numScans;
		scanDayData->scans = malloc(sizeof(ScanData) * numScans);

		//assuming we always start with a downscan
//SSG		direction = DOWNSCAN;
		//SSG
		if(tags[0] == DOWNSCAN || tags[0] == UPENDPOINT)
			direction = DOWNSCAN;
		if(tags[0] == UPSCAN || tags[0] == DOWNENDPOINT)
			direction = UPSCAN;
		//SSG
		i = 0;
		scanCount = 0;
		while (i < numRecords) 
		{
			int start, end;

			//set start and end flags for this scan
			do {
				start = i;
				i++;
			} while (i<numRecords && tags[i] != direction);

			do {
				end = i;
				i++;
			} while (i<=numRecords && tags[i-1] == direction);

			if (end-start > 10) { //TODO confirm the arbitrairy limit
				scanDayData->scans[scanCount].records = &daydata->records[start];
				scanDayData->scans[scanCount].num_records = end-start;
				scanCount++;
			} else {
//				printf("WARN: short scan being ignored start:%i end:%i\n", start, end);
				;
			}
//SSG			direction = (direction == DOWNSCAN) ? UPSCAN : DOWNSCAN; 
//SSG
			if(tags[i] == DOWNSCAN || tags[i] == UPENDPOINT)
				direction = DOWNSCAN;
			if(tags[i] == UPSCAN || tags[i] == DOWNENDPOINT)
				direction = UPSCAN;
//SSG
	
		}
		//SSG
		if(scanCount < numScans)
		{		
//		printf("DIAGNOSTIC: scanCount %d numScans %d\n",scanCount-1,numScans);
		//SSG
		scanDayData->numScans = scanCount;//SSG
		}
//		printf("DIAGNOSTIC: new numscans %d\n",scanDayData->numScans);//SSG
		free(tags);
	}
}


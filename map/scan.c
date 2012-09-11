#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "common.h"
#include "fluxdata.h"
#include "scan.h"

enum ScanTag {UPSCAN, DOWNSCAN, UPENDPOINT, DOWNENDPOINT};
//-----------------------------------------------------------------------
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

for(j=0; j<numRecords-1; j++) 
	{
	if(dataset[j].DEC > (dataset[j+1].DEC-0.001)) tags[j] = DOWNSCAN;
	else if(dataset[j].DEC < (dataset[j+1].DEC+0.001))	tags[j] = UPSCAN;
	}
for(j=0; j<numRecords; j++) 
	{
	if(dataset[j].DEC > decmax) tags[j] = UPENDPOINT;
	else if(dataset[j].DEC < decmin) tags[j] = DOWNENDPOINT;
	}
//ensure the last scan is correctly terminated with an endpoint
//	if (tags[numRecords-1] == DOWNSCAN) 
//		tags[numRecords-1] = DOWNENDPOINT;
//	if (tags[numRecords-1] == UPSCAN) 
//		tags[numRecords-1] = UPENDPOINT;
}
//-----------------------------------------------------------------------
static int count_scanlines(enum ScanTag tags[], int size)
{
int i;
int upcount = 0;
int downcount = 0;

for(i=0; i<size-1; i++) 
	{
	if(tags[i] == DOWNSCAN && tags[i+1] != DOWNSCAN) downcount++;
	if(tags[i] == UPSCAN && tags[i+1] != UPSCAN) upcount++;
	}	
return upcount + downcount + 1;
}
//-----------------------------------------------------------------------
void determine_scan_lines(FluxWappData * wappdata, float decmin, float decmax)
{
int d;

for(d=0; d<wappdata->numDays; d++) 
	{
	FluxDayData *daydata;
	ScanDayData *scanDayData;
	enum ScanTag *tags, direction;
	int i, numRecords, numScans, scanCount;
	daydata = &wappdata->daydata[d];
	scanDayData = &wappdata->scanDayData[d];
	numRecords = daydata->numRecords;
	tags = malloc(sizeof(enum ScanTag)*numRecords);
	tag_scanlines(tags, daydata->records, numRecords, decmin, decmax);
	numScans = count_scanlines(tags, numRecords);
	scanDayData->numScans = numScans;
	scanDayData->scans = malloc(sizeof(ScanData)*numScans);
	if(tags[0] == DOWNSCAN || tags[0] == UPENDPOINT) direction = DOWNSCAN;		
	if(tags[0] == UPSCAN || tags[0] == DOWNENDPOINT) direction = UPSCAN;
	i = 0;
	scanCount = 0;
	while(i < numRecords) 
		{
		int start, end;
		//set start and end flags for this scan
		do{start = i; i++;} while(i<numRecords && tags[i] != direction);
		do{end = i;	i++;} while(i<=numRecords && tags[i-1] == direction);
		if(end - start > 10) 
			{
			scanDayData->scans[scanCount].records = &daydata->records[start];
			scanDayData->scans[scanCount].num_records = end - start;
			scanCount++;
			} 
		if(tags[i] == DOWNSCAN || tags[i] == UPENDPOINT) direction = DOWNSCAN;
		if(tags[i] == UPSCAN || tags[i] == DOWNENDPOINT) direction = UPSCAN;
		}		
	if(scanCount < numScans) scanDayData->numScans = scanCount; 
	free(tags);
	}
}
//-----------------------------------------------------------------------

/*
These data structures deal with the Stokes flux for a particular
channel.  

FluxRecord - a single line in a fluxtime.dat file
FluxDayData - a collection of flux records for a particular day
FluxWappData - an entire wapp worth of data

*/

#ifndef _FLUXDATA_H
#define _FLUXDATA_H

#include "common.h"
#include <stdio.h>
#include <stdlib.h>

#define MJD_LEN 15 
#define WAPP_LEN 15 

typedef struct {
    double I;
    double Q;
    double U;
    double V;
} StokesParams;

typedef struct {
        float RA;
        float DEC;
        float AST;
        StokesParams stokes;
} FluxRecord;
//struct _ScanData;

struct _CrossingPoint {
	float RA;
	struct _ScanData *crossScan;
	int ref_pos;
	int cross_pos;
};
typedef struct _CrossingPoint CrossingPoint;

struct _ScanData {
	int num_records;
	FluxRecord *records; //pointer to start of scan
	int num_cross_points;
	struct _CrossingPoint crossPoints[MAX_NUM_DAYS];
};
typedef struct _ScanData ScanData;

typedef struct {
	//char mjd[MJD_LEN+1];
	int numScans;
	ScanData *scans;
} ScanDayData;


typedef struct {
        char mjd[MJD_LEN+1];
        int numRecords;
        float RAmin;
        float RAmax;
        FluxRecord * records;
} FluxDayData;

typedef struct {
        char wapp[WAPP_LEN+1];
        int numDays;
        FluxDayData * daydata;
		ScanDayData * scanDayData; //new
} FluxWappData;


/*
typedef struct {
	char wapp[WAPP_LEN+1];
	int numDays;
	ScanDayData updaydata[MAX_NUM_DAYS];
	ScanDayData downdaydata[MAX_NUM_DAYS];
}ScanWappData;
*/
int fluxrecord_read(FluxRecord * pRec, FILE * file);
int fluxrecord_write(FluxRecord * pRec, FILE * file);
int fluxwappdata_readavg(FluxWappData * wappdata);
int fluxwappdata_readchan(FluxWappData * wappdata, int chan);
int fluxwappdata_writeavg(FluxWappData * wappdata);
int fluxwappdata_writechan(FluxWappData * wappdata, int chan);
FluxWappData * fluxwappdata_alloc(const char *wapp, char **days, int numDays);
void fluxwappdata_free(FluxWappData * wappdata);

int fluxchanneldata_readchan(FluxWappData * wappdata, int chan);
//int scanwappdata_readchan(ScanWappData *wappdata, char * wapp, int chan, char ** days, int num_days);



#endif //_FLUXDATA_H


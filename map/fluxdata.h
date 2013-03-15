/*
These data structures deal with the Stokes flux for a particular
channel.  

FluxRecord - a single line in a fluxtime.dat file
FluxDayData - a collection of flux records for a particular day
FluxWappData - an entire wapp worth of data for all days (single channel)
*/

#ifndef _FLUXDATA_H
#define _FLUXDATA_H

#include "common.h"
#include <stdio.h>
#include <stdlib.h>

#define MJD_LEN 15 
#define WAPP_LEN 15 

//sguram
#define BASKETWEAVE 0 //added to distinguish between which files to read for map and clean programs
#define CLEAN 1

typedef struct {
    float I;
    float Q;
    float U;
    float V;
} StokesParams;

typedef struct {
	float RA;
	float DEC;
	float AST;
	StokesParams stokes;
	int count;
	float weight;
} FluxRecord;

struct _CrossingPoint {
	float RA;
	float DEC;
	struct _ScanData *crossScan;
	int ref_pos;
	int cross_pos;
	StokesParams diff;
};
typedef struct _CrossingPoint CrossingPoint;

struct _ScanData {
	int scantype;
	int num_records;
	FluxRecord *records; //pointer to start of scan
	int num_cross_points;
	struct _CrossingPoint crossPoints[MAX_NUM_DAYS];
};
typedef struct _ScanData ScanData;

typedef struct {
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
	ScanDayData * scanDayData;
} FluxWappData;


int fluxrecord_read_binary(FluxRecord * pRec, FILE * file);
int fluxdaydata_read_binary(const char *field, FluxDayData *daydata, FILE *infile, int beam, int chan, int day);
void fluxwappdata_readchan_binary(const char *field, int band, FluxWappData * wappdata, int chan, int id, int avg, float decmin, float decmax); // band 0
//void fluxwappdata_readchan_binary1(const char *field, FluxWappData * wappdata, int chan, int id, int avg, float decmin, float decmax); // band 1
FluxWappData * fluxwappdata_alloc(const char *wapp, char **days, int numDays);
void fluxwappdata_free(FluxWappData * wappdata);
int fluxwappdata_writechan(FluxWappData * wappdata, int chan);
int fluxwappdata_writechan_binary(FluxWappData * wappdata, int chan);
int fluxrecord_write(FluxRecord * pRec, FILE * file);
int fluxrecord_write_binary(FluxRecord * pRec, FILE * file);
int fluxdaydata_read_binary_single_file(const char *field, FluxDayData *daydata, FILE *infile, FILE *configfile, int beam, int chan, int day);

#endif //_FLUXDATA_H


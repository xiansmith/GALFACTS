#ifndef _RFI_H
#define _RFI_H

#include "common.h"
#include "spec.h"


typedef struct {
	double OffXX[MAX_CHANNELS]; 
	double OffXY[MAX_CHANNELS]; 
	double OffYX[MAX_CHANNELS]; 
	double OffYY[MAX_CHANNELS]; 
	double OnXX[MAX_CHANNELS]; 
	double OnXY[MAX_CHANNELS]; 
	double OnYX[MAX_CHANNELS]; 
	double OnYY[MAX_CHANNELS]; 
} PolDifferences;

typedef struct {
	double meanOffXX;
	double meanOffYY;
	double meanOnXX;
	double meanOnYY;
	double sigmaOffXX;
	double sigmaOffYY;
	double sigmaOnXX;
	double sigmaOnYY;
	double meanOffXY;
	double meanOffYX;
	double meanOnXY;
	double meanOnYX;
	double sigmaOffXY;
	double sigmaOffYX;
	double sigmaOnXY;
	double sigmaOnYX;
} PolStatistics;


void compute_stats_on_diffs(const SpecRecord * pRec, const PolDifferences * pDiff, PolStatistics * pStat, int lowchan, int highchan);
void compute_stats_on_stats(const PolStatistics * stats, int size, PolStatistics * pStatOut);
void compute_diffs(const SpecRecord * pRec, PolDifferences * pDiff, int lowchan, int highchan);
void print_means(FILE * file, const SpecRecord * pRec, const PolStatistics * pStat);
void print_sigmas(FILE * file, const SpecRecord * pRec, const PolStatistics * pStat);
void print_diffs(FILE * file, const SpecRecord * pRec, const PolDifferences * pDiff, int lowchan, int highchan);
void rfi_detection(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float sigmaThreshold, 
	int ignoreA_low, int ignoreA_high, int ignoreB_low, int ignoreB_high);
void rfi_spanning(SpecRecord dataset[], int size, int lowchan, int highchan, int span);
void rfi_blanking(SpecRecord dataset[], int numRecords, int lowchan, int highchan, int rfiTolerance);
void rfi_write(SpecRecord dataset[], int size, int lowchan, int highchan, float * freq);


void aerostat_rfi_blanking(SpecRecord dataset[], int size, int lowchan, int highchan);

#endif


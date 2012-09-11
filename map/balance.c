#include "balance.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsLib.h"
#include "grid.h"
#include "string.h"
#include "chebyshev.h"
//-------------------------------------------------------------------------------
static int line_intersect(float x1,float y1,float x2,float y2, float x3,float y3,float x4,float y4, float *xi, float *yi)
{
/*
   Determine if the line segment drawn from the given points intersect and return
   their intersection point (xi, yi).

   Theory:
   Say there are two points, one on each line:
   Pa = P1 + ua * (P2-P1) 
   Pb = P3 + ub * (P4-P3)
   Solve for where Pa = Pb gives two equations with two unknowns
   x1 + ua*(x2-x1) = x3 + ub*(x4 - x3)
   y1 + ua*(y2-y1) = y3 + ub*(y4 - y3)
   The code below solves this system for ua and ub.  Simple line equation is used to
   compute RA and DEC positions. The check for 0 > ua < 1 (and same for ub)
   is required to determine if the line segments overlap, not just the lines.
*/
float denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1);
float ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom; //unknown a
float ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom; //unknown b
//assumed that if the lines share a point, they do not intersect
if(isnan(ua) || isnan(ub)) return 0;
if(ua < 0 || ua > 1 || ub < 0 || ub > 1) return 0;
*xi = x1 + ua*(x2-x1); *yi = y1 + ua*(y2-y1);
return 1;
}
//-------------------------------------------------------------------------------------------------------------
static float swapf_tmp;
#define SWAPF(X,Y) swapf_tmp=*X; *X=*Y; *Y=swapf_tmp;
//-------------------------------------------------------------------------------------------------------------
static int scan_intersect(ScanData *ref, ScanData *curr, float *pRA, float *pDEC, int *refpos, int *currpos)
{
	int i, j;
	float x1, x2, x3, x4, y1, y2, y3, y4;
	float refMin, refMax, currMin, currMax;

	refMax = ref->records[0].RA;
	refMin = ref->records[ref->num_records-1].RA;
	if(refMax < refMin) SWAPF(&refMax, &refMin);
	currMax = curr->records[0].RA;
	currMin = curr->records[curr->num_records-1].RA;
	if(currMax < currMin) SWAPF(&currMax, &currMin);
	//quick check
	if(currMin > refMax || currMax < refMin) return 0;
	//exhaustive determination
	i = j = 1;
	do{
		if(fabs(ref->records[i+1].RA - curr->records[j].RA) < fabs(ref->records[i].RA - curr->records[j].RA)) i++;
		else if(fabs(curr->records[j+1].RA - ref->records[i].RA) < fabs(curr->records[j].RA - ref->records[i].RA)) j++;
		else if(fabs(ref->records[i+1].DEC - curr->records[j].DEC) < fabs(ref->records[i].DEC - curr->records[j].DEC)) i++;
		else if(fabs(curr->records[j+1].DEC - ref->records[i].DEC) < fabs(curr->records[j].DEC - ref->records[i].DEC)) j++;
		else break;
		}while((i+1 < ref->num_records) && (j+1 < curr->num_records));
	//so, now we have the two points that are on the near side of crossing.
	//make a line segment with the next point to make a line segment to determine
	//the acutal intersection point
	x1 = ref->records[i-1].RA; 
	y1 = ref->records[i-1].DEC; 
	x2 = ref->records[i+1].RA; 
	y2 = ref->records[i+1].DEC; 
	x3 = curr->records[j-1].RA; 
	y3 = curr->records[j-1].DEC; 
	x4 = curr->records[j+1].RA; 
	y4 = curr->records[j+1].DEC; 
	//determine intersection location
	if(!line_intersect(x1, y1, x2, y2, x3, y3, x4, y4, pRA, pDEC)) return 0;
	//out of range check; TODO: is this useful?
	if(*pRA > refMax || *pRA < refMin) return 0;
	//line length check to prevent interpolation across excessive distances (1 arc min)
	//TODO: could have an max_interpolation_width_in_arcmin and divide by 60 to use below
	if(sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) > 0.0167 || sqrt((x4-x3)*(x4-x3) + (y4-y3)*(y4-y3)) > 0.0167 ) return 0;
	//i and j are not guarnteed to be the closest to the actual crossing position.
	//make a better final determination by checking the one on either side
	if(fabs(*pRA - ref->records[i+1].RA) < fabs(*pRA - ref->records[i].RA)) *refpos = i+1;
	else if(fabs(*pRA - ref->records[i-1].RA) < fabs(*pRA - ref->records[i].RA)) *refpos = i-1;
	else *refpos = i;
	if(fabs(*pRA - curr->records[j+1].RA) < fabs(*pRA - curr->records[j].RA)) *currpos = j+1;
	else if(fabs(*pRA - curr->records[j-1].RA) < fabs(*pRA - curr->records[j].RA)) *currpos = j-1;
	else *currpos = j;

	return 1;
}
//-------------------------------------------------------------------------------------------------------------
void find_intersections(FluxWappData *wappdata)
{
int r, c;
int i, j, k;
int numDays;

numDays = wappdata->numDays;
for(r=0; r<numDays; r++) 
	{
	ScanDayData *refScanDay = &wappdata->scanDayData[r];
	for(i=0; i<refScanDay->numScans; i++) 
		{
		ScanData *refscan = &refScanDay->scans[i];
		refscan->num_cross_points = 0;
		if(refscan->num_records == 0) continue;
		for(c=0; c<numDays; c++) 
			{
			ScanDayData *currScanDay = &wappdata->scanDayData[c];
			if(r==c) continue;
			for(j=0; j<currScanDay->numScans; j++) 
				{
				float RA, DEC;
				int refpos, crosspos;
				ScanData *currscan = &currScanDay->scans[j];
				if(currscan->num_records <= 0) continue;
				if(scan_intersect(refscan, currscan, &RA, &DEC, &refpos, &crosspos)) 
					{
					CrossingPoint *crossPoint;
					if(refscan->num_cross_points > MAX_NUM_DAYS) 
						{
						printf("WARN: maximum number of crossing points breached\n");
						break;
						}
					if(refscan->records[refpos].stokes.I == 0.0 || currscan->records[crosspos].stokes.I == 0.0) break;
					crossPoint = &refscan->crossPoints[refscan->num_cross_points];
					crossPoint->RA = RA;
					crossPoint->DEC = DEC;
					crossPoint->crossScan = currscan;
					crossPoint->ref_pos = refpos;
					crossPoint->cross_pos = crosspos;
					refscan->num_cross_points++;
					}
				}
			}
		//the scans are not in a useful order, so lets sort them by RA using simple bubble sort

		for(j=0; j<refscan->num_cross_points; j++) 
			{
			int minpos = j;
			for(k=j+1; k<refscan->num_cross_points; k++) 
				{
				if(refscan->crossPoints[k].RA < refscan->crossPoints[minpos].RA) minpos = k;
				}
			CrossingPoint tmp = refscan->crossPoints[j];
			refscan->crossPoints[j] =  refscan->crossPoints[minpos];
			refscan->crossPoints[minpos] = tmp;
			}
			
		}
	}
}
//-------------------------------------------------------------------------------------------------------------
static float scan_weave(ScanDayData *daydata, int scan, float loop_gain, int apply, int order)
{
int x, i, j, k, p, num_delta = 0;
float RA, nRA, min, max, delta_sum=0.0, nsigma = 2.5;
float dI[MAX_NUM_DAYS*MAX_NUM_SCANS], dQ[MAX_NUM_DAYS*MAX_NUM_SCANS], dU[MAX_NUM_DAYS*MAX_NUM_SCANS], dV[MAX_NUM_DAYS*MAX_NUM_SCANS], dRA[MAX_NUM_DAYS*MAX_NUM_SCANS];
float cI[order+1], cQ[order+1], cU[order+1], cV[order+1];
float *xminmax; xminmax = (float*)malloc(2*sizeof(float));
ScanData * refscan;

refscan = &daydata->scans[scan];
for(j=0; j<refscan->num_cross_points; j++) 
	{
	float refI, refQ, refU, refV;
	float crossI, crossQ, crossU, crossV;
	CrossingPoint *crossPoint =  &refscan->crossPoints[j];
	ScanData *crossScan = crossPoint->crossScan;
	RA = crossPoint->RA;
	p = crossPoint->ref_pos;
	if(p >= 0 && p < refscan->num_records) 
		{
		if(isfinite(refscan->records[p].stokes.I)) 
			{
			refI = refscan->records[p].stokes.I;
			refQ = refscan->records[p].stokes.Q;
			refU = refscan->records[p].stokes.U;
			refV = refscan->records[p].stokes.V;
			}
		}
	if(!isfinite(refI)) continue;
	p = crossPoint->cross_pos;
	if(p >= 0 && p < crossScan->num_records) 
		{
		if(isfinite(crossScan->records[p].stokes.I)) 
			{
			crossI = crossScan->records[p].stokes.I;
			crossQ = crossScan->records[p].stokes.Q;
			crossU = crossScan->records[p].stokes.U;
			crossV = crossScan->records[p].stokes.V;
			}
		}
	if(!isfinite(crossI)) continue;
	dI[num_delta] = refI - crossI;
	dQ[num_delta] = refQ - crossQ;
	dU[num_delta] = refU - crossU;
	dV[num_delta] = refV - crossV;
	dRA[num_delta] = RA;
	num_delta++;
	}

if(num_delta > order && apply)
	{
	delta_sum = 0.0; for(k=0;k<num_delta;k++) delta_sum += dI[k]*dI[k]; delta_sum = sqrt(delta_sum/num_delta);
	refscan = &daydata->scans[scan];
	chebyshev_minmax(dRA, num_delta, &min, &max);
	chebyshev_normalize(dRA, num_delta, min, max);
	chebyshev_fit_bw_scan(dRA, dI, num_delta, nsigma, cI, order);
	chebyshev_fit_bw_scan(dRA, dQ, num_delta, nsigma, cQ, order);
	chebyshev_fit_bw_scan(dRA, dU, num_delta, nsigma, cU, order);
	chebyshev_fit_bw_scan(dRA, dV, num_delta, nsigma, cV, order);

	for(k=0; k<refscan->num_records; k++)
		{
		if(isfinite(refscan->records[k].stokes.I) && refscan->records[k].RA >= min && refscan->records[k].RA <= max)
			{
			RA = CNORMALIZE(refscan->records[k].RA, min, max);
			refscan->records[k].stokes.I -= chebyshev_eval(RA, cI, order) * loop_gain;
			refscan->records[k].stokes.Q -= chebyshev_eval(RA, cQ, order) * loop_gain;
			refscan->records[k].stokes.U -= chebyshev_eval(RA, cU, order) * loop_gain;
			refscan->records[k].stokes.V -= chebyshev_eval(RA, cV, order) * loop_gain;
			}
		}
	}
	
return sqrt(delta_sum);
}
//-------------------------------------------------------------------------------------------------------------
static float day_weave(ScanDayData *daydata, float loop_gain, int apply, int order)
{
int i, j, k, p, num_delta = 0;
float RA, nRA, min, max, delta_sum = 0.0, nsigma = 2.5;
float refI, refQ, refU, refV, crossI, crossQ, crossU, crossV; 
float dI[MAX_NUM_DAYS*MAX_NUM_SCANS], dQ[MAX_NUM_DAYS*MAX_NUM_SCANS], dU[MAX_NUM_DAYS*MAX_NUM_SCANS], dV[MAX_NUM_DAYS*MAX_NUM_SCANS], dRA[MAX_NUM_DAYS*MAX_NUM_SCANS];
float cI[order+1], cQ[order+1], cU[order+1], cV[order+1];
ScanData *refscan;

for(i=0; i<daydata->numScans; i++) 
	{
	refscan = &daydata->scans[i];
	int num_cross_points = refscan->num_cross_points;		
	for(j=0; j<num_cross_points; j++) 
		{
		CrossingPoint *crossPoint =  &refscan->crossPoints[j];
		ScanData *crossScan = crossPoint->crossScan;
		RA = crossPoint->RA;				
		p = crossPoint->ref_pos;
		if(p >= 0 && p < refscan->num_records) 
			{
			refI = refscan->records[p].stokes.I;
			refQ = refscan->records[p].stokes.Q;
			refU = refscan->records[p].stokes.U;
			refV = refscan->records[p].stokes.V;
			}
		if(!isfinite(refI)) continue;		
		p = crossPoint->cross_pos;
		if(p >= 0 && p < crossScan->num_records) 
			{
			crossI = crossScan->records[p].stokes.I;
			crossQ = crossScan->records[p].stokes.Q;
			crossU = crossScan->records[p].stokes.U;
			crossV = crossScan->records[p].stokes.V;
			}
		if(!isfinite(crossI)) continue;
		dI[num_delta] = refI - crossI;
		dQ[num_delta] = refQ - crossQ;
		dU[num_delta] = refU - crossU;
		dV[num_delta] = refV - crossV;
		dRA[num_delta] = RA;
		num_delta++;
		}
	}


if(num_delta > order && apply)
	{
	delta_sum = 0.0; for(i=0; i<num_delta; i++) delta_sum += dI[i]*dI[i]; delta_sum = sqrt(delta_sum/num_delta);
	chebyshev_minmax(dRA, num_delta, &min, &max);
	chebyshev_normalize(dRA, num_delta, min, max);
	chebyshev_fit_bw_day(dRA, dI, num_delta, nsigma, cI, order);
	chebyshev_fit_bw_day(dRA, dQ, num_delta, nsigma, cQ, order);
	chebyshev_fit_bw_day(dRA, dU, num_delta, nsigma, cU, order);
	chebyshev_fit_bw_day(dRA, dV, num_delta, nsigma, cV, order);

	for(i=0; i<daydata->numScans; i++) 
		{
		refscan = &daydata->scans[i];
		for(k=0; k<refscan->num_records; k++) 
			{
			if(isfinite(refscan->records[k].stokes.I) && refscan->records[k].RA >= min && refscan->records[k].RA <= max)
				{
				float RA = CNORMALIZE(refscan->records[k].RA, min, max);
				refscan->records[k].stokes.I -= chebyshev_eval(RA, cI, order) * loop_gain;
				refscan->records[k].stokes.Q -= chebyshev_eval(RA, cQ, order) * loop_gain;
				refscan->records[k].stokes.U -= chebyshev_eval(RA, cU, order) * loop_gain;
				refscan->records[k].stokes.V -= chebyshev_eval(RA, cV, order) * loop_gain;
				}
			}				
		}
	}
	
return sqrt(delta_sum);
}
//-------------------------------------------------------------------------------------------------------------
static void basket_weave(FluxWappData *wappdata, int day_order, int scan_order, float loop_gain, float loop_epsilon, int order)
{
int r, i, count, ord;
float chisqtmp, chisqglobal, chisqday, chisqglobalprev, globalchange;


//do day by day weaving
chisqglobalprev = INFINITY; count = 0; 
do{
	chisqglobal = 0.0;
	for(r=0; r<wappdata->numDays; r++) 
		{ 
		chisqtmp = day_weave(&wappdata->scanDayData[r], loop_gain, 1, order);
		chisqglobal += chisqtmp;
		}
	if(wappdata->numDays) chisqglobal /= wappdata->numDays;
	count++;
	globalchange = chisqglobalprev - chisqglobal;
	printf("Day iteration:%i global:%f change:%f prev:%f\n", count, chisqglobal, globalchange, chisqglobalprev);
	chisqglobalprev = chisqglobal;
	}while(globalchange > loop_epsilon && count < day_order);

//do the scan by scan weaving
chisqglobalprev = INFINITY; count = 0; 
do{
	chisqglobal = 0.0;
	for(r=0; r<wappdata->numDays; r++) 
		{ 
		ScanDayData *daydata = &wappdata->scanDayData[r];
		chisqday = 0.0;
		int scancount =0;
		for(i=0; i<daydata->numScans; i++)
			{
			chisqtmp = 0.0;
			if(daydata->numScans)
				{
				chisqtmp = scan_weave(daydata, i, loop_gain*0.7, 1, order);
				scancount++;
				}
			chisqday += chisqtmp;
			}
		if(scancount) chisqday /= scancount;
		chisqglobal += chisqday;
		}
	if(wappdata->numDays) chisqglobal /= wappdata->numDays;
	count++;
	globalchange = chisqglobalprev - chisqglobal;
	printf("Scan iteration:%i global:%f change:%f prev:%f\n", count, chisqglobal, globalchange, chisqglobalprev);
	chisqglobalprev = chisqglobal;
	}while(globalchange > loop_epsilon && count < scan_order);
}
//-------------------------------------------------------------------------------------------------------------
void balance_data(FluxWappData * wappdata, int day_order, int scan_order, float loop_gain, float loop_epsilon, int bw_order)
{
/* 
Main entry point to the balancing operation.
After this function returns, the entire frequency should be balanced day to day.
decmin and decmax set ranges where data outside these decs (in the endpoints) will be ignored.
loopgain - a multiplier that is used to limit the amount of change on each iteration.
loop_epsilon - the threshold where any percent changes less than this signal the termiation criteria of the iterations.
*/

printf("Performing basket weaving\n"); printf("Loop_gain: %f Loop_epsilon: %f\n", loop_gain, loop_epsilon);

basket_weave(wappdata, day_order, scan_order, loop_gain, loop_epsilon, bw_order);
}
//-------------------------------------------------------------------------------------------------------------

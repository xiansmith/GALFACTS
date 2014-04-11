#include "balance.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jsd/jsd_futil.h"
#include "jsd/jsd_util.h"
#include "jsd/jsd_fit.h"
#include "programs/fitsLib.h"
#include "grid.h"
#include "string.h"


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
static int line_intersect(double x1,double y1,double x2,double y2,
		double x3,double y3,double x4,double y4,
		double *xi, double *yi)
{
	double denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1);
	double ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom; //unknown a
	double ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom; //unknown b

	//assumed that if the lines share a point, they do not intersect
	if (isnan(ua) || isnan(ub))
		return 0;

	if (ua < 0 || ua > 1 || ub < 0 || ub > 1) {
		return 0;
	}

	*xi = x1 + ua*(x2-x1);
	*yi = y1 + ua*(y2-y1);

	return 1;
}


static float swapf_tmp;
#define SWAPF(X,Y) swapf_tmp=*X; *X=*Y; *Y=swapf_tmp;

static int scan_intersect(FILE * file, ScanData *ref, ScanData *curr, double *pRA, double *pDEC, int *refpos, int *currpos)
{
	int i, j;
	float x1, x2, x3, x4, y1, y2, y3, y4;
	float refMin, refMax, currMin, currMax;

	refMax = ref->records[0].RA;
	refMin = ref->records[ref->num_records-1].RA;
	if (refMax < refMin) SWAPF(&refMax, &refMin);

	currMax = curr->records[0].RA;
	currMin = curr->records[curr->num_records-1].RA;
	if (currMax < currMin) SWAPF(&currMax, &currMin);
	
	//quick check
	if (currMin > refMax || currMax < refMin) {
		return 0;
	}

	//exhaustive determination
	i = j = 1;
	do {
		if (fabs(ref->records[i+1].RA - curr->records[j].RA) < 
				fabs(ref->records[i].RA - curr->records[j].RA)) 
			i++;
		else if (fabs(curr->records[j+1].RA - ref->records[i].RA) < 
				fabs(curr->records[j].RA - ref->records[i].RA)) 
			j++;
		else if (fabs(ref->records[i+1].DEC - curr->records[j].DEC) < 
				fabs(ref->records[i].DEC - curr->records[j].DEC)) 
			i++;
		else if (fabs(curr->records[j+1].DEC - ref->records[i].DEC) < 
				fabs(curr->records[j].DEC - ref->records[i].DEC)) 
			j++;
		else 
			break;
	} while ((i+1 < ref->num_records) && (j+1 < curr->num_records));


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
	if (!line_intersect(x1, y1, x2, y2, x3, y3, x4, y4, pRA, pDEC)) {
		return 0;
	}
	

	//out of range check
	//TODO: is this useful?
	if (*pRA > refMax || *pRA < refMin) {
		return 0;
	}

	// line length check to prevent interpolation across excessive distances (1 arc min)
	//TODO: could have an max_interpolation_width_in_arcmin and divide by 60 to use below
	if (sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) > 0.0167 ||
		sqrt((x4-x3)*(x4-x3) + (y4-y3)*(y4-y3)) > 0.0167 ) {
		return 0;
	}

	//i and j are not guarnteed to be the closest to the actual crossing position.
	//make a better final determination by checking the one on either side
	if ( fabs(*pRA - ref->records[i+1].RA) < fabs(*pRA - ref->records[i].RA) ) *refpos = i+1;
	else if ( fabs(*pRA - ref->records[i-1].RA) < fabs(*pRA - ref->records[i].RA) ) *refpos = i-1;
	else *refpos = i;

	if ( fabs(*pRA - curr->records[j+1].RA) < fabs(*pRA - curr->records[j].RA) ) *currpos = j+1;
	else if ( fabs(*pRA - curr->records[j-1].RA) < fabs(*pRA - curr->records[j].RA) ) *currpos = j-1;
	else *currpos = j;

	//printf("found cross: %f %f\n", *RA, *DEC);
//	fprintf(file, "COLOUR BLUE\n"); //the cross lines
//	fprintf(file, "LINE W %f %f %f %f\n", x1, y1, x2, y2);
//	fprintf(file, "LINE W %f %f %f %f\n", x3, y3, x4, y4);
//	fprintf(file, "COLOUR GREEN\n"); //the intersection point
//	fprintf(file, "CROSS W %f %f 0.004 0.004\n", *pRA, *pDEC);
//	fprintf(file, "COLOUR YELLOW\n"); //two points that are closest to the intersection point
//	fprintf(file, "DOT W %f %f\n", curr->records[*currpos].RA, curr->records[*currpos].DEC );
//	fprintf(file, "DOT W %f %f\n", ref->records[*refpos].RA, ref->records[*refpos].DEC );

	return 1;
}


static void print_scandata(FILE * file, ScanData * scandata)
{
	int i;
	for (i=0; i<scandata->num_records; i++) {
		fprintf(file, "%f %f\n", scandata->records[i].RA,  scandata->records[i].DEC);
	}
}

static void find_intersections(FluxWappData *wappdata)
{
	int r, c;
	int i, j, k;
	int numDays;

	numDays = wappdata->numDays;

	for (r=0; r<numDays; r++) 
	{
		FILE * crossfile;
		char filename[32+1];

		FluxDayData *refday = &wappdata->daydata[r];
		ScanDayData *refScanDay = &wappdata->scanDayData[r];
		//printf("refday: %s\n", refday->mjd);

//		sprintf(filename, "cross%s.ann", refday->mjd);
//		sprintf(filename, "cross%s_beam%d.ann", refday->mjd, r%7); //ssg to fix for overwriting of files
//		crossfile = fopen(filename, "w");

		for (i=0; i<refScanDay->numScans; i++) 
		{
			ScanData *refscan = &refScanDay->scans[i];
			refscan->num_cross_points = 0;
//			fprintf(crossfile, "#scan %i\n", i);
			if (refscan->num_records == 0) continue;

			for (c=0; c<numDays; c++) 
			{
				//FluxDayData *currday = &wappdata->daydata[c];
				ScanDayData *currScanDay = &wappdata->scanDayData[c];
				if (r==c) continue;

				//fprintf(crossfile, "#%s\n", currday->mjd);

				for (j=0; j<currScanDay->numScans; j++) 
				{
					double RA, DEC;
					int refpos, crosspos;
					ScanData *currscan = &currScanDay->scans[j];
					if (currscan->num_records <= 0) continue;

					if (scan_intersect(crossfile, refscan, currscan, &RA, &DEC, &refpos, &crosspos)) 
					{
						CrossingPoint *crossPoint;

						if (refscan->num_cross_points > MAX_NUM_DAYS) {
							printf("WARN: maximum number of crossing points breached\n");
							break;
						}
						if (refscan->records[refpos].stokes.I == 0.0 
							|| currscan->records[crosspos].stokes.I == 0.0) {
							break;
						}


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
			for (j=0; j<refscan->num_cross_points; j++) {
				int minpos = j;
				for (k=j+1; k<refscan->num_cross_points; k++) {
					if (refscan->crossPoints[k].RA < refscan->crossPoints[minpos].RA) 
						minpos = k;
				}
				CrossingPoint tmp = refscan->crossPoints[j];
				refscan->crossPoints[j] =  refscan->crossPoints[minpos];
				refscan->crossPoints[minpos] = tmp;
			}
		}
//		fclose(crossfile);
	}
}


static void cross_line(FluxRecord records[], int size, int pos, double cI[], double cQ[], double cU[], double cV[], int order, double *xminmax)
{
	int k, p;
	//TODO: make these constants parameters
	const int stroke = 10; //number of points on either side
	const float nsigma = 3.0;
	double chisq;
	double yI[MAX_NUM_DAYS];
	double yQ[MAX_NUM_DAYS];
	double yU[MAX_NUM_DAYS];
	double yV[MAX_NUM_DAYS];
	double RA[MAX_NUM_DAYS];
	int count = 0;

	for (k=-stroke; k<=stroke; k++) {
		p = pos+k;
		if (p >= 0 && p < size) {
			if (isfinite(records[p].stokes.I)) {
				yI[count] = records[p].stokes.I;
				yQ[count] = records[p].stokes.Q;
				yU[count] = records[p].stokes.U;
				yV[count] = records[p].stokes.V;
				RA[count] = records[p].RA;
				count++;
			}
		}
	}
        double min, max;
        jsd_minmax(RA, count, &min, &max);
        xminmax[0] = min; xminmax[1] = max;
        jsd_normalize(RA, count, min, max);

	if (count > order) {
		jsd_poly_fit(RA, yI, count, nsigma, cI, order, &chisq);
		jsd_poly_fit(RA, yQ, count, nsigma, cQ, order, &chisq);
		jsd_poly_fit(RA, yU, count, nsigma, cU, order, &chisq);
		jsd_poly_fit(RA, yV, count, nsigma, cV, order, &chisq);
	} else {
		memset(cI, 0, sizeof(double));
		memset(cQ, 0, sizeof(double));
		memset(cU, 0, sizeof(double));
		memset(cV, 0, sizeof(double));
	}
	
}

//compute the average value of the array elements below low and high
//inclusive of array[low], exclusive of array[high]
static double arravg(double array[], int low, int high)
{
	double sum;
	int i, count;

	sum = 0.0;
	if (low < 0) low = 0;
	for (i=low; i<high; i++) {
		sum += array[i];
		count++;
	}
	return sum/count;
}


//TODO: make this a parameter as needed
#define CROSS_ORDER 2

static double scan_weave(ScanDayData *daydata, int scan, int order, float loop_gain, int apply)
{
	//TODO: make these constants parameters
	const float nsigma = 2.5;
	const int span = 0; //number of scans to include on either side of the ref scan

	int x, j;
	int k;
	double RA, DEC,nRA;
	double min, max;
	int num_delta;
	double delta_sum = 0.0;

	double dI[MAX_NUM_DAYS*MAX_NUM_SCANS], dQ[MAX_NUM_DAYS*MAX_NUM_SCANS], dU[MAX_NUM_DAYS*MAX_NUM_SCANS], dV[MAX_NUM_DAYS*MAX_NUM_SCANS];
	double dRA[MAX_NUM_DAYS*MAX_NUM_SCANS];
	double chisq[4];
	double cI[CROSS_ORDER+1], cQ[CROSS_ORDER+1], cU[CROSS_ORDER+1], cV[CROSS_ORDER+1];
	ScanData * refscan;
        double *xminmax; xminmax = (double*)malloc(2*sizeof(double));

	num_delta = 0;

	//do the curve fit to data across several scans (span*2+1)
	for (x=-span; x<=span; x++) 
	{
		int pos;

		if (scan+x < 0) continue;
		else if (scan+x >= daydata->numScans) continue;
		else pos = scan+x;

		refscan = &daydata->scans[pos];
		for (j=0; j<refscan->num_cross_points; j++) 
		{
			double refI, refQ, refU, refV;
			double crossI, crossQ, crossU, crossV;

			CrossingPoint *crossPoint =  &refscan->crossPoints[j];
			ScanData *crossScan = crossPoint->crossScan;

			//get RA and DEC from the crossing point structure
			RA = crossPoint->RA;
			DEC = crossPoint->DEC;

			//compute interpolated values for RA and DEC for the refscan
			cross_line(refscan->records, refscan->num_records, crossPoint->ref_pos, cI, cQ, cU, cV, CROSS_ORDER,xminmax);
                        nRA = NORMALIZE(RA, xminmax[0], xminmax[1]);
			refI = jsd_poly_eval(nRA, cI, CROSS_ORDER);
			refQ = jsd_poly_eval(nRA, cQ, CROSS_ORDER);
			refU = jsd_poly_eval(nRA, cU, CROSS_ORDER);
			refV = jsd_poly_eval(nRA, cV, CROSS_ORDER);
			if(!isfinite(refI)) continue;
			//compute interpolated values for RA and DEC for the crossscan
			cross_line(crossScan->records, crossScan->num_records, crossPoint->cross_pos, cI, cQ, cU, cV, CROSS_ORDER,xminmax);
                        nRA = NORMALIZE(RA, xminmax[0], xminmax[1]);
			crossI = jsd_poly_eval(nRA, cI, CROSS_ORDER);
			crossQ = jsd_poly_eval(nRA, cQ, CROSS_ORDER);
			crossU = jsd_poly_eval(nRA, cU, CROSS_ORDER);
			crossV = jsd_poly_eval(nRA, cV, CROSS_ORDER);
			if(!isfinite(crossI)) continue;

			dI[num_delta] = refI - crossI;
			dQ[num_delta] = refQ - crossQ;
			dU[num_delta] = refU - crossU;
			dV[num_delta] = refV - crossV;
			dRA[num_delta] = RA;
			num_delta++;
		}
	}
        free(xminmax);

	//ensure we have at least enough points
	if (num_delta > order) 
	{

		//constrain the endpoints
		int i;
		int num = 10; //number of points to constrain on each end.
		for (i=num_delta-num; i<num_delta; i++) {
			dI[i] = arravg(dI, i-num, i);
			dQ[i] = arravg(dQ, i-num, i);
			dU[i] = arravg(dU, i-num, i);
			dV[i] = arravg(dV, i-num, i);
		}
		for (i=num; i>=0; i--) {
			dI[i] = arravg(dI, i, i+num);
			dQ[i] = arravg(dQ, i, i+num);
			dU[i] = arravg(dU, i, i+num);
			dV[i] = arravg(dV, i, i+num);
		}

		//only apply the curve fit to the current scan
		//apply the curve fits
		if (apply) 
		{
			delta_sum = 0.0;
                        for(k=0;k<num_delta;k++)
                                delta_sum += dI[k]*dI[k];
                        delta_sum/=num_delta;

			//printf("deltasum:%lf num_delta:%d\n",delta_sum,num_delta);

			refscan = &daydata->scans[scan];

			jsd_minmax(dRA, num_delta, &min, &max);
			jsd_normalize(dRA, num_delta, min, max);

			//curve fit the dX values
			jsd_poly_fit(dRA, dI, num_delta, nsigma, cI, order, &chisq[0]);
			jsd_poly_fit(dRA, dQ, num_delta, nsigma, cQ, order, &chisq[1]);
			jsd_poly_fit(dRA, dU, num_delta, nsigma, cU, order, &chisq[2]);
			jsd_poly_fit(dRA, dV, num_delta, nsigma, cV, order, &chisq[3]);

			//jsd_print_poly(stdout, cI, order);

			for (k=0; k<refscan->num_records; k++) 
			{
				RA = NORMALIZE(refscan->records[k].RA, min, max);
				refscan->records[k].stokes.I -= jsd_poly_eval(RA, cI, order) * loop_gain;
				refscan->records[k].stokes.Q -= jsd_poly_eval(RA, cQ, order) * loop_gain;
				refscan->records[k].stokes.U -= jsd_poly_eval(RA, cU, order) * loop_gain;
				refscan->records[k].stokes.V -= jsd_poly_eval(RA, cV, order) * loop_gain;
			}

		}

		//sum up delta magnitudes and return result (for I only)!
/*		delta_sum = 0.0;
		for (k=0; k<num_delta; k++) {
			delta_sum += dI[k]*dI[k];
		}*/
		return delta_sum;
	}
	else
	{
		return 0.0;
	}
}

/*
static double mesh_weave(FluxWappData * wappdata)
{
	int h, i, j;
//	int k;
	double RA, DEC;
	double dI, dQ, dU, dV;
	FILE * weavefile;
	ScanDayData *daydata;
	double cI[CROSS_ORDER+1], cQ[CROSS_ORDER+1], cU[CROSS_ORDER+1], cV[CROSS_ORDER+1];

	weavefile = fopen("mesh.dat", "w");
	//fprintf(weavefile, "#RA DEC dI dQ dU dV\n");
	fprintf(weavefile, "#RA DEC dI\n");

	//for (h=0; h<wappdata->numDays; h++)
	for (h=0; h<wappdata->numDays; h+=7)
	{
		daydata = &wappdata->scanDayData[h];
		for (i=0; i<daydata->numScans; i++) 
		{
			ScanData *refscan = &daydata->scans[i];
			int num_cross_points = refscan->num_cross_points;

			for (j=0; j<num_cross_points; j++) 
			{

				double refI, refQ, refU, refV;
				double crossI, crossQ, crossU, crossV;
				CrossingPoint *crossPoint =  &refscan->crossPoints[j];
				ScanData *crossScan = crossPoint->crossScan;

				//get RA and DEC from the crossing point structure
				RA = crossPoint->RA;
				DEC = crossPoint->DEC;

				cross_line(refscan->records, refscan->num_records, crossPoint->ref_pos, cI, cQ, cU, cV, CROSS_ORDER);
				refI = jsd_poly_eval(RA, cI, CROSS_ORDER);
				refQ = jsd_poly_eval(RA, cQ, CROSS_ORDER);
				refU = jsd_poly_eval(RA, cU, CROSS_ORDER);
				refV = jsd_poly_eval(RA, cV, CROSS_ORDER);
				if (!isfinite(refI)) continue;

				cross_line(crossScan->records, crossScan->num_records, crossPoint->cross_pos, cI, cQ, cU, cV, CROSS_ORDER);
				crossI = jsd_poly_eval(RA, cI, CROSS_ORDER);
				crossQ = jsd_poly_eval(RA, cQ, CROSS_ORDER);
				crossU = jsd_poly_eval(RA, cU, CROSS_ORDER);
				crossV = jsd_poly_eval(RA, cV, CROSS_ORDER);
				if (!isfinite(crossI)) continue;

				dI = refI - crossI;
				dQ = refQ - crossQ;
				dU = refU - crossU;
				dV = refV - crossV;

				//fprintf(weavefile, "%g %g %g %g %g %g %g\n", RA, DEC, dI, dQ, dU, dV);
				fprintf(weavefile, "%lf %lf %lf %lf\n", RA, DEC, dI, 1.0);

			}
		}
	}

	fclose (weavefile);
}*/

static void apply_difference_corrections(ScanDayData *daydata, double *dRA, double *dI, double *dQ, double *dU, double *dV, int num_delta, int order, float loop_gain)
{
	int i, k;
	double min, max;
	const float nsigma = 2.5;
	double cI[MAX_NUM_DAYS], cQ[MAX_NUM_DAYS], cU[MAX_NUM_DAYS], cV[MAX_NUM_DAYS];
	double chisq;

	jsd_minmax(dRA, num_delta, &min, &max);
	jsd_normalize(dRA, num_delta, min, max);

	if(num_delta > order)
	{
		jsd_poly_fit(dRA, dI, num_delta, nsigma, cI, order, &chisq);
		jsd_poly_fit(dRA, dQ, num_delta, nsigma, cQ, order, &chisq);
		jsd_poly_fit(dRA, dU, num_delta, nsigma, cU, order, &chisq);
		jsd_poly_fit(dRA, dV, num_delta, nsigma, cV, order, &chisq);

		//jsd_print_poly(stdout, cI, order);
		//apply the dX values
		for (i=0; i<daydata->numScans; i++) 
		{
			ScanData *refscan = &daydata->scans[i];
			for (k=0; k<refscan->num_records; k++) 
			{
				double RA = NORMALIZE(refscan->records[k].RA, min, max);
				refscan->records[k].stokes.I -= jsd_poly_eval(RA, cI, order) * loop_gain;
				refscan->records[k].stokes.Q -= jsd_poly_eval(RA, cQ, order) * loop_gain;
				refscan->records[k].stokes.U -= jsd_poly_eval(RA, cU, order) * loop_gain;
				refscan->records[k].stokes.V -= jsd_poly_eval(RA, cV, order) * loop_gain;
			}
		}
	}
}

static double day_weave(ScanDayData *daydata, int order, float loop_gain, int apply)
{
	int i, j;
	double RA, DEC,nRA;
	double dI[MAX_NUM_DAYS*MAX_NUM_SCANS], dQ[MAX_NUM_DAYS*MAX_NUM_SCANS], dU[MAX_NUM_DAYS*MAX_NUM_SCANS], dV[MAX_NUM_DAYS*MAX_NUM_SCANS];
	//TODO: cX size is related to the order, not days
	double cI[MAX_NUM_DAYS], cQ[MAX_NUM_DAYS], cU[MAX_NUM_DAYS], cV[MAX_NUM_DAYS];
	double dRA[MAX_NUM_DAYS*MAX_NUM_SCANS];
        double *xminmax; xminmax = (double*)malloc(2*sizeof(double));
	int num_delta = 0;
	for (i=0; i<daydata->numScans; i++) 
	{
		ScanData *refscan = &daydata->scans[i];
		int num_cross_points = refscan->num_cross_points;

		for (j=0; j<num_cross_points; j++) 
		{

			double refI, refQ, refU, refV;
			double crossI, crossQ, crossU, crossV;
			CrossingPoint *crossPoint =  &refscan->crossPoints[j];
			ScanData *crossScan = crossPoint->crossScan;

			//get RA and DEC from the crossing point structure
			RA = crossPoint->RA;
			DEC = crossPoint->DEC;

			cross_line(refscan->records, refscan->num_records, crossPoint->ref_pos, cI, cQ, cU, cV, CROSS_ORDER,xminmax);
                        nRA = NORMALIZE(RA, xminmax[0], xminmax[1]);
			refI = jsd_poly_eval(nRA, cI, CROSS_ORDER);
			refQ = jsd_poly_eval(nRA, cQ, CROSS_ORDER);
			refU = jsd_poly_eval(nRA, cU, CROSS_ORDER);
			refV = jsd_poly_eval(nRA, cV, CROSS_ORDER);
			if (!isfinite(refI)) continue;

			cross_line(crossScan->records, crossScan->num_records, crossPoint->cross_pos, cI, cQ, cU, cV, CROSS_ORDER,xminmax);
                        nRA = NORMALIZE(RA, xminmax[0], xminmax[1]);
			crossI = jsd_poly_eval(nRA, cI, CROSS_ORDER);
			crossQ = jsd_poly_eval(nRA, cQ, CROSS_ORDER);
			crossU = jsd_poly_eval(nRA, cU, CROSS_ORDER);
			crossV = jsd_poly_eval(nRA, cV, CROSS_ORDER);
			if (!isfinite(crossI)) continue;

			dI[num_delta] = refI - crossI;
			dQ[num_delta] = refQ - crossQ;
			dU[num_delta] = refU - crossU;
			dV[num_delta] = refV - crossV;
			dRA[num_delta] = RA;
			num_delta++;
		}
	}

        free(xminmax);

	if (num_delta <= 0) {
		//printf("WARN: no crossing points for day %s\n", daydata->mjd);
		return 0.0;
	}

        double sum = 0.0;
        for (i=0; i<num_delta; i++) {
                sum += dI[i]*dI[i];
//              sum += dQ[i]*dQ[i];
        }


	if (apply) {
		apply_difference_corrections(daydata, dRA, dI, dQ, dU, dV, num_delta, order, loop_gain);
	}

	/* Difficult to determine what a good criteria for a good or bad
	day is.  chisq values are not valid, since bad data could have a good fit.
	using the average absolute difference instead.
	*/
	/*{
		double sum = 0.0;
		for (i=0; i<num_delta; i++) {
			sum += dI[i]*dI[i];
		}
		return sum/num_delta;
	}*/
	return sum/num_delta;
}

#define PERCENT_CHANGE(A,B) ((A-B)/B)

static void basket_weave(FluxWappData *wappdata, FILE * chisqfile, int day_order, int scan_order, float loop_gain, float loop_epsilon, MapMetaData * md, enum ShowProgress show_progress)
{
	int r, i;
	int count;
	double chisqtmp, chisqglobal, chisqday, chisqglobalprev;
	double globalchange;
	int ord;
	static header_param_list hpar;
	FILE *Iprogressfile,*Qprogressfile,*Uprogressfile,*Vprogressfile;
	float *dataI, *dataQ, *dataU, *dataV, *weight;
	int n1, n2, n3;
	if (show_progress) 
	{
		n1 = md->n1;
		n2 = md->n2;
		n3 = 0; //don't know how many planes in the cube, it will be set later below

		//start a fits cube
//		printf("Requesting malloc for %u bytes\n",n1 * n2 * sizeof (float));
		dataI = (float *) malloc (n1 * n2 * sizeof (float));
//		printf("Requesting malloc for %u bytes\n",n1 * n2 * sizeof (float));
		dataQ = (float *) malloc (n1 * n2 * sizeof (float));
//		printf("Requesting malloc for %u bytes\n",n1 * n2 * sizeof (float));
		dataU = (float *) malloc (n1 * n2 * sizeof (float));
//		printf("Requesting malloc for %u bytes\n",n1 * n2 * sizeof (float));
		dataV = (float *) malloc (n1 * n2 * sizeof (float));
//		printf("Requesting malloc for %u bytes\n",n1 * n2 * sizeof (float));
		weight = (float *) malloc (n1 * n2 * sizeof (float));

		init_header_param_list (&hpar);
		hpar.bitpix = -32;
		//hpar.num_axes = 2;
		hpar.num_axes = 3;
		hpar.naxis[0] = n1;
		hpar.naxis[1] = n2;
		hpar.naxis[2] = 0;
		sprintf (hpar.ctype[0], "RA---CAR");
		sprintf (hpar.ctype[1], "DEC--CAR");
		sprintf (hpar.ctype[2], "Iteration");
		hpar.crval[0] = md->RAcen;          /* hours */
		hpar.crval[1] = md->DECcen;               /* degrees */
		hpar.crval[2] = 0.0;
		hpar.crpix[0] = 0.5 + n1 / 2.0; /* image center in pixels */
		hpar.crpix[1] = 0.5 + n2 / 2.0;
		hpar.crpix[2] = 0.0;
		hpar.cdelt[0] = -md->cellsize;                     /* degrees */
		hpar.cdelt[1] = md->cellsize;                     /* degrees */
		hpar.cdelt[2] = 1.0;
		hpar.equinox = 2000.0;
		sprintf (hpar.bunit, "Kelvin");
		sprintf (hpar.telescope, "Arecibo");
		
		grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
		
		if(show_progress == ALL || show_progress == I)
		{
			sprintf (hpar.object, "I Basketweaving Progress");
			Iprogressfile = fopen("Iprogress.fits", "w");
			writefits_header(Iprogressfile, &hpar);
			writefits_plane(Iprogressfile, dataI, &hpar);
		//printf("plane %i\n", n3++);
		}
		if(show_progress == ALL || show_progress == Q)
		{		
			sprintf (hpar.object, "Q Basketweaving Progress");
			Qprogressfile = fopen("Qprogress.fits", "w");
			writefits_header(Qprogressfile, &hpar);
			writefits_plane(Qprogressfile, dataQ, &hpar);
		}
		if(show_progress == ALL || show_progress == U)
		{		
			sprintf (hpar.object, "U Basketweaving Progress");
			Uprogressfile = fopen("Uprogress.fits", "w");
			writefits_header(Uprogressfile, &hpar);
			writefits_plane(Uprogressfile, dataU, &hpar);
		}	
		if(show_progress == ALL || show_progress == V)
		{		
			sprintf (hpar.object, "V Basketweaving Progress");
			Vprogressfile = fopen("Vprogress.fits", "w");
			writefits_header(Vprogressfile, &hpar);
			writefits_plane(Vprogressfile, dataV, &hpar);
		}
	}

	//do day by day weaving
	chisqglobalprev = INFINITY;
	for (ord=0; ord<=day_order; ord++) 
	//ord = day_order;
	{
		count = 0;
		printf("Day weaving order %i\n", ord);
		do {
			chisqglobal = 0.0;
			for (r=0; r<wappdata->numDays; r++) 
			{ 
				chisqtmp = day_weave(&wappdata->scanDayData[r], ord, loop_gain, 1);
				chisqglobal += chisqtmp;
				//printf("weaved day:%i chisq:%f\n", r, chisqtmp);
			}
			chisqglobal /= wappdata->numDays;
			
			/*
			//done a weave to a set of days, write a progress plane
			if (show_progress) {
				printf("writing progress plane %i\n", n3++);
				grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
				writefits_plane(progressfile, dataQ, &hpar);
			}
			*/

			count++;
			globalchange = chisqglobalprev - chisqglobal;
			chisqglobalprev = chisqglobal;
			printf("Day iteration:%i global chisq:%g change:%g\n", count, chisqglobal, globalchange);
			fprintf(chisqfile, "%f %f\n", chisqglobal, globalchange);

		} while (fabs(globalchange) > loop_epsilon && count < 20);

		//done all weaves for this order, so write a progress plane
		if (show_progress) {
			printf("Writing progress plane %i\n", n3++);
			grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
			if(show_progress == ALL || show_progress == I)
				writefits_plane(Iprogressfile, dataI, &hpar);
			if(show_progress == ALL || show_progress == Q)
				writefits_plane(Qprogressfile, dataQ, &hpar);
			if(show_progress == ALL || show_progress == U)
				writefits_plane(Uprogressfile, dataU, &hpar);
			if(show_progress == ALL || show_progress == V)
				writefits_plane(Vprogressfile, dataV, &hpar);
		}
	}


	fprintf(chisqfile, "#Scan weaving\n");
	
	chisqglobalprev = INFINITY;
	//do the scan by scan weaving
	for (ord=0; ord<=scan_order; ord++) 
	//ord = scan_order;
	{
		count = 0;
		printf("Scan weaving order %i\n", ord);
//		fprintf(logfile,"Scan weaving order %i\n", ord);

		do {
			chisqglobal = 0.0;
			for (r=0; r<wappdata->numDays; r++) 
			{ 
				ScanDayData * daydata = &wappdata->scanDayData[r];
				chisqday = 0.0;
				int scancount =0;
				for (i=0; i<daydata->numScans; i++) {
					chisqtmp = 0.0;
					if(daydata->numScans)
					{
						chisqtmp = scan_weave(daydata, i, ord, loop_gain, 1);
						//printf("chisqtmp: %f\n",chisqtmp);
						scancount++;
					}
					chisqday += chisqtmp;
				}
				if(daydata->numScans)
					chisqday /= scancount;
				chisqglobal += chisqday;
				//printf("weaved day:%i chisq:%f\n", r, chisqday);

				/*
				if (show_progress) {
					printf("writing progress plane %i\n", n3++);
					grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
					writefits_plane(progressfile, dataQ, &hpar);
				}
				*/
			}
			chisqglobal /= wappdata->numDays;

			count++;
			globalchange = chisqglobalprev - chisqglobal;
			chisqglobalprev = chisqglobal;

			printf("INFO: scan iteration:%i global chisq:%g change:%g\n", count, chisqglobal, globalchange);
			fprintf(chisqfile, "%f %f\n", chisqglobal, globalchange);

		} while (globalchange > loop_epsilon && count < 10);

		//done all weaves for this order, so write a progress plane
		if (show_progress) 
		{
			printf("Writing progress plane %i\n", n3++);
			grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
			if(show_progress == ALL || show_progress == I)
				writefits_plane(Iprogressfile, dataI, &hpar);
			if(show_progress == ALL || show_progress == Q)
				writefits_plane(Qprogressfile, dataQ, &hpar);
			if(show_progress == ALL || show_progress == U)
				writefits_plane(Uprogressfile, dataU, &hpar);
			if(show_progress == ALL || show_progress == V)
				writefits_plane(Vprogressfile, dataV, &hpar);			
		}

	}


	// finalize the progress cube as needed
	if (show_progress) 
	{
		//free memory
		free(dataI);
		free(dataQ);
		free(dataU);
		free(dataV);
		free(weight);

		//close the fits cube
		if(show_progress == ALL || show_progress == I)
		{
			writefits_pad_end(Iprogressfile, &hpar);
			fclose(Iprogressfile);
			//change the fits header on disk to update the number of planes of the cube
			hpar.naxis[2] = n3;
			Iprogressfile = fopen("Iprogress.fits", "r+");
			writefits_header(Iprogressfile, &hpar);
			fclose(Iprogressfile);
		}
		if(show_progress == ALL || show_progress == Q)
		{
			writefits_pad_end(Qprogressfile, &hpar);
			fclose(Qprogressfile);
			//change the fits header on disk to update the number of planes of the cube
			hpar.naxis[2] = n3;
			Qprogressfile = fopen("Qprogress.fits", "r+");
			writefits_header(Qprogressfile, &hpar);
			fclose(Qprogressfile);
		}
		if(show_progress == ALL || show_progress == U)
		{
			writefits_pad_end(Uprogressfile, &hpar);
			fclose(Uprogressfile);
			//change the fits header on disk to update the number of planes of the cube
			hpar.naxis[2] = n3;
			Uprogressfile = fopen("Uprogress.fits", "r+");
			writefits_header(Uprogressfile, &hpar);
			fclose(Uprogressfile);
		}
		if(show_progress == ALL || show_progress == V)
		{
			writefits_pad_end(Vprogressfile, &hpar);
			fclose(Vprogressfile);
			//change the fits header on disk to update the number of planes of the cube
			hpar.naxis[2] = n3;
			Vprogressfile = fopen("Vprogress.fits", "r+");
			writefits_header(Vprogressfile, &hpar);
			fclose(Vprogressfile);
		}
	}
}



/* main entry point to the balancing operation.
 * After this function returns, the entire frequency should be balanced
 * day to day.
 * decmin and decmax set ranges where data outside these decs (in the endpoints) 
 * will be ignored.
 * loopgain - a multiplier that is used to limit the amount of change on each iteration.
 * loop_epsilon - the threshold where any percent changes less than this signal the 
 * termiation criteria of the iterations.
 */
void balance_data(FluxWappData * wappdata, MapMetaData *md, int day_order, int scan_order, float loop_gain, float loop_epsilon, int show_progress)
{
	FILE * chisqfile;

	chisqfile = fopen ("chisq.dat", "w");
	if (chisqfile == NULL) {
		printf("ERROR: unable to open chisq.dat\n");
		return;
	}

	printf("Finding crossing points\n");
	find_intersections(wappdata);

	printf("Performing basket weaving\n");
	printf("Loop_gain: %g Loop_epsilon: %g\n", loop_gain, loop_epsilon);

	//fprintf(chisqfile, "#day chisqmax, change, chisqglobal, globalchange\n");
	basket_weave(wappdata, chisqfile, day_order, scan_order, loop_gain, loop_epsilon, md, show_progress);
	//mesh_weave(wappdata);

	fclose (chisqfile);
}




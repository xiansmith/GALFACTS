#include "balance.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jsd/jsd_futil.h"
#include "jsd/jsd_util.h"
#include "jsd/jsd_fit.h"
#include "programs/fitsio.h"
#include "grid.h"
#include "stats.h"
#include "string.h"
#include <stdarg.h>

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
	
if (isnan(*pRA))
	printf("FOUND NULL RA!\n");

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
	fprintf(file, "COLOUR BLUE\n"); //the cross lines
	fprintf(file, "LINE W %f %f %f %f\n", x1, y1, x2, y2);
	fprintf(file, "LINE W %f %f %f %f\n", x3, y3, x4, y4);
	fprintf(file, "COLOUR GREEN\n"); //the intersection point
	fprintf(file, "CROSS W %f %f 0.004 0.004\n", *pRA, *pDEC);
	fprintf(file, "COLOUR YELLOW\n"); //two points that are closest to the intersection point
	fprintf(file, "DOT W %f %f\n", curr->records[*currpos].RA, curr->records[*currpos].DEC );
	fprintf(file, "DOT W %f %f\n", ref->records[*refpos].RA, ref->records[*refpos].DEC );

	return 1;
}


/*
static void print_scandata(FILE * file, ScanData * scandata)
{
	int i;
	for (i=0; i<scandata->num_records; i++) {
		fprintf(file, "%f %f\n", scandata->records[i].RA,  scandata->records[i].DEC);
	}
}
*/

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
		sprintf(filename, "cross%s_beam%d.ann", refday->mjd, r%7); //ssg to fix for overwriting of files
		crossfile = fopen(filename, "w");

		for (i=0; i<refScanDay->numScans; i++) 
		{
			ScanData *refscan = &refScanDay->scans[i];
			refscan->num_cross_points = 0;
			fprintf(crossfile, "#scan %i\n", i);
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
if (isnan(RA))
	printf("FOUND NULL RA!\n");
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
		fclose(crossfile);
	}
}



/*
Does a polynomial fit for the few points around a particular point.
records - the data array
size - size of the data array
pos - position within the data array to fit
c[IQUV] - the output fit parameters
order - the order of the fit
*/
static void cross_line(const FluxRecord records[], int size, int pos, double cI[], double cQ[], double cU[], double cV[], int order)
{
	int k, p;
	//TODO: make these constants parameters
	const int stroke = 3; //number of points on either side
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



//TODO: make this a parameter as needed
#define CROSS_ORDER 0


//compute a total chisq value for the entire day to compare fit
static double global_chisq(ScanDayData *daydata)
{
	//iterate over the days
	
	return 0.0;
}


void plot_display(double X[], double Y[], int size, double C[], int order)
{
	FILE *fh;
	int i;

	fh = fopen("/tmp/plot.cmd", "w");
	fprintf(fh, "set output '/tmp/plot.png'\n");
	fprintf(fh, "set term png colour\n");
	jsd_print_poly(fh, C, order);
	fprintf(fh, "plot '/tmp/plot.data' with linespoints, f%i(x)\n", order);
	fclose(fh);

	fh = fopen("/tmp/plot.data", "w");
	for (i=0; i<size; i++) {
		fprintf(fh, "%g %g\n", X[i], Y[i]);
	}
	fclose(fh);
		
	system("gnuplot /tmp/plot.cmd; display /tmp/plot.png");
}

static double scan_segment_weave(ScanData *scan, int start, int end, int order, float loop_gain, int min_points)
{
	double dI[MAX_NUM_DAYS*MAX_NUM_SCANS], dQ[MAX_NUM_DAYS*MAX_NUM_SCANS], dU[MAX_NUM_DAYS*MAX_NUM_SCANS], dV[MAX_NUM_DAYS*MAX_NUM_SCANS], dRA[MAX_NUM_DAYS*MAX_NUM_SCANS];
	//TODO: cX size is related to the order, not days
	double cI[MAX_NUM_DAYS], cQ[MAX_NUM_DAYS], cU[MAX_NUM_DAYS], cV[MAX_NUM_DAYS];
	const float nsigma = 2.0;
	int num_delta = 0;
	int j,k;
	double chisq;
	double min, max;
	int err;

//printf("scan_segment_weave: %i:%i\n", start, end);
	if (end-start < min_points) 
		return 0.0;

	for (j=start; j<end; j++) 
	{
		double refI, refQ, refU, refV;
		double crossI, crossQ, crossU, crossV;
		CrossingPoint *crossPoint =  &scan->crossPoints[j];
		ScanData *crossScan = crossPoint->crossScan;
		double RA = crossPoint->RA;
//if (isnan(RA))
//	printf("RA IS NULL!\n");

		cross_line(scan->records, scan->num_records, crossPoint->ref_pos, cI, cQ, cU, cV, CROSS_ORDER);
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

		dI[num_delta] = refI - crossI;
		dQ[num_delta] = refQ - crossQ;
		dU[num_delta] = refU - crossU;
		dV[num_delta] = refV - crossV;
		dRA[num_delta] = RA;
		num_delta++;
	}

	if (num_delta < 2) return 0.0;

	jsd_minmax(dRA, num_delta, &min, &max);
	jsd_normalize(dRA, num_delta, min, max);

	err = jsd_poly_fit(dRA, dI, num_delta, nsigma, cI, order, &chisq);
	if (err) plot_display(dRA, dI, num_delta, cI, order);
	jsd_poly_fit(dRA, dQ, num_delta, nsigma, cQ, order, &chisq);
	jsd_poly_fit(dRA, dU, num_delta, nsigma, cU, order, &chisq);
	jsd_poly_fit(dRA, dV, num_delta, nsigma, cV, order, &chisq);

	for (k = scan->crossPoints[start].ref_pos; k < scan->crossPoints[end].ref_pos; k++) 
	{
		double RA = NORMALIZE(scan->records[k].RA, min, max);
//if (isnan(RA))
//	printf("RA IS NULL!\n");
		scan->records[k].stokes.I -= jsd_poly_eval(RA, cI, order) * loop_gain;
		scan->records[k].stokes.Q -= jsd_poly_eval(RA, cQ, order) * loop_gain;
		scan->records[k].stokes.U -= jsd_poly_eval(RA, cU, order) * loop_gain;
		scan->records[k].stokes.V -= jsd_poly_eval(RA, cV, order) * loop_gain;
	}

	return compute_clean_mean(dI, num_delta, nsigma);
 
/*
	//apply changes to the data
	{
		//start and end pos in the data have to extend off the ends of the segment
		//when we are on the first or last segment of this scan
		int start_pos = (start<=0) ? 0 : scan->crossPoints[start].ref_pos;
		int end_pos = (end>=scan->num_cross_points-1) ? scan->num_records-1 : scan->crossPoints[end].ref_pos;

		//reject ourliers, and just computes an average delta for each stokes
		double dIavg = compute_clean_mean(dI, num_delta, nsigma);
		double dQavg = compute_clean_mean(dQ, num_delta, nsigma);
		double dUavg = compute_clean_mean(dU, num_delta, nsigma);
		double dVavg = compute_clean_mean(dV, num_delta, nsigma);

		if (!isfinite(dIavg+dQavg+dUavg+dVavg)) {
			printf("WARN: delta is NAN\n");
			return 0.0; 
		}

		//subtract a portion of the average from from the data
		for (k = start_pos; k < end_pos; k++) 
		{
			scan->records[k].stokes.I -= dIavg * loop_gain;
			scan->records[k].stokes.Q -= dQavg * loop_gain;
			scan->records[k].stokes.U -= dUavg * loop_gain;
			scan->records[k].stokes.V -= dVavg * loop_gain;
		}

		//return a chisq value
		return dIavg*dIavg;
	}
*/

}


/*
static double scan_segment_weave(ScanData *scan, int start, int end, int order, float loop_gain, int min_points)
{
	double dI[MAX_NUM_DAYS*MAX_NUM_SCANS], dQ[MAX_NUM_DAYS*MAX_NUM_SCANS], dU[MAX_NUM_DAYS*MAX_NUM_SCANS], dV[MAX_NUM_DAYS*MAX_NUM_SCANS];
	//TODO: cX size is related to the order, not days
	double cI[MAX_NUM_DAYS], cQ[MAX_NUM_DAYS], cU[MAX_NUM_DAYS], cV[MAX_NUM_DAYS];
	double dRA[MAX_NUM_DAYS*MAX_NUM_SCANS];
	const float nsigma = 3.0;
	int num_delta = 0;
	int j,k;
	double rv = 0.0;
	double min, max;
	double chisq;


	if (end-start < min_points) return 0.0;

	for (j=start; j<end; j++) 
	{
		double refI, refQ, refU, refV;
		double crossI, crossQ, crossU, crossV;
		CrossingPoint *crossPoint =  &scan->crossPoints[j];
		ScanData *crossScan = crossPoint->crossScan;

		//get RA and DEC from the crossing point structure
		double RA = crossPoint->RA;
		//double DEC = crossPoint->DEC;

		cross_line(scan->records, scan->num_records, crossPoint->ref_pos, cI, cQ, cU, cV, CROSS_ORDER);
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

		dI[num_delta] = refI - crossI;
		dQ[num_delta] = refQ - crossQ;
		dU[num_delta] = refU - crossU;
		dV[num_delta] = refV - crossV;
		dRA[num_delta] = RA;
		num_delta++;
	}

	//jsd_minmax(dRA, num_delta, &min, &max);
	//jsd_normalize(dRA, num_delta, min, max);

	jsd_poly_fit(dRA, dI, num_delta, nsigma, cI, order, &chisq);
	jsd_poly_fit(dRA, dQ, num_delta, nsigma, cQ, order, &chisq);
	jsd_poly_fit(dRA, dU, num_delta, nsigma, cU, order, &chisq);
	jsd_poly_fit(dRA, dV, num_delta, nsigma, cV, order, &chisq);

	for (k = scan->crossPoints[start].ref_pos; k < scan->crossPoints[end].ref_pos; k++) 
	{
		//double RA = NORMALIZE(scan->records[k].RA, min, max);
		double RA = scan->records[k].RA;
		scan->records[k].stokes.I -= jsd_poly_eval(RA, cI, order) * loop_gain;
		scan->records[k].stokes.Q -= jsd_poly_eval(RA, cQ, order) * loop_gain;
		scan->records[k].stokes.U -= jsd_poly_eval(RA, cU, order) * loop_gain;
		scan->records[k].stokes.V -= jsd_poly_eval(RA, cV, order) * loop_gain;
	}

	//split this segment into two segments and recurse
	rv += scan_segment_weave(scan, start, (end+start)/2, order, loop_gain, min_points);
	rv += scan_segment_weave(scan, (end+start)/2, end, order, loop_gain, min_points);

	return rv/2.0;
}
*/

/*
basket weaving a single days data
Iterates over each scan
Breaks down each scan into the specified number of non-overlappng segments
Performs scan segment weaving on each segment
Returns the average delta
*/
static double day_weave(ScanDayData *daydata, int segments, float loop_gain)
{
	int i,j;
	double delta = 0.0;
	int min_points = 3;
	int order=0;
printf("day_weave segments: %i\n", segments);
	for (i=0; i<daydata->numScans; i++) 
	{
		int num = daydata->scans[i].num_cross_points;
		for (j=0; j<segments; j++) 
		{
			int start = (num*j)/segments;
			int end = (num*(j+1))/segments;
			delta += scan_segment_weave(&daydata->scans[i], start, end, order, loop_gain, min_points);
		}
	}
	return delta/daydata->numScans;
}



#define PERCENT_CHANGE(A,B) ((A-B)/B)

static void basket_weave(FluxWappData *wappdata, FILE * chisqfile, int day_order, int scan_order, float loop_gain, float loop_epsilon, MapMetaData * md, int show_progress)
{
	int r;
	int count;
	double chisqtmp, chisqglobal, chisqglobalprev;
	double globalchange;
	int ord;
	static header_param_list hpar;
	FILE *progressfile;
	float *dataI, *dataQ, *dataU, *dataV, *weight;
	int n1, n2, n3;
	if (show_progress) 
	{
		n1 = md->n1;
		n2 = md->n2;
		n3 = 0; //don't know how many planes in the cube, it will be set later below

		//start a fits cube
		dataI = (float *) malloc (n1 * n2 * sizeof (float));
		dataQ = (float *) malloc (n1 * n2 * sizeof (float));
		dataU = (float *) malloc (n1 * n2 * sizeof (float));
		dataV = (float *) malloc (n1 * n2 * sizeof (float));
		weight = (float *) malloc (n1 * n2 * sizeof (float));

		init_header_param_list (&hpar);
		hpar.bitpix = -32;
		//hpar.num_axes = 2;
		hpar.num_axes = 3;
		hpar.naxis[0] = n1;
		hpar.naxis[1] = n2;
		hpar.naxis[2] = 0;
		sprintf (hpar.ctype[0], "RA---CAR");
		sprintf (hpar.ctype[1], "DEC---CAR");
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
		sprintf (hpar.object, "Q Basketweaving Progress");

		progressfile = fopen("Qprogress.fits", "w");
		writefits_header(progressfile, &hpar);

		grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
		writefits_plane(progressfile, dataQ, &hpar);
		printf("plane %i\n", n3++);
	}

	//do day by day weaving
	for (ord=0; ord<=day_order; ord++) 
	{
int segments = 0x1 << ord;
chisqglobalprev = INFINITY;
//ord is now segments
		count = 0;
		printf("Day weaving segments %i\n", segments);
		do {
			chisqglobal = 0.0;
			for (r=0; r<wappdata->numDays; r++) 
			{ 
				chisqtmp = day_weave(&wappdata->scanDayData[r], segments, loop_gain);
				chisqglobal += chisqtmp;
				//printf("weaved day:%i chisq:%f\n", r, chisqtmp);
			}
			chisqglobal /= wappdata->numDays;
			
			//done a weave to a set of days, write a progress plane
			//if (show_progress) {
			//	printf("writing progress plane %i\n", n3++);
			//	grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
			//	writefits_plane(progressfile, dataQ, &hpar);
			//}

			count++;
			globalchange = chisqglobalprev - chisqglobal;
			chisqglobalprev = chisqglobal;
			printf("day iteration:%i global chisq:%g change:%g\n", count, chisqglobal, globalchange);
			fprintf(chisqfile, "%f %f\n", chisqglobal, globalchange);

		//} while (globalchange > loop_epsilon && count < 30);
		} while (globalchange > loop_epsilon && segments < 100);

		//done all weaves for this order, so write a progress plane
		if (show_progress) {
			printf("writing progress plane %i\n", n3++);
			grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
			writefits_plane(progressfile, dataQ, &hpar);
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
		writefits_pad_end(progressfile, &hpar);
		fclose(progressfile);

		//change the fits header on disk to update the number of planes of the cube
		hpar.naxis[2] = n3;
		progressfile = fopen("Qprogress.fits", "r+");
		writefits_header(progressfile, &hpar);
		fclose(progressfile);
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

	printf("finding crossing points ...\n");
	find_intersections(wappdata);

	printf("performing basket weaving ...\n");
	printf("loop_gain: %g loop_epsilon: %g\n", loop_gain, loop_epsilon);

	fprintf(chisqfile, "#day chisqmax, change, chisqglobal, globalchange\n");
	basket_weave(wappdata, chisqfile, day_order, scan_order, loop_gain, loop_epsilon, md, show_progress);
	//mesh_weave(wappdata);

	fclose (chisqfile);

}




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
#include "string.h"

// this function computes the intersection of the sent lines
// and returns the intersection point, note that the function assumes
// the lines intersect. the function can handle vertical as well
// as horizontal lines. note the function isn't very clever, it simply
//applies the math
static void line_intersect(double x0,double y0,double x1,double y1,
		double x2,double y2,double x3,double y3,
		double *xi,double *yi)
{

	double a1,b1,c1, // constants of linear equations
	a2,b2,c2,
	det_inv,  // the inverse of the determinant of the coefficient matrix
	m1,m2;    // the slopes of each line

	// compute slopes, note the cludge for infinity, however, this will
	// be close enough

	if ((x1-x0)!=0)
		m1 = (y1-y0)/(x1-x0);
	else
		m1 = INFINITY; 

	if ((x3-x2)!=0)
		m2 = (y3-y2)/(x3-x2);
	else
		m2 = INFINITY;

	a1 = m1;
	a2 = m2;

	b1 = -1;
	b2 = -1;

	c1 = (y0-m1*x0);
	c2 = (y2-m2*x2);

	// compute the inverse of the determinate

	det_inv = 1/(a1*b2 - a2*b1);

	// use Kramers rule to compute xi and yi

	*xi=((b1*c2 - b2*c1)*det_inv);
	*yi=((a2*c1 - a1*c2)*det_inv);

} // end Intersect_Lines


static float swapf_tmp;
#define SWAPF(X,Y) swapf_tmp=*X; *X=*Y; *Y=swapf_tmp;

static int scan_intersect(FILE * file, ScanData *ref, ScanData *curr, float *RA, float *DEC, int *refpos, int *currpos)
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
	x2 = ref->records[i].RA; //TODO: should this be +1?
	y2 = ref->records[i].DEC; 
	x3 = curr->records[j-1].RA; 
	y3 = curr->records[j-1].DEC; 
	x4 = curr->records[j].RA; 
	y4 = curr->records[j].DEC; 

	


/*
	{
		float denom, ua, ub;

		denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1);
		if (denom == 0) return 0;

		ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom; 
		ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom; 

		*RA = x1 + ua * (x2 - x1);
		*DEC = y1 + ua * (y2 - y1);
	}
*/
	{
		double ar, ac, br, bc;

		br = (y1-y2)/(x1-x2);
		ar = y1 - br * x1;
		bc = (y3-y4)/(x3-x4);
		ac = y3 - bc * x3;

		*RA = -(ar-ac)/(br-bc);
		*DEC = ar + br * (*RA);
	}
	
	*refpos = i;
	*currpos = j;

	//out of range check
	if (*RA > refMax || *RA < refMin) {
		return 0;
	}

	//printf("found cross: %f %f\n", *RA, *DEC);
	fprintf(file, "COLOUR BLUE\n");
	fprintf(file, "LINE W %f %f %f %f\n", x1, y1, x2, y2);
	fprintf(file, "LINE W %f %f %f %f\n", x3, y3, x4, y4);
	fprintf(file, "COLOUR GREEN\n");
	fprintf(file, "CIRCLE W %f %f %f\n", *RA, *DEC, 0.001);


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
	int i, j;
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
		sprintf(filename, "cross%s_beam%d.ann", refday->mjd,numDays%7); //ssg to fix for overwriting of files
		crossfile = fopen(filename, "w");
		fprintf(crossfile, "COLOUR BLUE\n");

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
					float RA, DEC;
					int refpos, crosspos;
					ScanData *currscan = &currScanDay->scans[j];
					if (currscan->num_records == 0) continue;
					//SSG
					if (currscan->num_records < 0)
					{
						printf("OOPS...\n");
						continue;
					}
					//SSG
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

						fprintf(crossfile, "COLOUR RED\n");
						fprintf(crossfile, "CROSS W %f %f 0.004 0.004\n", RA, DEC);

						crossPoint = &refscan->crossPoints[refscan->num_cross_points];
						crossPoint->RA = RA;
						crossPoint->crossScan = currscan;
						crossPoint->ref_pos = refpos;
						crossPoint->cross_pos = crosspos;
						refscan->num_cross_points++;
					}
				}
			}
		}
		fclose(crossfile);
	}
}


static void cross_line(FluxRecord records[], int size, int pos, double cI[], double cQ[], double cU[], double cV[], int order)
{
	int k, p;
	int stroke = 3;
	int count = 0;
	double chisq;
	const float nsigma = 4.0;
	double yI[MAX_NUM_DAYS];
	double yQ[MAX_NUM_DAYS];
	double yU[MAX_NUM_DAYS];
	double yV[MAX_NUM_DAYS];
	double RA[MAX_NUM_DAYS];

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

#define CROSS_ORDER 2
static double scan_weave(ScanDayData *daydata, int scan, int order, float loop_gain, int apply)
{
	const float nsigma = 3.0;
	const int stroke = 1;

	int x, j;
	int k;
	double RA, DEC;
	double min, max;
	int num_delta;
	double delta_sum;

	double dI[MAX_NUM_DAYS*MAX_NUM_SCANS], dQ[MAX_NUM_DAYS*MAX_NUM_SCANS], dU[MAX_NUM_DAYS*MAX_NUM_SCANS], dV[MAX_NUM_DAYS*MAX_NUM_SCANS];
	//TODO: cX size is related to the order, 
	double dRA[MAX_NUM_DAYS*MAX_NUM_SCANS];
	double chisq[4];
	double cI[CROSS_ORDER+1], cQ[CROSS_ORDER+1], cU[CROSS_ORDER+1], cV[CROSS_ORDER+1];
	ScanData * refscan;


	num_delta = 0;

	for (x=-stroke; x<=stroke; x++) 
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
			FluxRecord *refRec, *crossRec;

			CrossingPoint *crossPoint =  &refscan->crossPoints[j];
			ScanData *crossScan = crossPoint->crossScan;

			refRec = &refscan->records[crossPoint->ref_pos];
			crossRec = &crossPoint->crossScan->records[crossPoint->cross_pos];

			line_intersect(refRec[-1].RA, refRec[-1].DEC, refRec[1].RA, refRec[1].DEC,
					crossRec[-1].RA, crossRec[-1].DEC, crossRec[1].RA, crossRec[1].DEC,
					&RA, &DEC);

			cross_line(refscan->records, refscan->num_records, crossPoint->ref_pos, cI, cQ, cU, cV, CROSS_ORDER);
			refI = jsd_poly_eval(RA, cI, CROSS_ORDER);
			refQ = jsd_poly_eval(RA, cQ, CROSS_ORDER);
			refU = jsd_poly_eval(RA, cU, CROSS_ORDER);
			refV = jsd_poly_eval(RA, cV, CROSS_ORDER);

			cross_line(crossScan->records, crossScan->num_records, crossPoint->cross_pos, cI, cQ, cU, cV, CROSS_ORDER);
			crossI = jsd_poly_eval(RA, cI, CROSS_ORDER);
			crossQ = jsd_poly_eval(RA, cQ, CROSS_ORDER);
			crossU = jsd_poly_eval(RA, cU, CROSS_ORDER);
			crossV = jsd_poly_eval(RA, cV, CROSS_ORDER);

			dI[num_delta] = refI - crossI;
			dQ[num_delta] = refQ - crossQ;
			dU[num_delta] = refU - crossU;
			dV[num_delta] = refV - crossV;
			dRA[num_delta] = RA;
			num_delta++;
		}
	}

	if (num_delta > order) 
	{
		//apply the curve fits
		if (apply) 
		{
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

		//sum up delta magnitudes and return result
		delta_sum = 0.0;
		for (k=0; k<num_delta; k++) {
			delta_sum += fabs(dI[k]);
		}
		return delta_sum/num_delta;
	}
	else
	{
		return 0.0;
	}
}


static double day_weave(ScanDayData *daydata, int order, float loop_gain, int apply)
{
	int i, j;
	int k;
	double RA, DEC;
	double min, max;

	double dI[MAX_NUM_DAYS*MAX_NUM_SCANS], dQ[MAX_NUM_DAYS*MAX_NUM_SCANS], dU[MAX_NUM_DAYS*MAX_NUM_SCANS], dV[MAX_NUM_DAYS*MAX_NUM_SCANS];
	//TODO: cX size is related to the order, not days
	double cI[MAX_NUM_DAYS], cQ[MAX_NUM_DAYS], cU[MAX_NUM_DAYS], cV[MAX_NUM_DAYS];
	double dRA[MAX_NUM_DAYS*MAX_NUM_SCANS];
	double chisq[4];
	const float nsigma = 3.0;

	int num_delta = 0;

	for (i=0; i<daydata->numScans; i++) 
	{
		ScanData *refscan = &daydata->scans[i];
		int num_cross_points = refscan->num_cross_points;

		for (j=0; j<num_cross_points; j++) 
		{

			double refI, refQ, refU, refV;
			double crossI, crossQ, crossU, crossV;
			FluxRecord *refRec, *crossRec;
			CrossingPoint *crossPoint =  &refscan->crossPoints[j];
			ScanData *crossScan = crossPoint->crossScan;

			refRec = &refscan->records[crossPoint->ref_pos];
			crossRec = &crossPoint->crossScan->records[crossPoint->cross_pos];

			line_intersect(refRec[-1].RA, refRec[-1].DEC, refRec[1].RA, refRec[1].DEC,
					       crossRec[-1].RA, crossRec[-1].DEC, crossRec[1].RA, crossRec[1].DEC,
							&RA, &DEC);

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

			dI[num_delta] = refI - crossI;
			dQ[num_delta] = refQ - crossQ;
			dU[num_delta] = refU - crossU;
			dV[num_delta] = refV - crossV;
			dRA[num_delta] = RA;
			num_delta++;

			//printf("dI=%f dQ=%f dU=%f dV=%f RA=%f\n", dI[j], dQ[j], dU[j], dV[j], RA[j]);
		}
	}
	if (num_delta <= 0) {
		//printf("WARN: no crossing points for day %s\n", daydata->mjd);
		return 0.0;
	}

	jsd_minmax(dRA, num_delta, &min, &max);
	jsd_normalize(dRA, num_delta, min, max);

	//curve fit the dX values
	jsd_poly_fit(dRA, dI, num_delta, nsigma, cI, order, &chisq[0]);
	jsd_poly_fit(dRA, dQ, num_delta, nsigma, cQ, order, &chisq[1]);
	jsd_poly_fit(dRA, dU, num_delta, nsigma, cU, order, &chisq[2]);
	jsd_poly_fit(dRA, dV, num_delta, nsigma, cV, order, &chisq[3]);



	if (apply) {
		//jsd_print_poly(stdout, cI, order);
		//apply the dX values
		for (i=0; i<daydata->numScans; i++) 
		{
			ScanData *refscan = &daydata->scans[i];
			if (!finite(jsd_poly_eval(RA, cI, order))) break; 
			for (k=0; k<refscan->num_records; k++) 
			{
				RA = NORMALIZE(refscan->records[k].RA, min, max);
				refscan->records[k].stokes.I -= jsd_poly_eval(RA, cI, order) * loop_gain;
				refscan->records[k].stokes.Q -= jsd_poly_eval(RA, cQ, order) * loop_gain;
				refscan->records[k].stokes.U -= jsd_poly_eval(RA, cU, order) * loop_gain;
				refscan->records[k].stokes.V -= jsd_poly_eval(RA, cV, order) * loop_gain;
			}
		}
	}
	

	/* Difficult to determine what a good criteria for a good or bad
	day is.  chisq values are not valid, since bad data could have a good fit.
	using the average difference instead.
	*/
	{
	double sum = 0.0;
	for (i=0; i<num_delta; i++) {
		sum += fabs(dI[i]);
	}
	return sum/num_delta;
	}
}


#define PERCENT_CHANGE(A,B) ((A-B)/B)

static void basket_weave(FluxWappData *wappdata, int order, float loop_gain, float loop_epsilon, MapMetaData * md, int show_progress)
{
	int r, i;
	int count;
	double chisqprev, chisqtmp, chisqmax;
	double change;
	int ord;
	int worst_day;

	static header_param_list hpar;
	FILE *progressfile;
	float *dataI, *dataQ, *dataU, *dataV, *weight;
	int n1, n2, n3;
	int scan_order = 0;

	if (show_progress) 
	{
		n1 = md->n1;
		n2 = md->n2;
		n3 = 0;

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
		sprintf (hpar.ctype[1], "Iteration");
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
		sprintf (hpar.object, "I Basketweaving Progress");

		progressfile = fopen("Iprogress.fits", "w");
		writefits_header(progressfile, &hpar);

		grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
		writefits_plane(progressfile, dataI, &hpar);
		printf("plane %i\n", n3++);
	}

	//do day by day weaving
	for (ord=0; ord<=order; ord++) 
	{
		chisqprev = INFINITY;
		count = 0;
		printf("Basket weaving order %i\n", ord);

		do {
			chisqmax = 0.0;
			for (r=0; r<wappdata->numDays; r++) { 
				chisqtmp = day_weave(&wappdata->scanDayData[r], ord, loop_gain, 0);
//printf("r: %i chisq: %f\n", r, chisqtmp);
				if (chisqtmp > chisqmax) {
					chisqmax = chisqtmp;
					worst_day = r;
				}
			}
			printf("weaving day:%i chisq:%f\n", worst_day, chisqmax);
			day_weave(&wappdata->scanDayData[worst_day], ord, loop_gain, 1);
			
			if (show_progress) {
				printf("writing progress plane %i\n", n3++);
				grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
				writefits_plane(progressfile, dataI, &hpar);
			}

			count++;
			change = chisqprev - chisqmax;
			chisqprev = chisqmax;
			printf("iteration: %i change: %g\n", count, change);

		} while (/*change > loop_epsilon && */count < 100);

		if (show_progress) {
			printf("writing progress plane %i\n", n3++);
			grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
			writefits_plane(progressfile, dataI, &hpar);
		}
	}


	//do the scan by scan weaving
	for (ord=0; ord<=scan_order; ord++) 
	{
		count = 0;
		printf("Scan weaving order %i\n", ord);

		do {

			chisqmax = 0.0;
			for (r=0; r<wappdata->numDays; r++) 
			{ 
				ScanDayData * daydata = &wappdata->scanDayData[r];

				for (i=0; i<daydata->numScans; i++) {
					chisqtmp = scan_weave(daydata, i, ord, loop_gain, 1);
					if (chisqtmp > chisqmax) {
						chisqmax = chisqtmp;
					}
				}
				if (show_progress) {
					printf("writing progress plane %i\n", n3++);
					grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
					writefits_plane(progressfile, dataI, &hpar);
				}
			}
			printf("chisqmax: %f\n", chisqmax);
			count++;

		} while (count < 5);
	}


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
		progressfile = fopen("Iprogress.fits", "r+");
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
void balance_data(FluxWappData * wappdata, MapMetaData *md, int order, float loop_gain, float loop_epsilon, int show_progress)
{
	printf("finding crossing points ...\n");
	find_intersections(wappdata);

	printf("performing basket weaving ...\n");
	printf("loop_gain: %f loop_epsilon: %f\n", loop_gain, loop_epsilon);
	basket_weave(wappdata, order, loop_gain, loop_epsilon, md, show_progress);

}




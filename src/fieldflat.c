#include "fieldflat.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "jsd/jsd_fit.h"

#define MAX_SCAN_SIZE 3000

/*
 * The intention here is to remove the large spacial frequencies from the data.
 * The data along each scan is fit with a low order polynomial, which is subtracted out.
 */
void flatten_field(FluxWappData * wappdata, float nsigma, int order)
{
	int r, i, j;
	int num_records;
	double RA;
	double min, max;
	double dI[MAX_SCAN_SIZE];
	double dRA[MAX_SCAN_SIZE];
	double chisq[4];
	double cI[20];

	for (r=0; r<wappdata->numDays; r++) 
	{ 
		ScanDayData * daydata = &wappdata->scanDayData[r];

		for (i=0; i<daydata->numScans; i++) 
		{
			ScanData * scan = &daydata->scans[i];
			num_records = scan->num_records;

			//printf("day: %i  scan: %i\n", r, i);
			
			//check for array overflow
			if (num_records > MAX_SCAN_SIZE) {
				printf("ERROR: MAX_SCAN_SIZE (%i) is smaller than num_records (%i)\n", MAX_SCAN_SIZE, num_records);
				return;
			}
			
			//create the data arrays for fitting
			for (j=0; j<num_records; j++) {
				dRA[j] = scan->records[j].RA;
				dI[j] = scan->records[j].stokes.I;
			}

			//fit the data on the scan with a polynomial
			jsd_minmax(dRA, num_records, &min, &max);
			jsd_normalize(dRA, num_records, min, max);
			jsd_poly_fit(dRA, dI, num_records, nsigma, cI, order, &chisq[0]);

			//subtract out the polynomial evaluation from the data
			for (j=0; j<num_records; j++) {
				RA = NORMALIZE(scan->records[j].RA, min, max);
				scan->records[j].stokes.I -= jsd_poly_eval(RA, cI, order);
			}
		}
	}
}



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fluxdata.h"
#include "map.h"
#include "grid.h"
#include "programs/fitsLib.h"

#define FWHM      2.0                    // FWHM of psf in arc minutes
//#define FWHM      1.5                    // FWHM of psf in arc minutes

//original PSF
//(exp(-2.772589*(offset/fwhm)*(offset/fwhm)));


static float* psf_map; //patch X patch in dimension
static int psf_map_size;
void init_psf_map(float fwhm, float cellsize, int beamwidths)
{
	int i,j;
	float scale = cellsize * 60.0; //scaled to arcminutes
	int radius = beamwidths*fwhm / scale;
	float maxdist = radius * scale;

	psf_map_size = radius*2+1;
	psf_map = malloc (sizeof(double) * psf_map_size*psf_map_size);
	for (i=0; i<psf_map_size; i++) {
		for (j=0; j<psf_map_size; j++) {
			int xpos = i-radius;
			int ypos = j-radius;
			double dist = sqrt(xpos*xpos + ypos*ypos) * scale;
			psf_map[j*psf_map_size+i] = (dist > maxdist) ? 0.0 : exp(-2.772589*(dist/fwhm)*(dist/fwhm));
		}
	}

	{
	header_param_list hpar;
	init_header_param_list (&hpar);  /* initialize parameter records */
	hpar.bitpix = -32;
	hpar.num_axes = 2;
	hpar.naxis[0] = psf_map_size;
	hpar.naxis[1] = psf_map_size;
	sprintf (hpar.ctype[0], "RA---CAR");
	sprintf (hpar.ctype[1], "DEC---CAR");
	hpar.crval[0] = 0;
	hpar.crval[1] = 0;
	hpar.crpix[0] = radius+1;
	hpar.crpix[1] = radius+1;
	hpar.cdelt[0] = -cellsize;
	hpar.cdelt[1] = cellsize;
	hpar.equinox = 2000.0;
	sprintf (hpar.bunit, "Weight");
	sprintf (hpar.object, "PSF");
	writefits_map ("psf_map.fits", psf_map, &hpar);
	}
}

static double * psf_lookup_table;
static double psf_step;
static int psf_size;
//static double dsqrarg;

//#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
//#define PSF(offset) ( exp(-2.772589 * DSQR(offset/fwhm)) )
//#define psf(offset) (psf_lookup_table[(int)(offset/psf_step)])

void init_psf_lookup_table(int size, float maxVal)
{
	int i;
	register float offset = 0.0;
	psf_step = maxVal / (float)size;
	psf_size = size;
	psf_lookup_table = malloc (sizeof(double) * size);
	for (i=0; i<size; i++) {
		offset += psf_step;
		psf_lookup_table[i] = (exp(-2.772589*(offset/FWHM)*(offset/FWHM)));
	}
}

//SSGstatic double psf(float offset) {
double psf(float offset) {
	return psf_lookup_table[(int)(offset/psf_step)];
}

static float dist(float dx, float dy, float scale) {
	return sqrt(dx*dx+dy*dy) * scale;
}

static void grid_fluxdaydata(const FluxDayData *daydata, const MapMetaData *md, float * dataI, float * dataQ, float * dataU, float * dataV, float * weight)
{
	int r, i, j;
	int n1, n2;

	//dimensions of the map
	n1 = md->n1;
	n2 = md->n2;

	for (r=0; r<daydata->numRecords; r++) 
	{
		int x, y;
		const FluxRecord * rec = &daydata->records[r];
		if (!finite(rec->stokes.I)) continue;

		//dont grid on data that is outside the map
		if ((rec->DEC > md->decmax) || (rec->DEC < md->decmin) ||
			(rec->RA > md->ramax) || (rec->RA < md->ramin))
			continue; 

		//find the pixel coords for this measurement
		x = n1-(rec->RA-md->ramin)/md->cellsize;
		y = (rec->DEC-md->decmin)/md->cellsize;

		for (i=0; i<psf_map_size; i++) 
		{
			for (j=0; j<psf_map_size; j++) 
			{
				//pixel coords for this particular weight
				int m = x+(i-psf_map_size/2);
				int n = y+(j-psf_map_size/2);

				if ( (m >= 0) && (m < n1) ) {
					if( (n >= 0) && (n < n2) ) {
						double spread = psf_map[j*psf_map_size+i];
						int indx = n*n1+m;
						dataI[indx] +=  rec->stokes.I * spread;
						dataQ[indx] +=  rec->stokes.Q * spread;          
						dataU[indx] +=  rec->stokes.U * spread;          
						dataV[indx] +=  rec->stokes.V * spread;
						weight[indx] += spread;
					}
				}
			} 
		} 
	}
}


/*
static void grid_fluxdaydata(const FluxDayData *daydata, const MapMetaData *md, float * dataI, float * dataQ, float * dataU, float * dataV, float * weight)
{
	int r, m, n;
	float x, y, scale;
	int n1, n2;
	int patch;

	n1 = md->n1;
	n2 = md->n2;
	patch = md->patch;
	scale = md->cellsize*60.0;

	for (r=0; r<daydata->numRecords; r++) 
	{
		const FluxRecord * rec = &daydata->records[r];
		if (!finite(rec->stokes.I)) continue;

		//dont grid on data that is outside the map
		if ((rec->DEC > md->decmax) || (rec->DEC < md->decmin) ||
			(rec->RA > md->ramax) || (rec->RA < md->ramin))
			continue; 

		x = n1-(rec->RA-md->ramin)/md->cellsize;
		y = (rec->DEC-md->decmin)/md->cellsize;

		for (m=-patch; m<=patch; m++) 
		{
			for (n=-patch; n<=patch; n++) 
			{
				int i, j;
				i = (int)(round(x-(float)m));
				j = (int)(round(y-(float)n));
				//distance = sqrt(dx*dx+dy*dy)*md->cellsize*60.0;
				if ( (i >= 0) && (i < n1) ) {
					if( (j >= 0) && (j < n2) ) {
						float dx = x - (float)i;
						float dy = y - (float)j;
						double distance = dist(dx, dy, scale);
						double spread = psf(distance);
						int indx = i+j*n1;
						dataI[indx] +=  rec->stokes.I * spread;
						dataQ[indx] +=  rec->stokes.Q * spread;          
						dataU[indx] +=  rec->stokes.U * spread;          
						dataV[indx] +=  rec->stokes.V * spread;
						weight[indx] += spread;
					}
				}
			} 
		} 
	}
}
*/

void grid_data(const FluxWappData *wappdata, const MapMetaData *md, float * dataI, float * dataQ, float * dataU, float * dataV, float *weight)
{
	int i, j, d;
	int n1, n2, n3;
	int numbytes;

	n1 = md->n1;
	n2 = md->n2;
	n3 = md->n3;

	numbytes = n1 * n2 * sizeof(float);
	memset(weight, 0, numbytes);
	memset(dataI, 0, numbytes);
	memset(dataQ, 0, numbytes);
	memset(dataU, 0, numbytes);
	memset(dataV, 0, numbytes);
	//TODO: currently ignoring the upscan/downscan functionality
	for (d=0; d<wappdata->numDays; d++) 
	{
		if(d%7 == 0)
			printf("Day %s: ",wappdata->daydata[d].mjd);
		//if (md->gridtype & GRID_UP_SCANS) {
		printf("beam %d ",d%7);
		fflush(stdout);
		if(d%7 == 6)
			printf("\n");
		grid_fluxdaydata(&wappdata->daydata[d], md, dataI, dataQ, dataU, dataV, weight);
		//}
	}
	printf("\n");

	//------------------- divide sums by weights ----------------------
	for(j=0;j<n2;j++) {
		for(i=0;i<n1;i++) {
			double w = weight[i+j*n1];
			if (w > 0.0) {
				dataI[i+j*n1] /= w;
				dataQ[i+j*n1] /= w;
				dataU[i+j*n1] /= w;
				dataV[i+j*n1] /= w;
				//printf("%5i %5i %8.6f %8.6f %8.6f %8.6f\n", i, j, 
				//	dataI[i+j*n1], dataQ[i+j*n1], dataU[i+j*n1], dataV[i+j*n1]);
			} else {
				dataI[i+j*n1] = NAN;
				dataQ[i+j*n1] = NAN;
				dataU[i+j*n1] = NAN;
				dataV[i+j*n1] = NAN;
			}
		}
	}

}




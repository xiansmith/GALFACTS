#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsio.h"
#include "map.h"

void write_fits_maps(const char * wapp, MapMetaData *md, float dataI[], float dataQ[], float dataU[], float dataV[])
{
        header_param_list hpar;
	char filename[31+1];

        init_header_param_list (&hpar);  /* initialize parameter records */
        hpar.bitpix = -32;
        hpar.num_axes = 2;
        hpar.naxis[0] = md->n1;
        hpar.naxis[1] = md->n2;
        sprintf (hpar.ctype[0], "RA---CAR");
        sprintf (hpar.ctype[1], "DEC---CAR");
        hpar.crval[0] = md->RAcen;          /* hours */
        hpar.crval[1] = md->DECcen;               /* degrees */
        hpar.crpix[0] = 0.5 + md->n1 / 2.0; /* image center in pixels */
        hpar.crpix[1] = 0.5 + md->n2 / 2.0;
        hpar.cdelt[0] = -md->cellsize;                     /* degrees */
        hpar.cdelt[1] = md->cellsize;                     /* degrees */
	hpar.equinox = 2000.0;
	//hpar.epoch = 2000.0;
        sprintf (hpar.bunit, "Kelvin");
        sprintf (hpar.telescope, "Arecibo");
		
//TODO: datamin, datamax, obsfreq, exptime, epoch,

        sprintf (hpar.object, "GALFACTS Test Region %s Stokes I", wapp);
		snprintf (filename, 31, "%s_Iavg.fits", wapp);
        writefits_map (filename, dataI, &hpar);

        sprintf (hpar.object, "GALFACTS Test Region %s Stokes Q", wapp);
		snprintf (filename, 31, "%s_Qavg.fits", wapp);
        writefits_map (filename, dataQ, &hpar);

        sprintf (hpar.object, "GALFACTS Test Region %s Stokes U", wapp);
		snprintf (filename, 31, "%s_Uavg.fits", wapp);
        writefits_map (filename, dataU, &hpar);

        sprintf (hpar.object, "GALFACTS Test Region %s Stokes V", wapp);
		snprintf (filename, 31, "%s_Vavg.fits", wapp);
        writefits_map (filename, dataV, &hpar);

}

static header_param_list hparI;
static header_param_list hparQ;
static header_param_list hparU;
static header_param_list hparV;
static header_param_list hparW;
static FILE * fitsI;
static FILE * fitsQ;
static FILE * fitsU;
static FILE * fitsV;
static FILE * fitsW;

static void create_header_param_list(header_param_list * hpar_ptr, MapMetaData *md)
{
	init_header_param_list (hpar_ptr);  /* initialize parameter records */
	hpar_ptr->bitpix = -32;
	hpar_ptr->num_axes = 3;
	hpar_ptr->naxis[0] = md->n1;
	hpar_ptr->naxis[1] = md->n2;
	hpar_ptr->naxis[2] = md->n3;
	sprintf (hpar_ptr->ctype[0], "RA---CAR");
	sprintf (hpar_ptr->ctype[1], "DEC---CAR");
	sprintf (hpar_ptr->ctype[2], "Frequency");
	hpar_ptr->crval[0] = md->RAcen;		/* degrees */
	hpar_ptr->crval[1] = md->DECcen;               /* degrees */
	hpar_ptr->crval[2] = md->fstart;               /* khz */
	hpar_ptr->crpix[0] = 0.5 + md->n1 / 2.0;	/* image center in pixels */
	hpar_ptr->crpix[1] = 0.5 + md->n2 / 2.0;
	hpar_ptr->crpix[2] = 0.0;
	hpar_ptr->cdelt[0] = -md->cellsize;                     /* degrees */
	hpar_ptr->cdelt[1] = md->cellsize;                     /* degrees */
	hpar_ptr->cdelt[2] = 100.0/256.0;	//TODO: magic numbers				/* Mhz */
	hpar_ptr->equinox = 2000.0;		//year of coordinate system
	//hpar_ptr->epoch = 2004.???;		//year observations were taken

	sprintf (hpar_ptr->bunit, "Kelvin");
	sprintf (hpar_ptr->telescope, "Arecibo");
}

void start_fits_cubes(const char * wapp, MapMetaData *md)
{
	char filename[31+1];

	create_header_param_list(&hparI, md);
	create_header_param_list(&hparQ, md);
	create_header_param_list(&hparU, md);
	create_header_param_list(&hparV, md);
	create_header_param_list(&hparW, md);

	sprintf (hparW.bunit, "Weight");

	sprintf (hparI.object, "%s Stokes I", md->title);
	sprintf (hparQ.object, "%s Stokes Q", md->title);
	sprintf (hparU.object, "%s Stokes U", md->title);
	sprintf (hparV.object, "%s Stokes V", md->title);
	sprintf (hparW.object, "%s Weight Cube", md->title);

	snprintf (filename, 31, "%s_Icube.fits", wapp);
	fitsI = fopen(filename, "w");
	snprintf (filename, 31, "%s_Qcube.fits", wapp);
	fitsQ = fopen(filename, "w");
	snprintf (filename, 31, "%s_Ucube.fits", wapp);
	fitsU = fopen(filename, "w");
	snprintf (filename, 31, "%s_Vcube.fits", wapp);
	fitsV = fopen(filename, "w");
	snprintf (filename, 31, "%s_Weightcube.fits", wapp);
	fitsW = fopen(filename, "w");

	if (fitsI == NULL || fitsQ == NULL || fitsU == NULL || fitsV == NULL || fitsW==NULL) {
		printf("ERROR: failed to open fits files\n");
	}

	writefits_header(fitsI, &hparI);
	writefits_header(fitsQ, &hparQ);
	writefits_header(fitsU, &hparU);
	writefits_header(fitsV, &hparV);
	writefits_header(fitsW, &hparW);
}

void write_fits_planes(float dataI[], float dataQ[], float dataU[], float dataV[], float dataW[])
{
	writefits_plane(fitsI, dataI, &hparI);
	writefits_plane(fitsQ, dataQ, &hparQ);
	writefits_plane(fitsU, dataU, &hparU);
	writefits_plane(fitsV, dataV, &hparV);
	writefits_plane(fitsW, dataW, &hparW);
}

void finish_fits_cubes()
{
	writefits_pad_end(fitsI, &hparI);
	writefits_pad_end(fitsQ, &hparQ);
	writefits_pad_end(fitsU, &hparU);
	writefits_pad_end(fitsV, &hparV);
	writefits_pad_end(fitsW, &hparW);

	fclose(fitsI);
	fclose(fitsQ);
	fclose(fitsU);
	fclose(fitsV);
	fclose(fitsW);

	fitsI = NULL;
	fitsQ = NULL;
	fitsU = NULL;
	fitsV = NULL;
	fitsW = NULL;
}



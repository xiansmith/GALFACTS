#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsLib.h"
#include "map.h"
#include "common.h"
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
        sprintf (hpar.ctype[1], "DEC--CAR");
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
	sprintf (hpar_ptr->ctype[1], "DEC--CAR");
	sprintf (hpar_ptr->ctype[2], "FREQ");
	hpar_ptr->crval[0] = md->RAcen;		/* degrees */
	hpar_ptr->crval[1] = 0.0; // SSG correction    md->DECcen;    /* degrees */
	hpar_ptr->crval[2] = md->fstart;  /* khz */   // for avg is favg (compute this = average of freq)
	
	hpar_ptr->crpix[0] = 0.5 + md->n1 / 2.0;	/* image center in pixels */
	hpar_ptr->crpix[1] = 0.5 + md->n2 / 2.0 - (md->DECcen / md->cellsize); // SSG correction for fits standard
	hpar_ptr->crpix[2] = 1.0; // fine for avg
	hpar_ptr->cdelt[0] = -md->cellsize;                     /* degrees */
	hpar_ptr->cdelt[1] = md->cellsize;                     /* degrees */
	//hpar_ptr->cdelt[2] = -172.032*1000000/MAX_CHANNELS;	//TODO: magic numbers				/* Mhz */
	hpar_ptr->cdelt[2] = md->df;
	
	hpar_ptr->equinox = 2000.0;		//year of coordinate system
	//hpar_ptr->epoch = 2004.???;		//year observations were taken

	sprintf (hpar_ptr->bunit, "Kelvin");
	sprintf (hpar_ptr->telescope, "Arecibo"); // added for VO compatibility
	sprintf (hpar_ptr->observer, "GALFACTS");
	sprintf (hpar_ptr->instrument, "ALFA");

	// for average images
	if ( (md->lowchan == 0) && (md->highchan == 1) )
	{
		FILE * BadChannelsFile; BadChannelsFile = fopen("BadChannels.list", "r");
		int *badchannels = (int*)calloc(MAX_CHANNELS, sizeof(int));

		int i = 0;
		int chan1, chan2;
		int chanCount = 0.0;
		float freqsum = 0.0;
		float fstart = md->fstart;
		float avgfreq = 0.0;
		int goodLowChan = 0;
		int goodHighChan = 0;
		int firstGoodChan = 1;

		if(BadChannelsFile != NULL)
		{
		int chan1, chan2;
		do
			{
			fscanf(BadChannelsFile, "%d %d\n", &chan1, &chan2);
			//printf("Excluded: chan1 = %d chan2 = %d\n", chan1, chan2);
			for(i=chan1; i<=chan2; i++)
				{
				badchannels[i] = 1;
				}
			}
		while(!feof(BadChannelsFile));
		fclose(BadChannelsFile);

		}

		i =0;

		for( i=md->avg_lowchan; i <= md->avg_highchan; i++ ) {
			if( badchannels[i] != 1 )
			{
				// record the first good low chan in case of channel exceptions
				if ( firstGoodChan == 1 ){
					goodLowChan = i;
					firstGoodChan = 0;
				}
				freqsum += fstart + (i * md->df);
				chanCount++;
			}
		}

		// find highest good channel for crval[2]
		firstGoodChan = 1;
		for( i=md->avg_highchan; i > md->avg_lowchan; i-- )
		{
			if( badchannels[i] != 1 )
			{
				// record the first good low chan in case of channel exceptions
				if ( firstGoodChan == 1 )
				{
					goodHighChan = i;
					firstGoodChan = 0;
				}
			}
		}


		avgfreq = freqsum / chanCount;

		hpar_ptr->crval[2] = avgfreq;
		hpar_ptr->cdelt[2] = md->df * (md->avg_highchan - md->avg_lowchan);  // for avg md->df * md->n3

		//printf("XXX For average channel got chanCount=%d fstart=%f freqsum=%f avgfreq=%f goodLowChan=%d goodHighChan=%d\n", chanCount, fstart, freqsum, avgfreq, goodLowChan, goodHighChan );

	}

}

void start_fits_cubes(const char * wapp, MapMetaData *md)
{
	char filename[64+1];

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


	snprintf (filename, 64, "images/%s_%04i_%04i_Icube.fits",md->title, md->lowchan, md->highchan-1 );
	fitsI = fopen(filename, "w");
	snprintf (filename, 64, "images/%s_%04i_%04i_Qcube.fits",md->title, md->lowchan, md->highchan-1 );
	fitsQ = fopen(filename, "w");
	snprintf (filename, 64, "images/%s_%04i_%04i_Ucube.fits",md->title, md->lowchan, md->highchan-1 );
	fitsU = fopen(filename, "w");
	snprintf (filename, 64, "images/%s_%04i_%04i_Vcube.fits",md->title, md->lowchan, md->highchan-1 );
	fitsV = fopen(filename, "w");
	snprintf (filename, 64, "images/%s_%04i_%04i_Weightcube.fits",md->title, md->lowchan, md->highchan-1 );
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



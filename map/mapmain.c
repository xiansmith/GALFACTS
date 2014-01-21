#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "common.h"
#include "fluxdata.h"
#include "balance.h"
#include "grid.h"
#include "map.h"
#include "scan.h"
#include "decdep.h"
#include "correctUV.h"
#include "pacorr.h"
#include "spec_corr.h"
// MAPMAIN with Chebyshev fitting

int multibeam; 
clock_t sclock;
//--------------------------------------------------------------------------------------------------------
static void print_usage(void)
{
	printf(
	"\n"
	"\twapp<n> - the beam number to create maps for (use multibeam for a multibeam map)\n"
	"\tfcenter - the center frequency of this wapp\n"
	"\tlowchan, highchan - the channel range to use\n"
	"\tramin, ramax - the Right assension range in degrees\n"
	"\tdecmin, decmax - the Declination range in degrees\n"
	"\tcell size - pixel size of map in arc minutes\n"
	"\tpatch radius - area in pixels to calculate point spread response (recommended 8)\n"
	"\tgridtype - upscans 1 downscans 2 both 3\n"
	"\tbalance day order - order of the polynomial fit in day-to-day balancing\n"
	"\tbalance scan order - order of the polynomial fit in scan-by-scan balancing\n"
	"\tbalance loopgain - the gain of each balance iteration (should be between 0 and 1)\n"
	"\tbalance epsilon - the percent change of simasq threshold where the iteraiton will stop\n"
	"\tshowprogress - 1 to create a basket weave progress cube, 0 otherwise\n"
	"\n");
}
//-------------------------------------------------------------------------------
static void create_fits_cube(FluxWappData * wappdata, char * wapp, MapMetaData * md)
{
	float *dataI, *dataQ, *dataU, *dataV, *weight; 
	int chan, n;
	int numbytes;

	printf("Map size: %d x %d\n", md->n1, md->n2);
	printf("DEC range (degrees): (%f %f)\n", md->decmin, md->decmax);
	printf("RA range (degrees): (%f %f)\n", md->ramin, md->ramax);
	printf("Map centre: %.2f %.2f\n", md->RAcen, md->DECcen); 
	printf("Channel range: (%i, %i]\n", md->lowchan, md->highchan);
	printf("Channel count: %i\n", md->n3);
	printf("Patch radius: %i\n", md->patch);
	printf("Beamwidths: %i\n", md->beamwidths);
	printf("Cell Size: %f degrees\n", md->cellsize);
	printf("fwhm: %f arcmin\n", md->fwhm);
	printf("Grid Type: %i\n", md->gridtype);
	printf("Center Frequency: %f\n", md->fcen); 
	printf("Start Frequency: %f\n", md->fstart);
 	printf("Balance Loop Gain: %f\n", md->balgain);
 	printf("Balance Epsilon: %f\n", md->balepsilon);
	printf("FWHM: %f\n", md->fwhm);

	numbytes = md->n1 * md->n2 * sizeof (float);
	dataI   = (float *) malloc (numbytes);
	dataQ   = (float *) malloc (numbytes);
	dataU   = (float *) malloc (numbytes);
	dataV   = (float *) malloc (numbytes);
	weight  = (float *) malloc (numbytes);
	
	if (!dataI || !dataQ || !dataU || !dataV || !weight){printf("ERROR: memory allocation of %i bytes failed !\n", numbytes); return;}

	init_psf_lookup_table(65537, 9.1, md->fwhm);
	init_psf_map(md->fwhm, md->cellsize, md->beamwidths);
	

// Read the gain calibration data
// GALFACTS_calibration.dat should be in the directory from where the script starts
// If GALFACTS_calibration.dat cannot be read then the gain correction is done by 
// the hard wired values from "beam_gain_calibration()"
// If GALFACTS_calibration.dat can be read (is present) then the calibration is done 
// with the values from the table using "beam_gain_calibration_table()"
// If chan<cal_low or chan>cal_high then the calibration uses the average values from the 
// first two lines of the table
	
	float cal_table[4098][7], gmean[7], cal_table_avg[md->n3+2][7];
	int i, j, cal_low, cal_high, cal_flag=0, Nc;

	FILE *cal_file;
    cal_file = fopen("GALFACTS_calibration.dat", "r");
    if(cal_file == NULL)
		{
		printf("Failed to open GALFACTS_calibration.dat\n"); cal_flag = 0;
		}
		else
			{
			cal_flag = 1;
			for(i=0; i<7; i++) fscanf(cal_file, "%f", &cal_table[0][i]); 
			for(i=0; i<7; i++) cal_table_avg[0][i] = cal_table[0][i]; 
			for(i=0; i<7; i++) fscanf(cal_file, "%f", &cal_table[1][i]); 
			fscanf(cal_file, "%d %d", &cal_low, &cal_high);
			for(j=0; j < cal_high - cal_low; j++)
				{
				for(i=0; i<7; i++) fscanf(cal_file, "%f", &cal_table[j+2][i]); 
				}		
			fclose(cal_file);			
			for(i=0; i<7; i++) cal_table[1][i] = 0;  Nc = 0;
			for(j=0; j < cal_high - cal_low; j++)
				{
				if(cal_low+j >= md->avg_lowchan && cal_low+j< md->avg_highchan)
					{
					for(i=0; i<7; i++) cal_table[1][i] += cal_table[j+2][i]; Nc++;
					}
				}					
			for(i=0; i<7; i++) cal_table[1][i] /= Nc;
			

			if(md->avg > 0)
				{
				for(j=2; j < cal_high - cal_low + 2 - md->avg; j++)
					{
					for(i=0; i<7; i++)
						{
						float s = 0.0;
						for(n=j; n<j+md->avg; n++)
							{
							s += cal_table[n][i];
							}
						cal_table[j][i] = s/md->avg;
						}
					}
				}
			}
//
	
	start_fits_cubes(wapp, md);
	
	float *cIc = (float*)malloc(wappdata->numDays*(md->dec_order+1) * sizeof(float));
	float *cQc = (float*)malloc(wappdata->numDays*(md->dec_order+1) * sizeof(float));
	float *cUc = (float*)malloc(wappdata->numDays*(md->dec_order+1) * sizeof(float));
	float *cVc = (float*)malloc(wappdata->numDays*(md->dec_order+1) * sizeof(float));
	
		chan = md->lowchan;
		while(chan < md->highchan)
			{
			printf("Channel: %i \n", chan);
			printf("Reading data ...\n"); 

			fluxwappdata_readchan_binary(md->field, md->band, wappdata, chan, BASKETWEAVE, md->avg, md->decmin, md->decmax);

			printf("Apply the feed coupling mueller matrix correction...\n");
			correct_UV(wappdata,chan,md);
			
			printf("Removing DEC dependence...\n"); 
			calculate_dec_dependence(wappdata, md->dec_order, chan, cIc, cQc, cUc, cVc, md->avg);
			
			printf("Apply the position angle correction...\n");
			pa_corr(wappdata,chan,md);

			printf("Apply the spectral Index correction...\n");
			spec_corr(wappdata,chan,md);

			printf("Beam gain calibration...\n"); 		
			if(cal_flag) beam_gain_calibration_table(wappdata, cal_low, cal_high, cal_table, chan); 
		
			printf("Determine scan lines ...\n"); 
			determine_scan_lines(wappdata, md->decmin, md->decmax);	
			printf("Finding crossing points ...\n"); 
			find_intersections(wappdata);							
			
			printf("Performing basketweaving...\n"); 
			balance_data(wappdata, md->day_iter, md->scan_iter, md->balgain, md->balepsilon, md->bw_order);
		
			//printf("Writing basketweaved time series ...\n");
			//fluxwappdata_writechan_binary(wappdata,chan);
			//fluxwappdata_writechan_binary_single(wappdata,chan);

			printf("Gridding ...\n"); 
			grid_data(wappdata, md, dataI, dataQ, dataU, dataV, weight);
			printf("Writing fits data ...\n"); 
			write_fits_planes(dataI, dataQ, dataU, dataV, weight);
			if(md->avg == 0) chan++; else chan += md->avg;
			
			for(i=0; i<wappdata->numDays; i++) 
				{
				free(wappdata->daydata[i].records);
				}
			}
	
	//printf("Finishing fits files\n");
	finish_fits_cubes();

	//free memory
	free(cIc); 
	free(cQc); 
	free(cUc); 
	free(cVc);
		
	free(dataI);
	free(dataQ);
	free(dataU);
	free(dataV);
	free(weight);

	printf("Done!\n"); 
}
//-------------------------------------------------------------------------------
int main(int argc, char * argv[])
{
// Hydrogen freq = 1420.4057, Hchan = 2752.63
clock_t time0 = clock();
FluxWappData *wappdata;
MapMetaData md;
char ** files;
int numDays;
/* Inputs to this program */
//TODO: it would be nice to make these defaults
//at the moment the values have no effect
char * wapp;
int pfit_type;
float pfit_lambda;
int day_iter = 0;
int scan_iter = 0;

/* Process command line arguments */ 
/* Convert args into more usable units were required */

if(argc != 25){printf("Usage: %s <parameters_list>\n", argv[0]); return EXIT_FAILURE;}

wapp = argv[1];
md.fcen = (float) atof(argv[2]);
md.lowchan = atoi(argv[3]);
md.highchan = atoi(argv[4]);
md.ramin = (float) atof(argv[5]); //convert to degrees ?
md.ramax = (float) atof(argv[6]); //convert to degrees ?
md.decmin = (float) atof(argv[7]);
md.decmax = (float) atof(argv[8]);
md.cellsize = (float) atof(argv[9]) / 60.0; //convert to degrees
md.patch = atoi(argv[10]); //TODO: remove this
md.beamwidths = atoi(argv[10]);
md.gridtype = atoi(argv[11]);
md.balgain = (float) atof(argv[12]);
md.balepsilon = (float) atof(argv[13]);
md.bw_order = atoi(argv[14]);
md.dec_order = atoi(argv[15]);
md.avg = atoi(argv[16]); 
md.avg_lowchan = atoi(argv[17]);
md.avg_highchan = atoi(argv[18]);
md.title = argv[19];
strcpy( md.field, argv[20]);
md.band = atoi(argv[21]);
md.fwhm = atoi(argv[22]); //2.0; //TODO: paramaterize this
md.day_iter = atoi(argv[23]);
md.scan_iter = atoi(argv[24]);



md.RAcen = (md.ramax + md.ramin)/2.0;
md.DECcen = (md.decmax + md.decmin)/2.0;
md.RArange =  md.ramax - md.ramin;
md.DECrange = md.decmax - md.decmin;
md.n1 = (int)(md.RArange/md.cellsize) + 1;
md.n2 = (int)(md.DECrange/md.cellsize) + 1;
if(md.avg == 0)
	{
	md.n3 = md.highchan - md.lowchan; 
	md.fstart = 0.5*(md.fcen - (md.avg_lowchan-(MAX_CHANNELS/2-1))*0.042 + md.fcen - (md.avg_highchan-(MAX_CHANNELS/2-1))*0.042)*1000000;
	//md.fstart = 0.5*(md.fcen + (md.avg_lowchan-(MAX_CHANNELS/2-1))*0.042 + md.fcen - (md.avg_highchan-(MAX_CHANNELS/2-1))*0.042)*1000000; // old calc
	//md.fstart = (md.fcen - (((MAX_CHANNELS/2-1) - 0.5 * (md.avg_highchan + md.avg_lowchan)) *0.042) * 1000000);   // average channel centre freq in Hz
	md.df = (md.avg_highchan - md.avg_lowchan)*42000; 
	}
else
	{
	md.n3 = (md.highchan - md.lowchan)/md.avg; 
	md.fstart = (md.fcen - (md.lowchan-(MAX_CHANNELS/2-1))*0.042)*1000000;//watchout for the sign for MOCK needs to change for WAPP
	md.df = -md.avg*172032000.0/MAX_CHANNELS; 
	}
		
numDays = get_date_dirs("./", &files);
if(numDays <= 0){printf("ERROR: could not find any date dirs\n"); return EXIT_FAILURE;}

if(!strcmp(wapp,"multibeam"))
	{
	numDays = numDays * 7;
	multibeam = 1;
	}
else multibeam = 0;
	
// allocate and initialize the wapp data days
wappdata = fluxwappdata_alloc(wapp, files, numDays);

printf("Creating a frequency cube\n");
create_fits_cube(wappdata, wapp, &md );
	
free(files);

printf("Total computation time = %f\n", ((double)clock() - time0)/CLOCKS_PER_SEC);

return EXIT_SUCCESS;
}
//-------------------------------------------------------------------------------


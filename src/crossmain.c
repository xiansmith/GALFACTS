#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "fluxdata.h"
//#include "jsd_futil.h"
#include "balance.h"
#include "grid.h"
#include "map.h"
#include "scan.h"
#include "decdependence.h"

int multibeam; //SSG
static void print_usage(const char * prog)
{

}


static void create_fits_cube(FluxWappData * wappdata, char * wapp, MapMetaData * md, int day_order, int scan_order, float balgain, float balepsilon, int show_progress)
{
	float *dataI, *dataQ, *dataU, *dataV, *weight; 
	int chan;
	int numbytes;


	//read channel data files for all the days
	printf("reading data ...\n");
	//scanwappdata_readchan(wappdata, wapp, chan, days, num_days);
	fluxwappdata_readchan(wappdata, chan);

	//determine scan lines
	printf("determine scan lines ...\n");
	determine_scan_lines(wappdata, md->decmin, md->decmax);

	//write out the balanced data
	printf("writing balanced data to file ...\n");
	fluxwappdata_writechan(wappdata, chan);


	printf("done!\n");
}

int main(int argc, char * argv[])
{
	FluxWappData *wappdata;
	MapMetaData md;

	char ** files;
	int numDays;
	int show_progress;//ssg
	/* Inputs to this program */
	//TODO: it would be nice to make these defaults
	//at the moment the values have no effect

	char * wapp;
	int day_order;
	int scan_order;
	float balgain;
	float balepsilon;
	int showprogress;

	/* Process command line arguments */ 
	/* Convert args into more usable units were required */

	if (argc != 18) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	} else {
		wapp = argv[1];
		md.fcen = (float) atof(argv[2]);
		md.lowchan = atoi(argv[3]);
		md.highchan = atoi(argv[4]);
		md.ramin = (float) atof(argv[5]) * 15; //convert to degrees
		md.ramax = (float) atof(argv[6]) * 15; //convert to degrees
		md.decmin = (float) atof(argv[7]);
		md.decmax = (float) atof(argv[8]);
		md.cellsize = (float) atof(argv[9]) / 60; //convert to degrees
		md.patch = atoi(argv[10]);
		md.gridtype = atoi(argv[11]);
		day_order = atoi(argv[12]);
		scan_order = atoi(argv[13]);
		balgain = (float) atof(argv[14]);
		balepsilon = (float) atof(argv[15]);
		showprogress = atoi(argv[16]);
		md.title = argv[17];

		md.RAcen = (md.ramax + md.ramin)/2.0;
		md.DECcen = (md.decmax + md.decmin)/2.0;
		md.RArange =  md.ramax - md.ramin;
		md.DECrange = md.decmax - md.decmin;
		md.n1 = (int)(md.RArange/md.cellsize) + 1;
		md.n2 = (int)(md.DECrange/md.cellsize) + 1;
		md.n3 = md.highchan - md.lowchan;
		md.fstart = md.fcen + (md.lowchan-127) * (100.0/256.0);
		
	}
	numDays = get_date_dirs("./", &files);
	if (numDays <= 0) {
		printf("ERROR: could not find any date dirs\n");
		return EXIT_FAILURE;
	}
	//SSG
	if (!strcmp(wapp,"beam8"))
	{
		numDays = numDays * 7;
		multibeam = 1;
	}
	else
		multibeam = 0;
	//SSG
	// allocate and initialize the wapp data days
	wappdata = fluxwappdata_alloc(wapp, files, numDays);
//	printf("After wapp alloc");
	printf("Creating a frequency cube\n");
	create_fits_cube(wappdata, wapp, &md, day_order, scan_order, balgain, balepsilon, showprogress);
	
	//printf("creating a avg map\n");
	//create_fits_map(wappdata, &md, fcen, patch, balorder, balgain, balepsilon, title);	
	
	free(files);
	//fluxwappdata_free(wappdata);
	//free(wappdata);

	return EXIT_SUCCESS;
}



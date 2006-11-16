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

static void print_usage(const char * prog)
{
	printf(
	"Usage: %s wapp<n> <fcenter> <lowchan> <highchan> <ramin> <ramax> "
	"<decmin> <decmax> <cell size> <patch radius> <balance order> <balance ref>\n"
	"\n"
	"\twapp<n> - the wapp number to create maps for\n"
	"\tfcenter - the center frequency of this wapp\n"
	"\tlowchan, highchan - the channel range to use\n"
	"\tramin, ramax - the Right assension range in hours\n"
	"\tdecmin, decmax - the Declination range in degrees\n"
	"\tcell size - cell size of map in arc minutes\n"
	"\tpatch radius - area in pixels to calculate point spread response\n"
	"\tbalance order - order of the polynomial fit in balancing\n"
	"\tbalance ref - the MJD of the day to use as the reference day\n"
	"\n"
	"eg: %s wapp1 1170.0 25 230 6.75 7.75 11.0 12.1 0.25 8 5 53108\n"
	"\n"
	"Call this program through the following names (hint: use softlinks)\n"
	"\tmapcube - to create cubes, on slice per channel\n"
	"\tmapavg - to create combined maps of all channels\n"
	"\n", prog, prog); 
}


static void create_fits_map(FluxWappData * wappdata, MapMetaData * md, float fcen, int patch, int balorder, char * balrefmjd)
{
	int n1, n2;
	float *dataI, *dataQ, *dataU, *dataV; 

	n1 = (int)(md->RArange/md->cellsize) + 1;
	n2 = (int)(md->DECrange/md->cellsize) + 1;

	printf("Map size: %d x %d\n", n1, n2);
	printf("Map centre: %.2f %.2f\n", md->RAcen, md->DECcen); 
	printf("Cell size: %7.5f\n", md->cellsize);
	printf("Channel range: (%i, %i]\n", md->lowchan, md->highchan);
	printf("Patch radius: %i\n", patch);
	printf("Center Frequency: %g\n", fcen); 

	dataI  = (float *) malloc (n1 * n2 * sizeof (float));
	dataQ  = (float *) malloc (n1 * n2 * sizeof (float));
	dataU  = (float *) malloc (n1 * n2 * sizeof (float));
	dataV  = (float *) malloc (n1 * n2 * sizeof (float));

	//read channel data files for all the days
	printf("reading data ...\n");
	fluxwappdata_readavg(wappdata);

	//perform balancing
	printf("performing balancing ...\n");
	balance_data(wappdata, balorder, balrefmjd);

	//write out the balanced data
	printf("writing balanced data to file ...\n");
	fluxwappdata_writeavg(wappdata);

	//perform gridding
	printf("performing gridding ...\n");
	grid_data(wappdata, md, dataI, dataQ, dataU, dataV, n1, n2, patch);

	printf("writing fits files ...");
	write_fits_maps(wappdata->wapp, dataI, dataQ, dataU, dataV, n1, n2, md->RAcen, md->DECcen, md->cellsize);

	//free working memory
	free(dataV);
	free(dataI);
	free(dataQ);
	free(dataU);

	printf("done!\n");

}


static void create_fits_cube(FluxWappData * wappdata, MapMetaData * md, float fcen, int patch, int balorder, char * balrefmjd)
{
	int n1, n2, n3;
	float *dataI, *dataQ, *dataU, *dataV; 
	float fstart;
	int chan;

	n1 = (int)(md->RArange/md->cellsize) + 1;
	n2 = (int)(md->DECrange/md->cellsize) + 1;
	n3 = md->highchan - md->lowchan;
	fstart = fcen + (md->lowchan-127) * (100.0/256.0);

	printf("Map size: %d x %d\n", n1, n2);
	printf("Map centre: %.2f %.2f\n", md->RAcen, md->DECcen); 
	printf("Cell size: %7.5f\n", md->cellsize);
	printf("Channel range: (%i, %i]\n", md->lowchan, md->highchan);
	printf("Channel count: %i\n", n3);
	printf("Patch radius: %i\n", patch);
	printf("Center Frequency: %g\n", fcen); 
	printf("Start Frequency: %g\n", fstart);

	dataI  = (float *) malloc (n1 * n2 * sizeof (float));
	dataQ  = (float *) malloc (n1 * n2 * sizeof (float));
	dataU  = (float *) malloc (n1 * n2 * sizeof (float));
	dataV  = (float *) malloc (n1 * n2 * sizeof (float));

	start_fits_cubes(wappdata->wapp, n1, n2, n3, md->RAcen, md->DECcen, fstart, md->cellsize);

	for (chan=md->lowchan; chan<md->highchan; chan++)
	{
		printf("Channel: %i\n", chan);

		//read channel data files for all the days
		printf("reading data ...\n");
		fluxwappdata_readchan(wappdata, chan);

		//perform balancing
		printf("performing balancing ...\n");
		balance_data(wappdata, balorder, balrefmjd);

		//write out the balanced data
		printf("writing balanced data to file ...\n");
		fluxwappdata_writechan(wappdata, chan);

		//perform gridding
		printf("performing gridding ...\n");
		grid_data(wappdata, md, dataI, dataQ, dataU, dataV, n1, n2, patch);

		//write fits data
		printf("writing fits data ...\n");
		write_fits_planes(dataI, dataQ, dataU, dataV);
	}

	//free memory
	free(dataI);
	free(dataQ);
	free(dataU);
	free(dataV);
	
	printf("finishing fits files ...");
	finish_fits_cubes();

	printf("done!\n");
}

int main(int argc, char * argv[])
{
	FluxWappData * wappdata;
	MapMetaData md;

	char ** files;
	int numDays;
	
	/* Inputs to this program */
	//TODO: it would be nice to make these defaults
	//at the moment the values have no effect

	char * prog;
	char * wapp;
	int patch;
	float fcen;
	int balorder;
	char * balref;

	/* Process command line arguments */ 
	/* Convert args into more usable units were required */
	
	if (argc != 13) {
printf("argc: %i\n", argc);
		print_usage(argv[0]);
		return EXIT_FAILURE;
	} else {
		prog = argv[0];
		if (!strstr(prog, "mapcube") && !strstr(prog, "mapavg")) {
			print_usage("[mapcube|mapavg]");
			return EXIT_FAILURE;
		}
		wapp = argv[1];
		fcen = (float) atof(argv[2]);
		md.lowchan = atoi(argv[3]);
		md.highchan = atoi(argv[4]);
		md.ramin = (float) atof(argv[5]) * 15; //convert to degrees
		md.ramax = (float) atof(argv[6]) * 15; //convert to degrees
		md.decmin = (float) atof(argv[7]);
		md.decmax = (float) atof(argv[8]);
		md.cellsize = (float) atof(argv[9]) / 60; //convert to degrees
		patch = atoi(argv[10]);
		balorder = atoi(argv[11]);
		balref = argv[12];

		md.RAcen = (md.ramax + md.ramin)/2.0;
		md.DECcen = (md.decmax + md.decmin)/2.0;
		md.RArange =  md.ramax - md.ramin;
		md.DECrange = md.decmax - md.decmin;

	}
	
	
	numDays = get_date_dirs("./", &files);
	if (numDays <= 0) {
		printf("ERROR: could not find any date dirs\n");
		return EXIT_FAILURE;
	}


	// allocate and initialize the wapp data days
	wappdata = fluxwappdata_alloc(wapp, files, numDays);


	if (strstr(prog, "mapcube")) {
		printf("creating a frequency cube\n");
		create_fits_cube(wappdata, &md, fcen, patch, balorder, balref);	
	} else if (strstr(prog, "mapavg")) {
		printf("creating band averaged maps\n");
		create_fits_map(wappdata, &md, fcen, patch, balorder, balref);	
	}
	
	free(files);
	fluxwappdata_free(wappdata);

	return EXIT_SUCCESS;
}



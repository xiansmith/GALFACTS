#include "common.h"
#include "spec.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>
#include <glob.h>

static void process_dataset(const char * datadirname, const char * datedir, const char * beamdir)
{	
	
	FILE * datafile, *outfile;
	int numRecords;
	int n, chan;
	int count;
	double linear = 0.0;

	char datafilename[128+1];
	char outfilename[128+1];
	char globpattern[128+1];
	glob_t globbuf;

	/* Open Datafile */
	globbuf.gl_offs = 1;
	sprintf(globpattern, "%s/%s/*.za_scan*.%s.%s.spec", datadirname, datedir, beamdir, datedir);
	glob(globpattern, 0, NULL, &globbuf);
	if (globbuf.gl_pathc < 1) {
		printf("ERROR: unable to find file with pattern %s\n", globpattern);
	}
	strncpy(datafilename, globbuf.gl_pathv[0], 128);
	//globfree(&globbuf);	

	if ( (datafile = fopen(datafilename, "rb") ) == NULL )
	{ 
		printf("ERROR: can't open data file '%s'\n", datafilename); 
		return;
	}

	sprintf(outfilename, "%s.hack", datafilename);
	if ( (outfile = fopen(outfilename, "wb") ) == NULL )
	{ 
		printf("ERROR: can't open data file '%s'\n", outfilename); 
		return;
	}

	//loop
	count = 0;
	while (1)
	{
		SpecPointingBlock block;
		PolSet calon;
		PolSet caloff;
		int itemsRead = 0;
	
		//read pointing block
		itemsRead = fread(&block, sizeof(SpecPointingBlock), 1, datafile);

		//read data
		itemsRead += fread(&(calon),  sizeof(PolSet), 1, datafile);
		itemsRead += fread(&(caloff), sizeof(PolSet), 1, datafile);

		if (itemsRead != 3) break;

		linear += 0.0002;
		//hack the spec file values
		for (chan=0; chan<MAX_CHANNELS; chan++) 
		{
			float val = 1.0;
			float cal = 0.1;
			if (strcmp(datedir, "53356") == 0) {
				val = 1.2;
			}
			if (strcmp(datedir, "53360") == 0) {
				val = 0.8;
			}

			if (strcmp(datedir, "53365") == 0) {
				val += linear;
			}
			caloff.xx[chan] = val;
			calon.xx[chan] = val+cal;
			caloff.xy[chan] = val;
			calon.xy[chan] = val+cal;
			caloff.yx[chan] = val;
			calon.yx[chan] = val+cal;
			caloff.yy[chan] = val+1;
			calon.yy[chan] = val+cal+1;
		}

		//write to the output
		itemsRead = fwrite(&block, sizeof(SpecPointingBlock), 1, outfile);
		itemsRead += fwrite(&calon,  sizeof(PolSet), 1, outfile);
		itemsRead += fwrite(&caloff, sizeof(PolSet), 1, outfile);

		count++;
	}

	fclose(datafile);
	fclose(outfile);

	/* done */
	printf("hacked %i records\n", count);

	return;
}


int main(int argc, char *argv[])
{
	int numdirs;
	char * datadirname;
	char ** datedirs;
	int mjd;
	char *beam;

	if (argc != 3) {
		printf("Usage: %s <datadir> <beam>\n", argv[0]);
		printf("eg: %s /n/swift2/galfacts/data/A1863/SPEC beam0\n", argv[0]);
		return EXIT_FAILURE;
	} 
	else 
	{
		datadirname = argv[1];
		beam = argv[2];
	}

	numdirs = get_date_dirs(datadirname, &datedirs);
	printf("Found %i data directories in %s\n", numdirs, datadirname);

	for (mjd=0; mjd<numdirs; mjd++) 
	{
		const char * datedir = datedirs[mjd];
		mode_t mode = S_IRWXU|S_IRWXG|S_IRWXO;
		mkdir(datedir, mode);
		chdir(datedir);
		mkdir(beam, mode);
		chdir(beam);
		process_dataset(datadirname, datedir, beam);
		chdir("..");
		chdir("..");
	}

	//TODO: free datedirs
	return EXIT_SUCCESS;
}


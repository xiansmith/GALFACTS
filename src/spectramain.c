#include "common.h"
#include "spec.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

static void write_spectra_files(SpecRecord dataset[], int numRecords, const char * date, const char * beam)
{
	int i;
	int chan;
	FILE *file;
	char filename[32+1];

	mode_t mode = S_IRWXU|S_IRWXG|S_IRWXO;
	mkdir("spectra", mode);
	chdir("spectra");
	mkdir(beam, mode);
	chdir(beam);
	mkdir(date, mode);
	chdir(date);
	
	for (i=0; i<numRecords; i++) 
	{
		SpecRecord *pRec = &dataset[i];

		sprintf(filename, "spectra%05i.dat", i);
		file = fopen(filename, "w");
		if (file == NULL) {
			printf("ERROR: failed to open file %s\n", filename);
		}


		fprintf(file, "#AST: %.2f    RA: %.6f   DEC: %.6f\n", pRec->AST, pRec->RA, pRec->DEC);
		fprintf(file, "#chan calonxx caloffxx calonyy caloffyy\n");

		for (chan=0; chan < 256; chan++) {
			fprintf(file, "%3i %.6f %.6f %.6f %.6f\n", chan, pRec->calon.xx[chan], pRec->caloff.xx[chan], pRec->calon.yy[chan], pRec->caloff.yy[chan]);
		}
		fclose(file);
	}
	chdir("../");
	chdir("../");
	chdir("../");

}
int main(int argc, char *argv[])
{
	char datafilename[128+1];
	char * datadir;
	int beam;
	char * position;
	char * date;
	int numRecords;
	FILE *datafile;
	SpecRecord * dataset;

	if (argc != 5) {
		printf("usage: %s <datadir> <position> <beam> <day>\n", argv[0]);
		return EXIT_FAILURE;
	} else {
		datadir = argv[1];
		position = argv[2];
		beam = atoi(argv[3]);
		date = argv[4];
	}
	sprintf(datafilename, "%s/%s/A1947.za_scan_%s.beam%i.%s.spec", datadir, date, position, beam, date);

		/* Read Datafile */
	if ( (datafile = fopen(datafilename, "rb") ) == NULL )
	{ 
		printf("ERROR: can't open data file '%s'\n", datafilename); 
		return EXIT_FAILURE;
	}

	printf("reading data file %s ... \n", datafilename);
	numRecords = read_datafile(datafile, &dataset, beam);
	fclose(datafile);
	printf("Read %i records\n", numRecords);

	if (numRecords <= 0) {
		printf("ERROR: there are no records!\n");
		return EXIT_FAILURE;
	}

	write_spectra_files(dataset, numRecords, date, beam);

	return EXIT_SUCCESS;
}

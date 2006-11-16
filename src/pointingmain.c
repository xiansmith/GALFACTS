#include "common.h"
#include <stdlib.h>
#include <stdio.h>



int main(int argc, char * argv[])
{
	char * filename;
	char outfilename[128+1];
	SpecRecord * dataset;
	FILE * datafile, * outfile;
	int i, n, size;
		


	for (i=1; i<argc; i++) 
	{
		filename = argv[i];

		datafile = fopen(filename, "r");
		if (datafile == NULL) {
			printf("ERROR: unable to open file %s\n", filename);
			continue;
		}

		/* read data file */
		size = read_datafile(datafile, &dataset);
		fclose(datafile);

		snprintf(outfilename, 128, "%s.pointing.dat", filename);

		outfile = fopen(outfilename, "w");
		if (outfile == NULL) {
			printf("ERROR: unable to open file %s\n", outfilename);
			continue;
		}

		fprintf(outfile, "# %s\n", filename);
		for (n=0; n<size; n++) {
			fprintf(outfile, "%8.6f %8.6f %7.2f\n", dataset[n].RA, dataset[n].DEC, dataset[n].AST);
		}


		free(dataset);
	}
	return EXIT_SUCCESS;
}


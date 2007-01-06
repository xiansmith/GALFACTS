#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "markdata.h"


/*
 * Reads the specified badfilename and marks entries in the given dataset
 * as bad (ie: flagBAD = 1) for the range in the datafile.
 * 
 * bad_datapoints.dat file format is:
 * # comment
 * blank lines allowed
 * <type> <low> <high>
 * type can be one of 'RA' 'DEC' or 'AST'
 * low and high determine the ranges to flag for the given type
 * negative values and inf allowed
 */
void mark_bad_datapoints(const char * badfilename, SpecRecord dataset[], int numRecords)
{
	int i;
	SpecRecord * pRec;
	FILE * badfile;
	char buf[BADFILE_LINE_LEN];
	float low, high;
	char type[8];
	
	badfile = fopen(badfilename, "r");
	if (badfile == NULL) {
		printf("INFO: No manually flagged bad datapoints file.\n");
		return;
	}

	//read line from file
	while (fgets(buf, 80, badfile) != NULL) 
	{
		if (buf[0] == '#') continue;
		if (buf[0] == '\n') continue;

		sscanf(buf, "%s %f %f", type, &low, &high);
		printf("bad datapoints in: %s (%f, %f)\n", type, low, high);

		if (strcmp("RA", type) == 0) {
			for (i=0; i<numRecords; i++) {
				pRec = &(dataset[i]);	
				if ((pRec->RA > low) && (pRec->RA < high)) {
					pRec->flagBAD = 1;
				}
			}
		}
		else if (strcmp("DEC", type) == 0) {
			for (i=0; i<numRecords; i++) {
				pRec = &(dataset[i]);	
				if ((pRec->DEC > low) && (pRec->DEC < high)) {
					pRec->flagBAD = 1;
				}
			}
		}
		else if (strcmp("AST", type) == 0) {
			for (i=0; i<numRecords; i++) {
				pRec = &(dataset[i]);	
				if ((pRec->AST > low) && (pRec->AST < high)) {
					pRec->flagBAD = 1;
				}
			}
		}
		else {
			printf("WARN: unrecognized bad type '%s'\n", type);
		}
	}

	fclose(badfile);
}



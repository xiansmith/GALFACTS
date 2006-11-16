#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "markdata.h"



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



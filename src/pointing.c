#include "pointing.h"
#include <stdlib.h>
#include <stdio.h>

void create_annotations(SpecRecord dataset[], int size)
{
	int n;
	FILE * annfile;
	const char * annfilename = "pointing.ann";
	double h, secderiv;
	int found;

	annfile = fopen(annfilename, "w");
	if (annfile == NULL) {
		printf("ERROR: unable to open '%s' for writing\n", annfilename);
		return;
	}

	fprintf(annfile, "#Annotations\n");

	found = 0;
	fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", dataset[0].RA, dataset[0].DEC, dataset[0].AST);
	for (n=1; n<size-1; n++)
	{
		if (dataset[n].flagBAD) {
			fprintf(annfile, "COLOUR %s\n", "RED"); 
		} else {
			fprintf(annfile, "COLOUR %s\n", "GREEN");
		}
		fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", 
			dataset[n].RA, dataset[n].DEC, dataset[n].AST);

		//write a red circle around pointing rages of change that are too high
		h = (dataset[n+1].AST - dataset[n-1].AST) / 2; //should be 0.2 
		secderiv = (dataset[n-1].DEC - 2*dataset[n].DEC + dataset[n+1].DEC) / (h*h); //second derivitive
		if (fabs(secderiv) > 0.08) {
			if (!found) {
				found = 1; 
				fprintf(annfile, "COLOUR %s\n", "RED"); 
				fprintf(annfile, "CIRCLE W %7.4f %7.4f %7.4f #%7.2f\n", 
					dataset[n].RA, dataset[n].DEC, 0.025, dataset[n].AST);
			}
		} else {
			found = 0;
		}
	}
	fprintf(annfile, "DOT W %7.4f %7.4f #%7.2f\n", dataset[n].RA, dataset[n].DEC, dataset[n].AST);
	
	fclose(annfile);
}


#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include "jsd/jsd_futil.h"
//this programs reads a standard data file (fluxtime,balance,clean,beammodel) with 7 fields per line and writes out a new file with 3 fields (RA,DEC,Stokes_I)

int main(int argc, char *argv[])
{
	char * outputfilename, *inputfilename;
	int j,linecount;
	FILE * inputfile, *outputfile;
	float RA,DEC,AST;
	double I,Q,U,V; 
	char header[81];
	if(argc != 3)
	{
		printf("Usage: convert <inputfilename> <outputfilename>\n");
		return EXIT_FAILURE;
	}
	else
	{
		inputfilename = argv[1];
		outputfilename = argv[2];
	}
	inputfile = fopen(inputfilename, "r");
	if(inputfile == NULL)
	{	
		printf("Error: Cannot open input file.\n");
		exit(EXIT_FAILURE);
	}	
	outputfile = fopen(outputfilename, "w");
	if(outputfile == NULL)
	{	
		printf("Error: Cannot open output file.\n");
		exit(EXIT_FAILURE);
	}	
        fgets(header, 80, inputfile);
	linecount = jsd_line_count(inputfile);
	fprintf(outputfile,"#RA DEC I\n");
	printf("Converting files...\n");
	for(j = 0;j < linecount;j++)
	{
		fscanf(inputfile,"%f %f %f %lf %lf %lf %lf",&RA,&DEC,&AST,&I,&Q,&U,&V);
		fprintf(outputfile,"%f %f %lf\n",RA,DEC,I);
	}
	fclose(outputfile);
	fclose(inputfile);
	printf("Done.\n");
	return(1);
}

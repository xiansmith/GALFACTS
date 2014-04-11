#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "jsd/jsd_futil.h"
//#include "jsd/jsd_fit.h"

void leakage_pattern(char *,char *,int,int,float,FILE*,float,float);

int main(int argc,char *argv[]) {
	//char *mjd;
	int i,beam;
        char ** files;
	char *subdir;
        int numDays;
	int chan;
	float RAh,RAm,RAs,DECd,DECm,DECs,RA,DEC,radius;
	subdir = argv[1];
	chan = atoi(argv[2]);
	RAh = atof(argv[3]);
	RAm = atof(argv[4]);
	RAs = atof(argv[5]);
	DECd = atof(argv[6]);
	DECm = atof(argv[7]);
	DECs = atof(argv[8]);
	radius = atof(argv[9]);
	
	RA = (RAh+RAm/60.0+RAs/3600.0)*15.0;
	DEC = DECd+DECm/60.0+DECs/3600.0;

	printf("location: %f %f\n",RA,DEC);

        numDays = get_date_dirs(subdir, &files);
        if (numDays <= 0) {
                printf("ERROR: could not find any date dirs\n");
                return EXIT_FAILURE;
        }
//    	mjd = argv[1];
//    	beam = atoi(argv[2]);
    	//agg = atoi(argv[3]);

  	char outfilename[81];
	FILE *outfile;
	
	for(beam=0;beam<7;beam++)
	{
		printf("\nBEAM %d\n\n",beam);
		sprintf(outfilename,"loc_%02.0f%02.0f%02.2f_%02.0f%02.0f%02.2f_b%d_%04i.dat",RAh,RAm,RAs,DECd,DECm,DECs,beam,chan);
		outfile = fopen(outfilename,"w");
		fprintf(outfile,"#dRA dDEC I Q U V\n");
		for(i=0;i<numDays;i++)
		    	leakage_pattern(subdir,files[i],beam,chan,radius,outfile,RA,DEC);
		fclose(outfile);
	}
  	return 0;
}

void leakage_pattern(char *subdir,char *mjd,int beam,int chan,float radius,FILE* outfile,float RAc,float DECc) 
{
  	char fluxfilename[80+1],header[80+1];
	FILE * fluxfile;
	int count;
	sprintf(fluxfilename,"%s/%s/beam%d/balance%04i.dat",subdir,mjd,beam,chan);
	fluxfile = fopen(fluxfilename,"r");
	if(fluxfile != NULL)
	{
	    	count = jsd_line_count(fluxfile)-1;
		printf("%s has %d records\n",fluxfilename,count);
		if(count > 0)
		{
			fgets(header,80,fluxfile);
	
		    	while (!feof(fluxfile))
			{
				float I,Q,U,V,RA,DEC,AST;
				fscanf(fluxfile,"%f %f %f %f %f %f %f",&RA, &DEC, &AST,&I, &Q, &U, &V);
				if(fabs(DEC-DECc) < radius && fabs(RA-RAc) < radius)
					fprintf(outfile,"%2.8f %2.8f %2.8f %2.8f %2.8f %2.8f\n",(RA-RAc),(DEC -DECc),I,Q,U,V);
			}
		}
		fclose(fluxfile);
	}
} 

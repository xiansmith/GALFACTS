#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "jsd/jsd_futil.h"
//#include "jsd/jsd_fit.h"

void quleakagefit(char *,char *,int,int,float,FILE*);

int main(int argc,char *argv[]) {
	//char *mjd;
	int i,beam;
        char ** files;
	char *subdir;
        int numDays;
	int chan;
	float cutoff;
	subdir = argv[1];
	chan = atoi(argv[2]);
	cutoff = atof(argv[3]);
        numDays = get_date_dirs(subdir, &files);
        if (numDays <= 0) {
                printf("ERROR: could not find any date dirs\n");
                return EXIT_FAILURE;
        }
//    	mjd = argv[1];
//    	beam = atoi(argv[2]);
    	//agg = atoi(argv[3]);

  	char outfilename[41];
	FILE *outfile;
	
	for(beam=0;beam<7;beam++)
	{
		printf("\nBEAM %d\n\n",beam);
		sprintf(outfilename,"beam%d_%04i_QU.dat",beam,chan);
		outfile = fopen(outfilename,"w");
		fprintf(outfile,"#I Q U\n");
		for(i=0;i<numDays;i++)
		    	quleakagefit(subdir,files[i],beam,chan,cutoff,outfile);
		fclose(outfile);
	}
  	return 0;
}

void quleakagefit(char *subdir,char *mjd,int beam,int chan,float cutoff,FILE* outfile) 
{
  	char fluxfilename[80+1],header[80+1];
	FILE * fluxfile;
  	//int i;
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
			//double *I,*V,C[2],chisq;
			//float nsigma = 2.5;
		    	//int k = 0;
		
			//I = (double*)malloc(sizeof(double)*count);
			//V = (double*)malloc(sizeof(double)*count);
	
			float t1;
			//double t2;
	
		    	while (!feof(fluxfile))
			{
				float I,Q,U;
				//fscanf(fluxfile,"%f %f %f %lf %lf %lf %lf",&t1, &t1, &t1,&I[k], &t2, &t2, &V[k]);
				fscanf(fluxfile,"%f %f %f %f %f %f %f",&t1, &t1, &t1,&I, &Q, &U, &t1);
				if(I > cutoff)
					fprintf(outfile,"%2.8f %2.8f %2.8f\n",I,Q,U);
	//			printf("I:%lf V:%f\n",I[k],V[k]);
		//		k++;
			}
		}
		fclose(fluxfile);
	}
//	jsd_poly_fit(I,V,count,nsigma,C,1,&chisq);
//	if(chisq > 50.0)
//	printf("Day %s Beam %d Fit c = %2.8lf m = %2.8lf chisq = %2.8lf\n",mjd,beam,C[0],C[1],chisq);
} 

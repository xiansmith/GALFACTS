#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "jsd/jsd_futil.h"
#include "jsd/jsd_fit.h"

void vleakagefit(char *,int beam);

int main(int argc,char *argv[]) {
	char *mjd;
	int beam,agg,i,b;
        char ** files;
        int numDays;
        numDays = get_date_dirs("./", &files);
        if (numDays <= 0) {
                printf("ERROR: could not find any date dirs\n");
                return EXIT_FAILURE;
        }
//    	mjd = argv[1];
//    	beam = atoi(argv[2]);
    	//agg = atoi(argv[3]);
	for(b=0;b<7;b++)
	{
		printf("\nBEAM %d\n\n",b);
		for(i=0;i<numDays;i++)
		    	vleakagefit(files[i],b);
//		    	vleakagefit(files[i],0);
	}
  	return 0;
}

void vleakagefit(char *mjd,int beam) 
{
  	char fluxfilename[40+1],header[80+1];
	FILE * fluxfile;
  	int i,count;
	sprintf(fluxfilename,"%s/beam%d/balance0001.dat",mjd,beam);
	fluxfile = fopen(fluxfilename,"r");
    	count = jsd_line_count(fluxfile)-1;
//	printf("%s has %d records\n",fluxfilename,count);
	fgets(header,80,fluxfile);
	double *I,*V,C[2],chisq;
	float nsigma = 2.5;
    	int k = 0;

	I = (double*)malloc(sizeof(double)*count);
	V = (double*)malloc(sizeof(double)*count);
	
	float t1;
	double t2;

    	while (!feof(fluxfile) && k<count)
	{
		fscanf(fluxfile,"%f %f %f %lf %lf %lf %lf",&t1, &t1, &t1,&I[k], &t2, &t2, &V[k]);
//		printf("I:%lf V:%f\n",I[k],V[k]);
		k++;
	}
	jsd_poly_fit(I,V,count,nsigma,C,1,&chisq);
//	if(chisq > 50.0)
	printf("Day %s Beam %d Fit c = %2.8lf m = %2.8lf chisq = %2.8lf\n",mjd,beam,C[0],C[1],chisq);
	fclose(fluxfile);
} 

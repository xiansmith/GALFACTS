#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

void plotPolydiffs(char *,int beam,int agg);

int main(int argc,char *argv[]) {
	char *mjd;
	int beam,agg;
    mjd = argv[1];
    beam = atoi(argv[2]);
    agg = atoi(argv[3]);
  plotPolydiffs(mjd,beam,agg);
  return 0;
}

#define ORDER 5
void plotPolydiffs(char *mjd,int beam,int agg) 
{
  FILE *gnuplotPipe,*polyFile;
  char polyFileName[40+1];
  float x[ORDER+1],x1[ORDER+1],x2[ORDER+1];
  int i;

if(agg)
{
  char ** files;
  int numDays;

   numDays = get_date_dirs("./", &files);
        if (numDays <= 0) {
                printf("ERROR: could not find any date dirs\n");
         //       return EXIT_FAILURE;
        }

int j,b;

  gnuplotPipe = popen("gnuplot","w");
      fprintf(gnuplotPipe,"set term post color\n");
      fprintf(gnuplotPipe,"set xrange [-1:1]\n");
      printf("set output \'decpolydiffsagg.ps\'\n");
      fprintf(gnuplotPipe,"set output \'decpolydiffsagg.ps\'\n");
      fprintf(gnuplotPipe,"set key off\n");
      fprintf(gnuplotPipe,"set xlabel \'Normalized DEC Range\'\n");
      fprintf(gnuplotPipe,"set title \'Stokes I Difference Polynomials\'\n");
      fprintf(gnuplotPipe,"set ylabel \'K\'\n");
    	  fprintf(gnuplotPipe,"plot ");
 for(j =0;j<numDays;j++)
 {
   for(b=0;b<7;b++)
   {	
 	 for(i =0;i<=ORDER;i++)
  	{
		x1[i] = 0.0;
		x2[i] = 0.0;
  	}

	  //sprintf(polyFileName,"polyfile1200.dat");
 	 sprintf(polyFileName,"%s/beam%d/polyfile1200.dat",files[j],b);
 	 polyFile = fopen(polyFileName,"r");
	if(polyFile == NULL)
		continue;
  	for(i =0;i<=ORDER;i++)
 	 {
		fscanf(polyFile,"%f ",&x1[i]);
//		printf("%f ",x1[i]);
  	}
 //	 printf("\n");
 	 fclose(polyFile);

	 // sprintf(polyFileName,"polyfile1201.dat");
 	 sprintf(polyFileName,"%s/beam%d/polyfile1201.dat",files[j],b);
	  polyFile = fopen(polyFileName,"r");
	if(polyFile == NULL)
		continue;
 	 for(i =0;i<=ORDER;i++)
 	 {
		fscanf(polyFile,"%f ",&x2[i]);
//		printf("%f ",x2[i]);
		x[i] = x1[i] - x2[i];
  	}
  //	printf("\n");
  	fclose(polyFile);

  	if (gnuplotPipe) {
    	  //fprintf(gnuplotPipe,"plot ");
     	 for(i = 0;i <ORDER;i++)
      	{
		fprintf(gnuplotPipe,"%f*x**%d+",x[i],i);
//		printf("%f*x**%d+",x[i],i);
      	}
    	  //fprintf(gnuplotPipe,"%f*x**%d",x[i],i);
    	  fprintf(gnuplotPipe,"%f*x**%d with lines,",x[i],i);
//      	printf("%f*x**%d\n",x[i],i);
      
      	//fprintf(gnuplotPipe," with lines\n");
 	 } else {        
      	printf("gnuplot not found...");    
  	}
  }
 }
      	fprintf(gnuplotPipe,"0*x**0 with lines\n");
      fflush(gnuplotPipe);
      fprintf(gnuplotPipe,"exit \n");    
}
else
 {
  gnuplotPipe = popen("gnuplot","w");
      fprintf(gnuplotPipe,"set term post color\n");
      fprintf(gnuplotPipe,"set xrange [-1:1]\n");
      printf("set output \'decpolydiffs%sb%d.ps\'\n",mjd,beam);
      fprintf(gnuplotPipe,"set output \'decpolydiffs%sb%d.ps\'\n",mjd,beam);
      fprintf(gnuplotPipe,"set key off\n");
      fprintf(gnuplotPipe,"set xlabel \'Normalized DEC Range\'\n");
      fprintf(gnuplotPipe,"set title \'Stokes I Difference Polynomial\'\n");
      fprintf(gnuplotPipe,"set ylabel \'K\'\n");
    	  fprintf(gnuplotPipe,"plot ");
 	
	 for(i =0;i<=ORDER;i++)
  	{
		x1[i] = 0.0;
		x2[i] = 0.0;
  	}

	  sprintf(polyFileName,"polyfile1200.dat");
 	 polyFile = fopen(polyFileName,"r");
//	if(polyFile == NULL)
//		continue;
  	for(i =0;i<=ORDER;i++)
 	 {
		fscanf(polyFile,"%f ",&x1[i]);
//		printf("%f ",x1[i]);
  	}
 //	 printf("\n");
 	 fclose(polyFile);

	  sprintf(polyFileName,"polyfile1201.dat");
	  polyFile = fopen(polyFileName,"r");
//	if(polyFile == NULL)
//		continue;
 	 for(i =0;i<=ORDER;i++)
 	 {
		fscanf(polyFile,"%f ",&x2[i]);
//		printf("%f ",x2[i]);
		x[i] = x1[i] - x2[i];
  	}
  //	printf("\n");
  	fclose(polyFile);

  	if (gnuplotPipe) {
    	  //fprintf(gnuplotPipe,"plot ");
     	 for(i = 0;i <ORDER;i++)
      	{
		fprintf(gnuplotPipe,"%f*x**%d+",x[i],i);
//		printf("%f*x**%d+",x[i],i);
      	}
    	  //fprintf(gnuplotPipe,"%f*x**%d",x[i],i);
    	  fprintf(gnuplotPipe,"%f*x**%d with lines,",x[i],i);
      
 	 } else {        
      	printf("gnuplot not found...");    
  	}
      	fprintf(gnuplotPipe,"0*x**0 with lines\n");
      fflush(gnuplotPipe);
      fprintf(gnuplotPipe,"exit \n");    
 }
} 

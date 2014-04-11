#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "jsd/jsd_fit.h"
#include "leakagefit.h"
#include "fluxdata.h"

void vleakagefit(FluxWappData * wappdata,int chan) 
{
  	int i,count;
	double *I,*V,C[2],chisq;
	float nsigma = 2.5;
    	int k = 0;


        for (i=0; i<wappdata->numDays; i++) 
        {
                FluxDayData * daydata;

                daydata = &wappdata->daydata[i];
		count = daydata->numRecords;
		I = (double*)malloc(sizeof(double)*count);
		V = (double*)malloc(sizeof(double)*count);
                for (k=0; k<count; k++) 
                {			
                        I[k] = (double)daydata->records[k].stokes.I;
                        V[k] = (double)daydata->records[k].stokes.V;
                }
		jsd_poly_fit(I,V,count,nsigma,C,1,&chisq);
	//	if(chisq > 50.0)
	//	printf("Day %s Beam %d Fit c = %2.8lf m = %2.8lf chisq = %2.8lf\n",daydata->mjd,i%7,C[0],C[1],chisq);
                for (k=0; k<count; k++) 
                {			
                        daydata->records[k].stokes.V -= (daydata->records[k].stokes.I*C[1] + C[0]);
                        //V[k] = (double)daydata->records[k].stokes.V;
                }
        }
	
} 

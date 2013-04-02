#include <stdlib.h>
#include <math.h>
#include "calibrate.h"
#include "rfi.h"
#include "chebyshev.h"
#include <string.h>
//----------------------------------------------------------------------------------------------------------
//pre-compute the raw cal for every channel
//void compute_raw_cal(SpecRecord dataset[], int size)

void compute_raw_cal(SpecRecord dataset[], int size, int lowchan, int highchan)
{
	int n, i;

//To make plots of the fits
        FILE * outfile = fopen("rawbandcal.dat","w");
        if(outfile == NULL)
        {
                printf("Cant open rawbandcal.dat\n");
                exit(1);
        }
        PolAvg avgcal; 

	for(n=0; n<size; n++)
	{
                int count =0;
                memset(&avgcal, 0, sizeof(PolAvg));
                if(dataset[n].flagBAD) continue;

                for(i=lowchan; i<highchan; i++)
                {

                        if(dataset[n].flagRFI[i] == RFI_NONE)
                        {
                               dataset[n].cal.xx[i] = dataset[n].calon.xx[i] - dataset[n].caloff.xx[i];
                               dataset[n].cal.yy[i] = dataset[n].calon.yy[i] - dataset[n].caloff.yy[i];
                               dataset[n].cal.xy[i] = dataset[n].calon.xy[i] - dataset[n].caloff.xy[i];
                               dataset[n].cal.yx[i] = dataset[n].calon.yx[i] - dataset[n].caloff.yx[i];
                               avgcal.xx += dataset[n].cal.xx[i];
                               avgcal.yy += dataset[n].cal.yy[i];
                               avgcal.xy += dataset[n].cal.xy[i];
                               avgcal.yx += dataset[n].cal.yx[i];
                               count++;
                        }
                 }
                 fprintf(outfile, "%05i %7.6f %7.6f %7.6f %7.6f\n",n,\
                 avgcal.xx/count, avgcal.yy/count, avgcal.xy/count, avgcal.yx/count);
                 printf("%f %%\r", (n + 1)*100.0/size);
	}

	printf("\n");
        fclose(outfile);
}
//----------------------------------------------------------------------------------------------------------
void linear_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int RFIF)
{
	int n, chan, order = 1;
	float C[order + 1]; 
	float nsigma = 3.0;
	float sumsq;
	float *XRA, *Yxx, *Yyy, *Yxy, *Yyx;
	float min, max;
	size_t count;

	XRA = (float*) malloc(sizeof(float) * size);
	Yxx = (float*) malloc(sizeof(float) * size);
	Yyy = (float*) malloc(sizeof(float) * size);
	Yxy = (float*) malloc(sizeof(float) * size);
	Yyx = (float*) malloc(sizeof(float) * size);

	for(chan=lowchan; chan<highchan; chan++) 
		{
		count = 0;
		for(n=0; n<size; n++) 
			{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
				{
				XRA[count] = dataset[n].RA;
				Yxx[count] = dataset[n].cal.xx[chan];
				Yyy[count] = dataset[n].cal.yy[chan];
				Yxy[count] = dataset[n].cal.xy[chan];
				Yyx[count] = dataset[n].cal.yx[chan];
				count++;
				}
			}

		// normalize the X values for the curve fit 	
		chebyshev_minmax(XRA, count, &min, &max);
		chebyshev_normalize(XRA, count, min, max);

		//XX
		chebyshev_fit_bw(XRA, Yxx, count, nsigma, C, order);
		for(n=0; n<size; n++) 
			{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
				{
				dataset[n].cal.xx[chan] = chebyshev_eval(CNORMALIZE(dataset[n].RA,min,max), C, order);
				}
			}

		//YY
		chebyshev_fit_bw(XRA, Yyy, count, nsigma, C, order);		
		for(n=0; n<size; n++) 
			{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
				{
				dataset[n].cal.yy[chan] = chebyshev_eval(CNORMALIZE(dataset[n].RA,min,max), C, order);
				}
			}
		
		//XY
		chebyshev_fit_bw(XRA, Yxy, count, nsigma, C, order);
		for(n=0; n<size; n++) 
			{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
				{
				dataset[n].cal.xy[chan] = chebyshev_eval(CNORMALIZE(dataset[n].RA,min,max), C, order);
				}
			}

		//YX
		chebyshev_fit_bw(XRA, Yyx, count, nsigma, C, order);		
		for(n=0; n<size; n++) 
			{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
				{
				dataset[n].cal.yx[chan] = chebyshev_eval(CNORMALIZE(dataset[n].RA,min,max), C, order);
				}
			}
		
		//printf("%f %%\r", (chan - lowchan + 1)*100.0/(highchan - lowchan));
		}
	printf("\n");
	free(XRA);
	free(Yxx);
	free(Yyy);
	free(Yxy);
	free(Yyx);
}
//----------------------------------------------------------------------------------------------------------
void smooth_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int window)
{
	int n, chan, dchan = highchan - lowchan;
	float mean[4], sigma[4], tmp;
	float *XRA, *Yxx, *Yyy, *Yxy, *Yyx, *Fxx, *Fyy, *Fxy, *Fyx, *w;
//	FILE *fxx, *fyy, *fxy, *fyx;

	Yxx = (float*) malloc(sizeof(float) * size);
	Yyy = (float*) malloc(sizeof(float) * size);
	Yxy = (float*) malloc(sizeof(float) * size);
	Yyx = (float*) malloc(sizeof(float) * size);
	
	Fxx = (float*) malloc(sizeof(float) * dchan);
	Fyy = (float*) malloc(sizeof(float) * dchan);
	Fxy = (float*) malloc(sizeof(float) * dchan);
	Fyx = (float*) malloc(sizeof(float) * dchan);

	w = (float*) malloc(sizeof(float) * (window+1));
	
	for(n=0; n<=window; n++) {tmp = n*2.0/window; w[n] = exp(-0.5*tmp*tmp);}
	tmp = w[0]; for(n=1; n<=window; n++) tmp += 2*w[n];
	for(n=0; n<=window; n++) w[n] /= tmp;

	for(chan=lowchan; chan<highchan; chan++) 
		{
		mean[0] = mean[1] = mean[2] = mean[3] = 0;
		int count = 0;
		for(n=0; n<size; n++) 
			{
			mean[0] += dataset[n].cal.xx[chan];
			mean[1] += dataset[n].cal.yy[chan];
			mean[2] += dataset[n].cal.xy[chan];
			mean[3] += dataset[n].cal.yx[chan];
			count++;
			}
		//avoid divide by zero
		if(count)
		{
			mean[0] /= size; 
			mean[1] /= size; 
			mean[2] /= size; 
			mean[3] /= size;
		}
		sigma[0] = sigma[1] = sigma[2] = sigma[3] = 0; 
		count = 0;
		for(n=0; n<size; n++) 
			{
			sigma[0] += (dataset[n].cal.xx[chan] - mean[0])*(dataset[n].cal.xx[chan] - mean[0]);
			sigma[1] += (dataset[n].cal.yy[chan] - mean[1])*(dataset[n].cal.yy[chan] - mean[1]);
			sigma[2] += (dataset[n].cal.xy[chan] - mean[2])*(dataset[n].cal.xy[chan] - mean[2]);
			sigma[3] += (dataset[n].cal.yx[chan] - mean[3])*(dataset[n].cal.yx[chan] - mean[3]);
			count++;
			}
		if (count)
		{
			sigma[0] = sqrt(sigma[0]/size); 
			sigma[1] = sqrt(sigma[1]/size); 
			sigma[2] = sqrt(sigma[2]/size); 
			sigma[3] = sqrt(sigma[3]/size);
		}
		
		Fxx[chan - lowchan] = Fyy[chan - lowchan] = Fxy[chan - lowchan] = Fyx[chan - lowchan] = 0;
		for(n=0; n<size; n++) 
			{
			if(fabs(dataset[n].cal.xx[chan] - mean[0]) > 3.0*sigma[0]) dataset[n].cal.xx[chan] = mean[0];
			if(fabs(dataset[n].cal.yy[chan] - mean[1]) > 3.0*sigma[1]) dataset[n].cal.yy[chan] = mean[1];
			if(fabs(dataset[n].cal.xy[chan] - mean[2]) > 3.0*sigma[2]) dataset[n].cal.xy[chan] = mean[2];
			if(fabs(dataset[n].cal.yx[chan] - mean[3]) > 3.0*sigma[3]) dataset[n].cal.yx[chan] = mean[3];
			Fxx[chan - lowchan] += dataset[n].cal.xx[chan];
			Fyy[chan - lowchan] += dataset[n].cal.yy[chan];
			Fxy[chan - lowchan] += dataset[n].cal.xy[chan];
			Fyx[chan - lowchan] += dataset[n].cal.yx[chan];
			}
		Fxx[chan - lowchan] /= size;
		Fyy[chan - lowchan] /= size;
		Fxy[chan - lowchan] /= size;
		Fyx[chan - lowchan] /= size;
		}

	mean[0] = mean[1] = mean[2] = mean[3] = 0;
	for(n=0; n<size; n++) 
		{
		Yxx[n] = Yyy[n] = Yxy[n] = Yyx[n] = 0;
		for(chan=lowchan; chan<highchan; chan++) 
			{
			Yxx[n] += dataset[n].cal.xx[chan];
			Yyy[n] += dataset[n].cal.yy[chan];
			Yxy[n] += dataset[n].cal.xy[chan];
			Yyx[n] += dataset[n].cal.yx[chan];
			}
		Yxx[n] /= dchan; 
		Yyy[n] /= dchan; 
		Yxy[n] /= dchan; 

		Yyx[n] /= dchan;
		mean[0] += Yxx[n];
		mean[1] += Yyy[n];
		mean[2] += Yxy[n];
		mean[3] += Yyx[n];
		}
	mean[0] /= size; 
	mean[1] /= size; 
	mean[2] /= size; 
	mean[3] /= size;	
	
	for(n=0; n<size; n++) 
		{
		Yxx[n] /= mean[0];
		Yyy[n] /= mean[1];
		Yxy[n] /= mean[2];
		Yyx[n] /= mean[3];
		}	

/*
	fxx=fopen("acxx.dat", "w");
	fyy=fopen("acyy.dat", "w");
	fxy=fopen("acxy.dat", "w");
	fyx=fopen("acyx.dat", "w");

	for(n=0; n<size; n++) 
		{
		fprintf(fxx, "%f %f\n", dataset[n].RA, Yxx[n]);
		fprintf(fyy, "%f %f\n", dataset[n].RA, Yyy[n]);
		fprintf(fxy, "%f %f\n", dataset[n].RA, Yxy[n]);
		fprintf(fyx, "%f %f\n", dataset[n].RA, Yyx[n]);
		}
	fclose(fxx);
	fclose(fyy);
	fclose(fxy);
	fclose(fyx);
*/
/*	
	moving_average_filter(Yxx, size, window);
	moving_average_filter(Yyy, size, window);
	moving_average_filter(Yxy, size, window);
	moving_average_filter(Yyx, size, window);
*/

/*	gaussian_filter(Yxx, size, w, window);
	gaussian_filter(Yyy, size, w, window);
	gaussian_filter(Yxy, size, w, window);
	gaussian_filter(Yyx, size, w, window);
*/

	//In this case window is actually the number of iterations for diffusion smoothing
	//Carried over the variable name from the other fitting routines that actually 
	//used a smoothing window.
        diffusion_filter(Yxx, size, window);
        diffusion_filter(Yyy, size, window);
        diffusion_filter(Yxy, size, window);
        diffusion_filter(Yyx, size, window);

	//apply the moving average filter to reduce noise
	moving_average_filter(Yxx, size, window);
	moving_average_filter(Yyy, size, window);
	moving_average_filter(Yxy, size, window);
	moving_average_filter(Yyx, size, window);
/*	
	fxx=fopen("scxx.dat", "w");
	fyy=fopen("scyy.dat", "w");
	fxy=fopen("scxy.dat", "w");
	fyx=fopen("scyx.dat", "w");
	for(n=0; n<size; n++) 
		{
		fprintf(fxx, "%f %f\n", dataset[n].RA, Yxx[n]);
		fprintf(fyy, "%f %f\n", dataset[n].RA, Yyy[n]);
		fprintf(fxy, "%f %f\n", dataset[n].RA, Yxy[n]);
		fprintf(fyx, "%f %f\n", dataset[n].RA, Yyx[n]);
		}
	fclose(fxx);
	fclose(fyy);
	fclose(fxy);
	fclose(fyx);
*/	
				
	for(chan=lowchan; chan<highchan; chan++) 
		{
		for(n=0; n<size; n++) 
			{
			dataset[n].cal.xx[chan] = Yxx[n]*Fxx[chan-lowchan];
			dataset[n].cal.yy[chan] = Yyy[n]*Fyy[chan-lowchan];
			dataset[n].cal.xy[chan] = Yxy[n]*Fxy[chan-lowchan];
			dataset[n].cal.yx[chan] = Yyx[n]*Fyx[chan-lowchan];
			}
		//printf("%f %%\r", (chan - lowchan + 1)*100.0/(highchan - lowchan));
		}
	printf("\n");

//For making plots of the fits
        FILE * outfile2 = fopen("smoothcal.dat","w");
        if(outfile2 == NULL)
        {
                printf("Can't open smoothcal.dat\n");
                exit(1);
        }
        PolAvg avgcal;
        int i;
        for(n=0; n<size; n++)
        {
                int count = 1;
                memset(&avgcal, 0, sizeof(PolAvg));
                if(dataset[n].flagBAD) continue;
                for(i=lowchan; i<highchan; i++)
                {
                       if(dataset[n].flagRFI[i] == RFI_NONE)
                       {
                                avgcal.xx += dataset[n].cal.xx[i];
                                avgcal.yy += dataset[n].cal.yy[i];
                                avgcal.xy += dataset[n].cal.xy[i];
                                avgcal.yx += dataset[n].cal.yx[i];
                                count++;
                       }
                }
                fprintf(outfile2, "%05i %7.6f %7.6f %7.6f %7.6f\n",n,\
                avgcal.xx/count, avgcal.yy/count, avgcal.xy/count, avgcal.yx/count);
                printf("%f %%\r", (n + 1)*100.0/size);
        }

        fclose(outfile2);

	free(Yxx);
	free(Yyy);
	free(Yxy);
	free(Yyx);
	free(Fxx);
	free(Fyy);
	free(Fxy);
	free(Fyx);
	free(w);
}
//----------------------------------------------------------------------------------------------------------

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
                //if(dataset[n].flagBAD) continue;

                for(i=0; i<MAX_CHANNELS; i++)
                {

                        //if(dataset[n].flagRFI[i] == RFI_NONE)
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
		
                 //printf("%f %%\r", (n + 1)*100.0/size);
	}

	printf("\n");
        fclose(outfile);
	return;
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
	//size_t count;
	int count;

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
	return;
}
//----------------------------------------------------------------------------------------------------------
//void smooth_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int window)
void smooth_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int window, int cwindow)
{
	//fix this loop need to run for all 4096 channels
	lowchan = 0; highchan = MAX_CHANNELS;

	int n, chan, dchan = highchan - lowchan;
	float mean[4], sigma[4], tmp;
	float *XRA, *Yxx, *Yyy, *Yxy, *Yyx, *Fxx, *Fyy, *Fxy, *Fyx, *w;
	float *Fxxs, *Fyys, *Fxys, *Fyxs;

	Yxx = (float*) malloc(sizeof(float) * size);
	Yyy = (float*) malloc(sizeof(float) * size);
	Yxy = (float*) malloc(sizeof(float) * size);
	Yyx = (float*) malloc(sizeof(float) * size);
	
	float **Sxx,**Syy,**Syx,**Sxy;

	Sxx = (float**) malloc(sizeof(float *) * size);
	Syy = (float**) malloc(sizeof(float *) * size);
	Sxy = (float**) malloc(sizeof(float *) * size);
	Syx = (float**) malloc(sizeof(float *) * size);
	for(n=0; n<size; n++) 
	{
		Sxx[n] = (float*) malloc(sizeof(float) * MAX_CHANNELS);
		Syy[n] = (float*) malloc(sizeof(float) * MAX_CHANNELS);
		Sxy[n] = (float*) malloc(sizeof(float) * MAX_CHANNELS);
		Syx[n] = (float*) malloc(sizeof(float) * MAX_CHANNELS);
	}

	Fxx = (float*) malloc(sizeof(float) * dchan);
	Fyy = (float*) malloc(sizeof(float) * dchan);
	Fxy = (float*) malloc(sizeof(float) * dchan);
	Fyx = (float*) malloc(sizeof(float) * dchan);

	Fxxs = (float*) malloc(sizeof(float) * dchan);
	Fyys = (float*) malloc(sizeof(float) * dchan);
	Fxys = (float*) malloc(sizeof(float) * dchan);
	Fyxs = (float*) malloc(sizeof(float) * dchan);

	if(Yxx == NULL || Yyy == NULL || Yxy == NULL || Yyx == NULL)
	{
		printf("malloc failed !\n");
		exit(0);
	}
	if(Fxx == NULL || Fyy == NULL || Fxy == NULL || Fyx == NULL)
	{
		printf("malloc failed !\n");
		exit(0);
	}
	if(Fxxs == NULL || Fyys == NULL || Fxys == NULL || Fyxs == NULL)
	{
		printf("malloc failed !\n");
		exit(0);
	}


        FILE * outfile2 = fopen("smoothcal1.dat","w");
        if(outfile2 == NULL)
        {
                printf("Can't open smoothcal1.dat\n");
                exit(0);
        }

	for(chan=lowchan; chan<highchan; chan++) 
	{
		mean[0] = mean[1] = mean[2] = mean[3] = 0;
		int count[MAX_CHANNELS] = {0};
		for(n=0; n<size; n++) 
		{
			//if(dataset[n].flagBAD) continue;
                        //if(dataset[n].flagRFI[chan] == RFI_NONE)
			{
				mean[0] += dataset[n].cal.xx[chan];
				mean[1] += dataset[n].cal.yy[chan];
				mean[2] += dataset[n].cal.xy[chan];
				mean[3] += dataset[n].cal.yx[chan];
				count[chan]++;
			}
		}

		//avoid divide by zero
		if(count[chan])
		{
			mean[0] /= count[chan]; 
			mean[1] /= count[chan]; 
			mean[2] /= count[chan]; 
			mean[3] /= count[chan];
		}
	
		fprintf(outfile2,"%f %f %f %f\n",mean[0],mean[1],mean[2],mean[3]);		

		/*sigma[0] = sigma[1] = sigma[2] = sigma[3] = 0; 
		count[chan] = 0;
		for(n=0; n<size; n++) 
		{
			if(dataset[n].flagBAD) continue;
                        if(dataset[n].flagRFI[chan] == RFI_NONE)
			{
			sigma[0] += (dataset[n].cal.xx[chan] - mean[0])*(dataset[n].cal.xx[chan] - mean[0]);
			sigma[1] += (dataset[n].cal.yy[chan] - mean[1])*(dataset[n].cal.yy[chan] - mean[1]);
			sigma[2] += (dataset[n].cal.xy[chan] - mean[2])*(dataset[n].cal.xy[chan] - mean[2]);
			sigma[3] += (dataset[n].cal.yx[chan] - mean[3])*(dataset[n].cal.yx[chan] - mean[3]);
			count[chan]++;
			}
		}
		if(count[chan])
		{
			sigma[0] = sqrt(sigma[0]/count[chan]); 
			sigma[1] = sqrt(sigma[1]/count[chan]); 
			sigma[2] = sqrt(sigma[2]/count[chan]); 
			sigma[3] = sqrt(sigma[3]/count[chan]);
		}*/
		
		/*for(n=0; n<size; n++) 
		{
			//if(dataset[n].flagBAD || dataset[n].flagRFI[chan] != RFI_NONE)
			//if(dataset[n].flagBAD)
			{
				dataset[n].cal.xx[chan] = mean[0];
				dataset[n].cal.yy[chan] = mean[1];
				dataset[n].cal.xy[chan] = mean[2];
				dataset[n].cal.yx[chan] = mean[3];
			}
		}*/
	}
	fclose(outfile2);

	mean[0] = mean[1] = mean[2] = mean[3] = 0;
	sigma[0] = sigma[1] = sigma[2] = sigma[3] = 0; 
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
		sigma[0] += (Yxx[n] - mean[0])*(Yxx[n] - mean[0]);
		sigma[1] += (Yyy[n] - mean[1])*(Yyy[n] - mean[1]);
		sigma[2] += (Yxy[n] - mean[2])*(Yxy[n] - mean[2]);
		sigma[3] += (Yyx[n] - mean[3])*(Yyx[n] - mean[3]);
	}

	sigma[0] = sqrt(sigma[0]/size); 
	sigma[1] = sqrt(sigma[1]/size); 
	sigma[2] = sqrt(sigma[2]/size); 
	sigma[3] = sqrt(sigma[3]/size); 

	for(n=0; n<size; n++) 
	{
		if(fabs(Yxx[n] - mean[0]) > 3.0*sigma[0]) Yxx[n] = mean[0];
		if(fabs(Yyy[n] - mean[1]) > 3.0*sigma[1]) Yyy[n] = mean[1];
		if(fabs(Yxy[n] - mean[2]) > 3.0*sigma[2]) Yxy[n] = mean[2];
		if(fabs(Yyx[n] - mean[3]) > 3.0*sigma[3]) Yyx[n] = mean[3];
	}
	
	for(n=0; n<size; n++) 
	{
		Yxx[n] /= mean[0];
		Yyy[n] /= mean[1];
		Yxy[n] /= mean[2];
		Yyx[n] /= mean[3];
	}

	//In this case window is actually the number of iterations for diffusion smoothing
	//Carried over the variable name from the other fitting routines that actually 
	//used a smoothing window.
        diffusion_filter(Yxx, size, window);
        diffusion_filter(Yyy, size, window);
        diffusion_filter(Yxy, size, window);
        diffusion_filter(Yyx, size, window);

/*        FILE * outfile2 = fopen("smoothcal1.dat","w");
        if(outfile2 == NULL)
        {
                printf("Can't open smoothcal1.dat\n");
                exit(0);
        }

        for(n=0; n<size; n++)
        {
                if(dataset[n].flagBAD) continue;
                //fprintf(outfile2, "%05i %7.6f %7.6f %7.6f %7.6f\n",n,Yxx[n]*Fxxm,Yyy[n]*Fyym,Yxy[n]*Fxym,Yyx[n]*Fyxm);
                fprintf(outfile2, "%05i %7.6f %7.6f %7.6f %7.6f\n",n,Yxx[n],Yyy[n],Yxy[n],Yyx[n]);
                //printf("%f %%\r", (n + 1)*100.0/size);
	}
	fclose(outfile2);
*/
	//apply the moving average filter to reduce noise
	moving_average_filter(Yxx, size, 500);
	moving_average_filter(Yyy, size, 500);
	moving_average_filter(Yxy, size, 500);
	moving_average_filter(Yyx, size, 500);

	for(chan=lowchan; chan<highchan; chan++)
 	{
		Fxx[chan - lowchan] = 0.0;
		Fyy[chan - lowchan] = 0.0;
		Fxy[chan - lowchan] = 0.0;
		Fyx[chan - lowchan] = 0.0;
	}

	int cc[MAX_CHANNELS] = {0};
	for(n=0; n<size; n++) 
	{
		for(chan=lowchan; chan<highchan; chan++) 
		{

			if(n == 0)
			{
				cc[chan] = 0;
				int nn;
				for(nn=0; nn<cwindow; nn++) 
				{
					//if(!dataset[nn].flagBAD && dataset[nn].flagRFI[chan] == RFI_NONE)
					{

						Fxx[chan - lowchan] += dataset[nn].cal.xx[chan];
						Fyy[chan - lowchan] += dataset[nn].cal.yy[chan];
						Fxy[chan - lowchan] += dataset[nn].cal.xy[chan];
						Fyx[chan - lowchan] += dataset[nn].cal.yx[chan];
						cc[chan]++;
					}
				}
				if(cc[chan])
				{
					Fxx[chan - lowchan] /= cc[chan];
					Fyy[chan - lowchan] /= cc[chan];
					Fxy[chan - lowchan] /= cc[chan];
					Fyx[chan - lowchan] /= cc[chan];
				}
				else
				{
					printf("Rec %d chan %d count zero.\n",n,chan);
					Fxx[chan - lowchan] = Fxx[chan-lowchan-1];
					Fyy[chan - lowchan] = Fyy[chan-lowchan-1];
					Fxy[chan - lowchan] = Fxy[chan-lowchan-1];
					Fyx[chan - lowchan] = Fyx[chan-lowchan-1];
				}
			}
			else if(n < cwindow/2 && n > 0)
			{
				;
			}
			else if(n >= cwindow/2 && n < size - cwindow/2)
			{
				int a = n - cwindow/2;
				int b = n + cwindow/2;
				Fxx[chan - lowchan] = Fxx[chan - lowchan]*cc[chan] - dataset[a].cal.xx[chan] + dataset[b].cal.xx[chan];
				Fyy[chan - lowchan] = Fyy[chan - lowchan]*cc[chan] - dataset[a].cal.yy[chan] + dataset[b].cal.yy[chan];
				Fxy[chan - lowchan] = Fxy[chan - lowchan]*cc[chan] - dataset[a].cal.xy[chan] + dataset[b].cal.xy[chan];
				Fyx[chan - lowchan] = Fyx[chan - lowchan]*cc[chan] - dataset[a].cal.yx[chan] + dataset[b].cal.yx[chan];
				Fxx[chan - lowchan] = Fxx[chan - lowchan]/cc[chan];
				Fyy[chan - lowchan] = Fyy[chan - lowchan]/cc[chan];
				Fxy[chan - lowchan] = Fxy[chan - lowchan]/cc[chan];
				Fyx[chan - lowchan] = Fyx[chan - lowchan]/cc[chan];
				/*
				//if(!dataset[a].flagBAD && dataset[a].flagRFI[chan] == RFI_NONE)
				{
					if(cc[chan])
					{
						Fxx[chan - lowchan] = Fxx[chan - lowchan]*cc[chan] - (dataset[a].cal.xx[chan]);
						Fyy[chan - lowchan] = Fyy[chan - lowchan]*cc[chan] - (dataset[a].cal.yy[chan]);
						Fxy[chan - lowchan] = Fxy[chan - lowchan]*cc[chan] - (dataset[a].cal.xy[chan]);
						Fyx[chan - lowchan] = Fyx[chan - lowchan]*cc[chan] - (dataset[a].cal.yx[chan]);
						cc[chan]--;
						if(cc[chan])
						{
							Fxx[chan - lowchan] = Fxx[chan - lowchan]/cc[chan];
							Fyy[chan - lowchan] = Fyy[chan - lowchan]/cc[chan];
							Fxy[chan - lowchan] = Fxy[chan - lowchan]/cc[chan];
							Fyx[chan - lowchan] = Fyx[chan - lowchan]/cc[chan];
						}
						else
						{
							Fxx[chan - lowchan] = Fxx[chan-lowchan-1];
							Fyy[chan - lowchan] = Fyy[chan-lowchan-1];
							Fxy[chan - lowchan] = Fxy[chan-lowchan-1];
							Fyx[chan - lowchan] = Fyx[chan-lowchan-1];
						}
					}
				}
				//if(!dataset[b].flagBAD && dataset[b].flagRFI[chan] == RFI_NONE)
				{
					if(cc[chan])
					{
						Fxx[chan-lowchan] = Fxx[chan-lowchan]*cc[chan];
						Fyy[chan-lowchan] = Fyy[chan-lowchan]*cc[chan];
						Fxy[chan-lowchan] = Fxy[chan-lowchan]*cc[chan];
						Fyx[chan-lowchan] = Fyx[chan-lowchan]*cc[chan];
					}	
					cc[chan]++;
					Fxx[chan - lowchan] = (Fxx[chan - lowchan] + dataset[b].cal.xx[chan])/cc[chan];
					Fyy[chan - lowchan] = (Fyy[chan - lowchan] + dataset[b].cal.yy[chan])/cc[chan];
					Fxy[chan - lowchan] = (Fxy[chan - lowchan] + dataset[b].cal.xy[chan])/cc[chan];
					Fyx[chan - lowchan] = (Fyx[chan - lowchan] + dataset[b].cal.yx[chan])/cc[chan];
				}*/
			}
			else
			{
				;
			}

			if(Fxx[chan-lowchan] == 0.0)
			{
				Fxxs[chan-lowchan] = Fxx[chan - lowchan-1];
				Fyys[chan-lowchan] = Fyy[chan - lowchan-1];
				Fxys[chan-lowchan] = Fxy[chan - lowchan-1];
				Fyxs[chan-lowchan] = Fyx[chan - lowchan-1];
			}
			else
			{
				Fxxs[chan-lowchan] = Fxx[chan - lowchan];
				Fyys[chan-lowchan] = Fyy[chan - lowchan];
				Fxys[chan-lowchan] = Fxy[chan - lowchan];
				Fyxs[chan-lowchan] = Fyx[chan - lowchan];
			}
		}		

		//moving_average_filter(Fxxs, dchan, 100);
		//moving_average_filter(Fyys, dchan, 100);
		//moving_average_filter(Fxys, dchan, 100);
		//moving_average_filter(Fyxs, dchan, 100);

		for(chan=lowchan; chan<highchan; chan++) 
		{
			Sxx[n][chan] = Yxx[n]*Fxxs[chan-lowchan];
			Syy[n][chan] = Yyy[n]*Fyys[chan-lowchan];
			Sxy[n][chan] = Yxy[n]*Fxys[chan-lowchan];
			Syx[n][chan] = Yyx[n]*Fyxs[chan-lowchan];

		}
                if(!(n%1000))
                {
                        write_band(Fxxs,Fyys, Fxys, Fyxs,n,lowchan,highchan);
                }


	}
	printf("\n");

        outfile2 = fopen("smoothcal2.dat","w");
        if(outfile2 == NULL)
        {
                printf("Can't open smoothcal2.dat\n");
                exit(0);
        }
        for(n=0; n<size; n++)
        {
		if(dataset[n].flagBAD)
			continue;
                fprintf(outfile2, "%05i %7.6f %7.6f %7.6f %7.6f\n",n,Yxx[n],Yyy[n],Yxy[n],Yyx[n]);
		for(chan=lowchan; chan<highchan; chan++) 
		{
			dataset[n].cal.xx[chan] = Sxx[n][chan];
			dataset[n].cal.yy[chan] = Syy[n][chan];
			dataset[n].cal.xy[chan] = Sxy[n][chan];
			dataset[n].cal.yx[chan] = Syx[n][chan];
		}
	}
	fclose(outfile2);

/*        PolAvg avgcal;
        int i;
	outfile2 = fopen("smoothcal.dat","w");
        for(n=0; n<size; n++)
        {
                int count = 0;
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
                //if(dataset[n].flagBAD) continue;
                //fprintf(outfile2, "%05i %7.6f %7.6f %7.6f %7.6f\n",n,Yxx[n],Yyy[n],Yxy[n],Yyx[n]);
                printf("%f %%\r", (n + 1)*100.0/size);
        }

        fclose(outfile2);
*/
	free(Yxx);
	free(Yyy);
	free(Yxy);
	free(Yyx);
	free(Fxx);
	free(Fyy);
	free(Fxy);
	free(Fyx);
	free(Fxxs);
	free(Fyys);
	free(Fxys);
	free(Fyxs);
	for(n=0; n<size; n++) 
	{
		free(Sxx[n]);
		free(Syy[n]);
		free(Sxy[n]);
		free(Syx[n]);
	}
	free(Sxx);
	free(Syy);
	free(Sxy);
	free(Syx);
	return;
}
//----------------------------------------------------------------------------------------------------------
void write_band(float* Fxx,float* Fyy, float *Fxy, float *Fyx,int r, int low, int high)
{
        char filename[64];
        int i;
        sprintf(filename,"bandw_%04i.dat",r);
        FILE * tcalfile = fopen(filename,"w");
        if(tcalfile == NULL)
        {
                printf("Can't open band.dat\n");
                exit(1);
        }
        for (i=low;i<high;i++)
        {
                        fprintf(tcalfile,"%2.6f %2.6f %2.6f %2.6f\n",Fxx[i-low],Fyy[i-low],Fxy[i-low],Fyx[i-low]);
                        //printf("%2.6f %2.6f\n",tcalxx[i],tcalyy[i]);
                        if(Fxx[i-low] > 1000.0 || Fyy[i-low] > 1000.0 || Fxy[i-low] > 1000.0 || Fyx[i-low] > 1000.0)
                        {
                                printf("High Value\n");
                        }
        }

        fclose(tcalfile);
	return;
}


#include <stdlib.h>
#include <math.h>
#include "calibrate.h"
#include "rfi.h"
#include "chebyshev.h"
#include "denoising.h"
#include <string.h>
//----------------------------------------------------------------------------------------------------------
//pre-compute the raw cal for every channel
//void compute_raw_cal(SpecRecord dataset[], int size)

void compute_raw_cal(SpecRecord dataset[], int size, int lowchan, int highchan)
{
	int n, i;

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

                for(i=0; i<MAX_CHANNELS; i++)
                //for(i=lowchan; i<highchan; i++)
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
void linear_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int RFIF,int fwindow)
{
	int n, chan, order = 1;
	float C[order + 1]; 
	float nsigma = 2.0;
	float *XRA, *Yxx, *Yyy, *Yxy, *Yyx;
	float *Cxx, *Cyy, *Cxy, *Cyx;
	float min, max;
	int count;

	XRA = (float*)calloc(size,sizeof(float));
	Yxx = (float*)calloc(size,sizeof(float));
	Yyy = (float*)calloc(size,sizeof(float));
	Yxy = (float*)calloc(size,sizeof(float));
	Yyx = (float*)calloc(size,sizeof(float));

	Cxx = (float*)calloc((highchan-lowchan),sizeof(float));
	Cyy = (float*)calloc((highchan-lowchan),sizeof(float));
	Cxy = (float*)calloc((highchan-lowchan),sizeof(float));
	Cyx = (float*)calloc((highchan-lowchan),sizeof(float));

	for(chan=lowchan; chan<highchan; chan++) 
	{
		count = 0;
		
		for(n=0; n<size; n++) 
		{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
			{
				XRA[count] = dataset[n].RA;
				Yxx[count] += dataset[n].cal.xx[chan];
				Yyy[count] += dataset[n].cal.yy[chan];
				count++;
				Cxx[chan-lowchan] += dataset[n].cal.xx[chan];
				Cyy[chan-lowchan] += dataset[n].cal.yy[chan];
				Cxy[chan-lowchan] += dataset[n].cal.xy[chan];
				Cyx[chan-lowchan] += dataset[n].cal.yx[chan];
			}
		}

		if(count)
		{
			Cxx[chan-lowchan]/=count;
			Cyy[chan-lowchan]/=count;
			Cxy[chan-lowchan]/=count;
			Cyx[chan-lowchan]/=count;
		}
	}

	//printf("Count is %d\n",count);
	FILE *f = fopen("avgtimecal.raw","w");
	for(chan=lowchan; chan<highchan; chan++) 
	{
		 fprintf(f,"%d %f %f %f %f\n",chan,Cxx[chan-lowchan],Cyy[chan-lowchan],Cxy[chan-lowchan],Cyx[chan-lowchan]);
	}
	fclose(f);

	diffusion_filter(Cxx, highchan - lowchan, fwindow);
	diffusion_filter(Cyy, highchan - lowchan, fwindow);
	diffusion_filter(Cxy, highchan - lowchan, fwindow);
	diffusion_filter(Cyx, highchan - lowchan, fwindow);

	float avgxx=0.0,avgyy=0.0,avgxy=0.0,avgyx=0.0;

	for(chan=lowchan; chan<highchan; chan++) 
	{
		avgxx += Cxx[chan-lowchan];
		avgyy += Cyy[chan-lowchan];
	}

	avgxx/=(highchan-lowchan);
	avgyy/=(highchan-lowchan);

	f = fopen("avgtimecal.dat","w");
	for(chan=lowchan; chan<highchan; chan++) 
	{
		 Cxx[chan-lowchan]/=avgxx;
		 Cyy[chan-lowchan]/=avgyy;
		 fprintf(f,"%d %f %f %f %f\n",chan,Cxx[chan-lowchan],Cyy[chan-lowchan],Cxy[chan-lowchan],Cyx[chan-lowchan]);
	}
	fclose(f);

	for(n=0; n<count; n++) 
	{
		
		Yxx[n] /=(highchan-lowchan);
		Yyy[n] /=(highchan-lowchan);
	}

	if(count)
	{
		chebyshev_minmax(XRA, count, &min, &max);
		chebyshev_normalize(XRA, count, min, max);
	}

	//XX
	C[0] = 0.0;
	C[1] = 0.0;
	if(count)
		chebyshev_fit_bw(XRA, Yxx, count, nsigma, C, order);
	printf("XX M %2.6f C %2.6f\n",C[1],C[0]);

	for(chan=lowchan; chan<highchan; chan++) 
	{
		for(n=0; n<size; n++) 
		{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
			{
				dataset[n].cal.xx[chan] = chebyshev_eval(CNORMALIZE(dataset[n].RA,min,max), C, order)*Cxx[chan-lowchan];
			}
		}
	}

	//YY
	C[0] = 0.0;
	C[1] = 0.0;
	if(count)
		chebyshev_fit_bw(XRA, Yyy, count, nsigma, C, order);		
	printf("YY M %2.6f C %2.6f\n",C[1],C[0]);
	for(chan=lowchan; chan<highchan; chan++) 
	{
		for(n=0; n<size; n++) 
		{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
			{
				dataset[n].cal.yy[chan] = chebyshev_eval(CNORMALIZE(dataset[n].RA,min,max), C, order)*Cyy[chan-lowchan];
			}
		}
	}

	//XY
	for(chan=lowchan; chan<highchan; chan++) 
	{
		for(n=0; n<size; n++) 
		{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
			{
				dataset[n].cal.xy[chan] = Cxy[chan-lowchan];
			}
		}
	}

	//YX
	for(chan=lowchan; chan<highchan; chan++) 
	{
		for(n=0; n<size; n++) 
		{
			if(dataset[n].flagBAD) continue;
			if(!RFIF || dataset[n].flagRFI[chan] == RFI_NONE) 
			{
				dataset[n].cal.yx[chan] = Cyx[chan-lowchan];
			}
		}
	}

	//printf("\n");
	free(XRA);
	free(Yxx);
	free(Yyy);
	free(Yxy);
	free(Yyx);
	free(Cxx);
	free(Cyy);
	free(Cxy);
	free(Cyx);
//	fclose(f);
	return;
}
//----------------------------------------------------------------------------------------------------------
//void smooth_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int window)
//void simple_smooth_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int twindow, int fwindow)
void simple_smooth_cal(SpecRecord dataset[], int start, int end, int records, int lowchan, int highchan, int twindow, int fwindow)
{
	//fix this loop need to run for all 4096 channels
	lowchan = 0; highchan = MAX_CHANNELS;

	int n, chan, dchan = highchan - lowchan;
	int size = end-start;
	float mean[4], tmp;
	static float mean_t[4] = {0};
	float sigma[4] = {0};
	float *XRA, *Yxx, *Yyy, *Yxy, *Yyx;
	static float *Fxx, *Fyy, *Fxy, *Fyx, *w;

	//Y has band avg series
	Yxx = (float*) malloc(sizeof(float) * size);
	Yyy = (float*) malloc(sizeof(float) * size);
	Yxy = (float*) malloc(sizeof(float) * size);
	Yyx = (float*) malloc(sizeof(float) * size);
	
	float **Sxx,**Syy,**Syx,**Sxy;

	//Fxx is time avg band shape
	if(0 == start)
	{
		Fxx = (float*) malloc(sizeof(float) * dchan);
		Fyy = (float*) malloc(sizeof(float) * dchan);
		Fxy = (float*) malloc(sizeof(float) * dchan);
		Fyx = (float*) malloc(sizeof(float) * dchan);
		if(Fxx == NULL || Fyy == NULL || Fxy == NULL || Fyx == NULL)
		{
			printf("malloc failed !\n");
			exit(0);
		}
	}

	if(Yxx == NULL || Yyy == NULL || Yxy == NULL || Yyx == NULL)
	{
		printf("malloc failed !\n");
		exit(0);
	}

        FILE * outfile2 = fopen("avgtimecal.dat","w");
        if(outfile2 == NULL)
        {
                printf("Can't open avgtimecal.dat\n");
                exit(0);
        }

	for(chan=lowchan; chan<highchan; chan++) 
	{
		mean[0] = mean[1] = mean[2] = mean[3] = 0;
		int count[MAX_CHANNELS] = {0};
		//for(n=0; n<size; n++) 
		for(n=start; n<end; n++) 
		{
			if(dataset[n].flagBAD) continue;
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

		//for(n=0; n<size; n++) 
		for(n=start; n<end; n++) 
		{
			//if(dataset[n].flagBAD || dataset[n].flagRFI[chan] != RFI_NONE)
			if(dataset[n].flagBAD)
			{
				dataset[n].cal.xx[chan] = mean[0];
				dataset[n].cal.yy[chan] = mean[1];
				dataset[n].cal.xy[chan] = mean[2];
				dataset[n].cal.yx[chan] = mean[3];
			}
		}
	}
	fclose(outfile2);

//	mean[0] = mean[1] = mean[2] = mean[3] = 0;

	int cnt = 0;
	if(0 == start)
	{
		float sumxx,sumyy,sumxy,sumyx;
		for(n=0; n<records; n++) 
		{
			if(dataset[n].flagBAD) continue;
			sumxx = 0.0,sumyy = 0.0,sumxy = 0.0;sumyx = 0.0;
			for(chan=lowchan; chan<highchan; chan++) 
			{
				sumxx += dataset[n].cal.xx[chan];
				sumyy += dataset[n].cal.yy[chan];
				sumxy += dataset[n].cal.xy[chan];
				sumyx += dataset[n].cal.yx[chan];
			}

			sumxx /= dchan; 
			sumyy /= dchan; 
			sumxy /= dchan; 
			sumyx /= dchan;
			mean_t[0] += sumxx;
			mean_t[1] += sumyy;
			mean_t[2] += sumxy;
			mean_t[3] += sumyx;
			cnt++;
		}

		mean_t[0] /= cnt; 
		mean_t[1] /= cnt; 
		mean_t[2] /= cnt; 
		mean_t[3] /= cnt;
	
	}

	cnt = 0;	
	for(n=start; n<end; n++) 
	{
		if(dataset[n].flagBAD) continue;
		Yxx[cnt] = Yyy[cnt] = Yxy[cnt] = Yyx[cnt] = 0.0;
		for(chan=lowchan; chan<highchan; chan++) 
		{
			Yxx[cnt] += dataset[n].cal.xx[chan];
			Yyy[cnt] += dataset[n].cal.yy[chan];
			Yxy[cnt] += dataset[n].cal.xy[chan];
			Yyx[cnt] += dataset[n].cal.yx[chan];
		}

		Yxx[cnt] /= dchan; 
		Yyy[cnt] /= dchan; 
		Yxy[cnt] /= dchan; 
		Yyx[cnt] /= dchan;

		sigma[0] += (Yxx[cnt] - mean_t[0])*(Yxx[cnt] - mean_t[0]);
		sigma[1] += (Yyy[cnt] - mean_t[1])*(Yyy[cnt] - mean_t[1]);
		sigma[2] += (Yxy[cnt] - mean_t[2])*(Yxy[cnt] - mean_t[2]);
		sigma[3] += (Yyx[cnt] - mean_t[3])*(Yyx[cnt] - mean_t[3]);

		cnt++;
	}

	sigma[0] = sqrt(sigma[0]/cnt);
	sigma[1] = sqrt(sigma[1]/cnt);
	sigma[2] = sqrt(sigma[2]/cnt);
	sigma[3] = sqrt(sigma[3]/cnt);

	for(n=0; n<cnt; n++)
	{
		if(fabs(Yxx[n] - mean_t[0]) > 3.0*sigma[0]) Yxx[n] = mean_t[0];
		if(fabs(Yyy[n] - mean_t[1]) > 3.0*sigma[1]) Yyy[n] = mean_t[1];
		if(fabs(Yxy[n] - mean_t[2]) > 3.0*sigma[2]) Yxy[n] = mean_t[2];
		if(fabs(Yyx[n] - mean_t[3]) > 3.0*sigma[3]) Yyx[n] = mean_t[3];
	}

	for(n=0; n<cnt; n++) 
	{
		Yxx[n] /= mean_t[0];
		Yyy[n] /= mean_t[1];
		Yxy[n] /= mean_t[2];
		Yyx[n] /= mean_t[3];
	}

	//To avoid window parameter from overflowing
	if(twindow > (cnt/2 - 1))
		twindow = cnt/2 -1;

	//apply the moving average filter to reduce noise
	moving_average_filter(Yxx, cnt, twindow);
	moving_average_filter(Yyy, cnt, twindow);
	moving_average_filter(Yxy, cnt, twindow);
	moving_average_filter(Yyx, cnt, twindow);

	for(chan=lowchan; chan<highchan; chan++)
 	{
		Fxx[chan - lowchan] = 0.0;
		Fyy[chan - lowchan] = 0.0;
		Fxy[chan - lowchan] = 0.0;
		Fyx[chan - lowchan] = 0.0;
	}

	int cc[MAX_CHANNELS] = {0};

	if(0 == start)
	{
		for(chan=lowchan; chan<highchan; chan++) 
		{
			for(n=0; n<records; n++) 
			{
				//if(!dataset[n].flagBAD && dataset[n].flagRFI[chan] == RFI_NONE)
				{
					Fxx[chan - lowchan] += dataset[n].cal.xx[chan];
					Fyy[chan - lowchan] += dataset[n].cal.yy[chan];
					Fxy[chan - lowchan] += dataset[n].cal.xy[chan];
					Fyx[chan - lowchan] += dataset[n].cal.yx[chan];
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

		diffusion_filter(Fxx, dchan, fwindow);
		diffusion_filter(Fyy, dchan, fwindow);
		diffusion_filter(Fxy, dchan, fwindow);
		diffusion_filter(Fyx, dchan, fwindow);

		write_band(Fxx,Fyy, Fxy, Fyx,n,lowchan,highchan);
	}

	if(0 == start)
	        outfile2 = fopen("smoothcal.dat","w");
	else
	        outfile2 = fopen("smoothcal.dat","a+");
        if(outfile2 == NULL)
        {
                printf("Can't open avgbandcal.dat\n");
        }
	cnt = 0;
        for(n=start; n<end; n++)
        {
		for(chan=lowchan; chan<highchan; chan++) 
		{
			dataset[n].cal.xx[chan] = Yxx[cnt]*Fxx[chan-lowchan];
			dataset[n].cal.yy[chan] = Yyy[cnt]*Fyy[chan-lowchan];
			dataset[n].cal.xy[chan] = Yxy[cnt]*Fxy[chan-lowchan];
			dataset[n].cal.yx[chan] = Yyx[cnt]*Fyx[chan-lowchan];
		}
		if(dataset[n].flagBAD)
			continue;
                fprintf(outfile2, "%05i %7.6f %7.6f %7.6f %7.6f\n",n,Yxx[cnt],Yyy[cnt],Yxy[cnt],Yyx[cnt]);
		cnt++;
	}
	fclose(outfile2);

	free(Yxx);
	free(Yyy);
	free(Yxy);
	free(Yyx);
	if(records == end)
	{
		free(Fxx);
		free(Fyy);
		free(Fxy);
		free(Fyx);
	}
	return;
}
//----------------------------------------------------------------------------------------------------------
void rolling_smooth_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int twindow, int maxiter)
{
	//fix this loop need to run for all 4096 channels
	lowchan = 0; highchan = MAX_CHANNELS;

	int n, chan, dchan = highchan - lowchan;
	float mean[4], sigma[4], tmp;
	float *XRA, *Yxx, *Yyy, *Yxy, *Yyx, *Fxx, *Fyy, *Fxy, *Fyx, *w;
	float *Fxxs, *Fyys, *Fxys, *Fyxs;

	//Y has band avg series
	Yxx = (float*) malloc(sizeof(float) * size);
	Yyy = (float*) malloc(sizeof(float) * size);
	Yxy = (float*) malloc(sizeof(float) * size);
	Yyx = (float*) malloc(sizeof(float) * size);
	
	float **Sxx,**Syy,**Syx,**Sxy;

	//S has smooth band avg series
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

	//Fxx is time avg band shape
	Fxx = (float*) malloc(sizeof(float) * dchan);
	Fyy = (float*) malloc(sizeof(float) * dchan);
	Fxy = (float*) malloc(sizeof(float) * dchan);
	Fyx = (float*) malloc(sizeof(float) * dchan);

	//smooth time avg band shape
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


        FILE * outfile2 = fopen("avgtimecal.dat","w");
        if(outfile2 == NULL)
        {
                printf("Can't open avgtimecal.dat\n");
                exit(0);
        }

	for(chan=lowchan; chan<highchan; chan++) 
	{
		mean[0] = mean[1] = mean[2] = mean[3] = 0;
		int count[MAX_CHANNELS] = {0};
		for(n=0; n<size; n++) 
		{
			if(dataset[n].flagBAD) continue;
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
		f(count[chan])
		{
			sigma[0] = sqrt(sigma[0]/count[chan]); 
			sigma[1] = sqrt(sigma[1]/count[chan]); 
			sigma[2] = sqrt(sigma[2]/count[chan]); 
			sigma[3] = sqrt(sigma[3]/count[chan]);
		}*/
		
		for(n=0; n<size; n++) 
		{
			//if(dataset[n].flagBAD || dataset[n].flagRFI[chan] != RFI_NONE)
			if(dataset[n].flagBAD)
			{
				dataset[n].cal.xx[chan] = mean[0];
				dataset[n].cal.yy[chan] = mean[1];
				dataset[n].cal.xy[chan] = mean[2];
				dataset[n].cal.yx[chan] = mean[3];
			}
		}
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
		if(fabs(Yxx[n] - mean[0]) > 2.5*sigma[0]) Yxx[n] = mean[0];
		if(fabs(Yyy[n] - mean[1]) > 2.5*sigma[1]) Yyy[n] = mean[1];
		if(fabs(Yxy[n] - mean[2]) > 2.5*sigma[2]) Yxy[n] = mean[2];
		if(fabs(Yyx[n] - mean[3]) > 2.5*sigma[3]) Yyx[n] = mean[3];
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
        diffusion_filter(Yxx, size, maxiter);
        diffusion_filter(Yyy, size, maxiter);
        diffusion_filter(Yxy, size, maxiter);
        diffusion_filter(Yyx, size, maxiter);

	//apply the moving average filter to reduce noise
	//moving_average_filter(Yxx, size, 500);
	//moving_average_filter(Yyy, size, 500);
	//moving_average_filter(Yxy, size, 500);
	//moving_average_filter(Yyx, size, 500);

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
				for(nn=0; nn<twindow; nn++) 
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
			else if(n < twindow/2 && n > 0)
			{
				;
			}
			else if(n >= twindow/2 && n < size - twindow/2)
			{
				int a = n - twindow/2;
				int b = n + twindow/2;
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

		//moving_average_filter(Fxxs, dchan, 36);
		//moving_average_filter(Fyys, dchan, 36);
		//moving_average_filter(Fxys, dchan, 36);
		//moving_average_filter(Fyxs, dchan, 36);

		diffusion_filter(Fxxs, dchan, 20);
		diffusion_filter(Fyys, dchan, 20);
		diffusion_filter(Fxys, dchan, 20);
		diffusion_filter(Fyxs, dchan, 20);

		for(chan=lowchan; chan<highchan; chan++) 
		{
			Sxx[n][chan] = Yxx[n]*Fxxs[chan-lowchan];
			Syy[n][chan] = Yyy[n]*Fyys[chan-lowchan];
			Sxy[n][chan] = Yxy[n]*Fxys[chan-lowchan];
			Syx[n][chan] = Yyx[n]*Fyxs[chan-lowchan];

		}
                if(!(n%10000))
                {
                        write_band(Fxxs,Fyys, Fxys, Fyxs,n,lowchan,highchan);
                }
	

	}

	printf("\n");

        outfile2 = fopen("smoothcal.dat","w");
        if(outfile2 == NULL)
        {
                printf("Can't open avgbandcal.dat\n");
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

/*        PolAvg avgcal;i
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
        sprintf(filename,"band_%04i.new",r);
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

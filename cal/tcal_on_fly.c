#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <glob.h>

#include "common.h"
#include "rfi.h"
#include "jsd_futil.h"
#include "denoising.h"
//--------------------------------------------------------------------------------------------------------
//int multibeam;

//--------------------------------------------------------------------------------------------------------
/*void compute_raw_cal(SpecRecord dataset[], int size, int lowchan, int highchan)
{
        int n, i;

        for(n=0; n<size; n++)
                {
                for(i=lowchan; i<highchan; i++) 
                        {
                                        dataset[n].cal.xx[i] = dataset[n].calon.xx[i] - dataset[n].caloff.xx[i];
                                        dataset[n].cal.yy[i] = dataset[n].calon.yy[i] - dataset[n].caloff.yy[i];
                                        //dataset[n].cal.xy[i] = dataset[n].calon.xy[i] - dataset[n].caloff.xy[i];
                                        //dataset[n].cal.yx[i] = dataset[n].calon.yx[i] - dataset[n].caloff.yx[i];
                        }
                        printf("%f %%\r", (n + 1)*100.0/size);
                }
}*/
//--------------------------------------------------------------------------------------------------------
void compute_tcal(SpecRecord dataset[], int size, int lowchan, int highchan, float hif, float hiband, float freq[], int *badchannels, float* tcalxx,float* tcalyy)
{
        int n, i;
        FILE * tcalfile = fopen("Tcal.dat","w");
        if(tcalfile == NULL)
        {
                printf("Can't open Tcal.dat\n");
                exit(1);
        }

        PolAvg avgcal[MAX_CHANNELS],avgdata[MAX_CHANNELS];
//	float tcalxx[MAX_CHANNELS],tcalyy[MAX_CHANNELS];
        int count =0;
        memset(&avgcal, 0,MAX_CHANNELS*sizeof(PolAvg));
        memset(&avgdata, 0, MAX_CHANNELS*sizeof(PolAvg));
        //memset(&tcalxx, 1.0, MAX_CHANNELS*sizeof(float));
        //memset(&tcalyy, 1.0, MAX_CHANNELS*sizeof(float));
	float minf=hif - hiband, maxf=hif + hiband;


        for(n=0; n<size; n++)
        {
		dataset[n].flagBAD = 0;
	}

        float mean[4], sigma[4];
	int chan;
        for(chan=lowchan; chan<highchan; chan++)
                {
		if(badchannels[chan])
			continue;

                mean[0] = mean[1] = mean[2] = mean[3] = 0;
                for(n=0; n<size; n++)
                        {
                        //mean[0] += dataset[n].cal.xx[chan];
                        mean[1] += dataset[n].cal.yy[chan];
                        //mean[2] += dataset[n].cal.xy[chan];
                        //mean[3] += dataset[n].cal.yx[chan];
                        mean[0] += dataset[n].caloff.xx[chan];
                        //mean[1] += dataset[n].caloff.yy[chan];
                        //mean[2] += dataset[n].caloff.xy[chan];
                        //mean[3] += dataset[n].caloff.yx[chan];i
			//printf("%f %f %f %f\n",dataset[n].cal.xx[chan],dataset[n].cal.yy[chan], dataset[n].caloff.xx[chan], dataset[n].caloff.yy[chan]);
                        }
                mean[0] /= size;
                mean[1] /= size;
                //mean[2] /= size;
                //mean[3] /= size;

                sigma[0] = sigma[1] = sigma[2] = sigma[3] = 0;
                for(n=0; n<size; n++)
                        {
                        //sigma[0] += (dataset[n].cal.xx[chan] - mean[0])*(dataset[n].cal.xx[chan] - mean[0]);
                        sigma[1] += (dataset[n].cal.yy[chan] - mean[1])*(dataset[n].cal.yy[chan] - mean[1]);
                        //sigma[2] += (dataset[n].cal.xy[chan] - mean[2])*(dataset[n].cal.xy[chan] - mean[2]);
                        //sigma[3] += (dataset[n].cal.yx[chan] - mean[3])*(dataset[n].cal.yx[chan] - mean[3]);
                        sigma[0] += (dataset[n].caloff.xx[chan] - mean[0])*(dataset[n].caloff.xx[chan] - mean[0]);
                        //sigma[1] += (dataset[n].caloff.yy[chan] - mean[1])*(dataset[n].caloff.yy[chan] - mean[1]);
                        //sigma[2] += (dataset[n].caloff.xy[chan] - mean[2])*(dataset[n].caloff.xy[chan] - mean[2]);
                        //sigma[3] += (dataset[n].caloff.yx[chan] - mean[3])*(dataset[n].caloff.yx[chan] - mean[3]);
                        }
                sigma[0] = sqrt(sigma[0]/size);
                sigma[1] = sqrt(sigma[1]/size);
                //sigma[2] = sqrt(sigma[2]/size);
                //sigma[3] = sqrt(sigma[3]/size);
		//printf("%d %f %f %f %f\n",chan,mean[0],mean[1],sigma[0],sigma[1]);

                for(n=0; n<size; n++)
                        {
                        //if(fabs(dataset[n].cal.xx[chan] - mean[0]) > 1.0*sigma[0]) dataset[n].flagBAD = 1;
                        if(fabs(dataset[n].cal.yy[chan] - mean[1]) > 4.0*sigma[1])  dataset[n].flagBAD = 1;
                        //if(fabs(dataset[n].cal.xy[chan] - mean[2]) > 1.0*sigma[2]) dataset[n].flagBAD = 1;
                        //if(fabs(dataset[n].cal.yx[chan] - mean[3]) > 1.0*sigma[3]) dataset[n].flagBAD = 1;
                        if(fabs(dataset[n].caloff.xx[chan] - mean[0]) > 2.5*sigma[0]) dataset[n].flagBAD = 1;
                        //if(fabs(dataset[n].caloff.yy[chan] - mean[1]) > 2.0*sigma[1]) dataset[n].flagBAD = 1;
                        //if(fabs(dataset[n].caloff.xy[chan] - mean[2]) > 2.0*sigma[2]) dataset[n].flagBAD = 1;
                        //if(fabs(dataset[n].caloff.yx[chan] - mean[3]) > 2.0*sigma[3]) dataset[n].flagBAD = 1;
                        }
                      //printf("%f %%\r", (chan - lowchan + 1)*100.0/(highchan-lowchan));
		}


        for(n=0; n<size; n++)
                {
        	if(dataset[n].flagBAD)
		{
			//printf("size %d\n",n);
			continue;
		}
                for(i=lowchan; i<highchan; i++)
                        {

                              if(dataset[n].flagRFI[i] == RFI_NONE)
                              {
                                        avgcal[i].xx += dataset[n].cal.xx[i];
                                        avgcal[i].yy += dataset[n].cal.yy[i];
                                        //avgcal[i].xy += dataset[n].cal.xy[i];
                                        //avgcal[i].yx += dataset[n].cal.yx[i];
                                        //avgdata[i].xx += dataset[n].calon.xx[i]+dataset[n].caloff.xx[i]-dataset[n].cal.xx[i];
                                        //avgdata[i].yy += dataset[n].calon.yy[i]+dataset[n].caloff.yy[i]-dataset[n].cal.yy[i];
                                        //avgdata[i].xy += dataset[n].calon.xy[i]+dataset[n].caloff.xy[i]-dataset[n].cal.xy[i];
                                        //avgdata[i].yx += dataset[n].calon.yx[i]+dataset[n].caloff.yx[i]-dataset[n].cal.yx[i];
                                        avgdata[i].xx += dataset[n].caloff.xx[i];
                                        avgdata[i].yy += dataset[n].caloff.yy[i];
                                        //avgdata[i].xy += dataset[n].caloff.xy[i];
                                        //avgdata[i].yx += dataset[n].caloff.yx[i];
					if(i == lowchan)
						count++;	
                              }
                        }
                        printf("%f %%\r", (n + 1)*100.0/size);
                }

	PolAvg avgavgcal, avgavgdata;
        memset(&avgavgcal, 0,sizeof(PolAvg));
        memset(&avgavgdata, 0,sizeof(PolAvg));
	int counter = 0;
        for (i=lowchan;i<highchan; i++)
        {
		if(!badchannels[i])
		{
		avgcal[i].xx /= count;
		avgcal[i].yy /= count;
		//avgcal[i].xy /= count;
		//avgcal[i].yx /= count;
		avgavgcal.xx += avgcal[i].xx;
		avgavgcal.yy += avgcal[i].yy;
		//avgavgcal.xy += avgcal[i].xy;
		//avgavgcal.yx += avgcal[i].yx;

		avgdata[i].xx /= count;
		avgdata[i].yy /= count;
		//avgdata[i].xy /= count;
		//avgdata[i].yx /= count;
		avgavgdata.xx += avgdata[i].xx;
		avgavgdata.yy += avgdata[i].yy;
		//avgavgdata.xy += avgdata[i].xy;
		//avgavgdata.yx += avgdata[i].yx;
		counter++;
		}
        }

	avgavgcal.xx /= counter;
	avgavgcal.xx /= counter;
	//avgavgcal.yy /= (lowchan-highchan);
	//avgavgcal.yy /= (lowchan-highchan);
	//avgavgcal.xy /= (lowchan-highchan);
	//avgavgcal.yx /= (lowchan-highchan);

	avgavgdata.xx /= counter;
	avgavgdata.xx /= counter;
	//avgavgdata.yy /= (lowchan-highchan);
	//avgavgdata.yy /= (lowchan-highchan);
	//avgavgdata.xy /= (lowchan-highchan);
	//avgavgdata.yx /= (lowchan-highchan);

	int min,max;
        for (i=0;i<MAX_CHANNELS; i++)
        {
		if(freq[i] > minf)
			max = i;
		if(freq[i] > maxf)
			min = i;
	}

	printf("%d %d\n",min,max);

        for (i=0;i<MAX_CHANNELS; i++)
        {

		if(i < lowchan)
		{
			//fprintf(tcalfile,"%2.6f %2.6f\n",1.0,1.0);
			tcalxx[i] = 1.0;
			tcalyy[i] = 1.0;
		}
		else if(badchannels[i] || i >= highchan)
		{
			tcalxx[i] = tcalxx[i-1];
			tcalyy[i] = tcalyy[i-1];
		}
		else
		{
			if(i > min && i < max)
			{       
				float xx,yy;
				tcalxx[i] = 0.5*((avgcal[min].xx*avgavgdata.xx)/(avgdata[min].xx*avgavgcal.xx)+\
				(avgcal[max].xx*avgavgdata.xx)/(avgdata[max].xx*avgavgcal.xx));                         
				tcalyy[i] = 0.5*((avgcal[min].yy*avgavgdata.yy)/(avgdata[min].yy*avgavgcal.yy)+\
				(avgcal[max].yy*avgavgdata.yy)/(avgdata[max].yy*avgavgcal.yy));                         
				//xx = 0.5*((avgcal[min].xx)/(avgavgcal.xx)+\
				(avgcal[max].xx)/(avgavgcal.xx));                         
				//yy = 0.5*((avgcal[min].yy)/(avgavgcal.yy)+\
				(avgcal[max].yy)/(avgavgcal.yy));                         
				//xx = 0.5*((avgcal[min].xy*avgavgdata.xy)/(avgdata[min].xy*avgavgcal.xy)+\
				(avgcal[max].xy*avgavgdata.xy)/(avgdata[max].xy*avgavgcal.xy));                         
				//yy = 0.5*((avgcal[min].yx*avgavgdata.yx)/(avgdata[min].yx*avgavgcal.yx)+\
				(avgcal[max].yx*avgavgdata.yx)/(avgdata[max].yx*avgavgcal.yx));                         
				//fprintf(tcalfile,"%2.6f %2.6f\n",tcal[i].xx,tcal[i].yy);
				//fprintf(tcalfile,"%2.6f %2.6f\n",xx,yy);
			}
			else
			{
				tcalxx[i] = (avgcal[i].xx*avgavgdata.xx)/(avgdata[i].xx*avgavgcal.xx);
				tcalyy[i] = (avgcal[i].yy*avgavgdata.yy)/(avgdata[i].yy*avgavgcal.yy);
				//fprintf(tcalfile,"%2.6f %2.6f\n",tcal[i].xx,tcal[i].yy);
				//fprintf(tcalfile,"%2.6f %2.6f\n",(avgcal[i].xx)/(avgavgcal.xx),\
				(avgcal[i].yy)/(avgavgcal.yy));
				//fprintf(tcalfile,"%2.6f %2.6f\n",(avgcal[i].xy*avgavgdata.xy)/(avgdata[i].xy*avgavgcal.xy),\
				(avgcal[i].yx*avgavgdata.yx)/(avgdata[i].yx*avgavgcal.yx));
			}
		}
	}

        moving_average_filter(tcalxx, MAX_CHANNELS, 40);
        moving_average_filter(tcalyy, MAX_CHANNELS, 40);
        //diffusion_filter(tcalxx, MAX_CHANNELS, 4000);		
        //diffusion_filter(tcalyy, MAX_CHANNELS, 4000);		

        for (i=0;i<MAX_CHANNELS; i=i+1)
        {
			fprintf(tcalfile,"%2.6f %2.6f\n",tcalxx[i],tcalyy[i]);
			//printf("%2.6f %2.6f\n",tcalxx[i],tcalyy[i]);
	}

        fclose(tcalfile);
}
//--------------------------------------------------------------------------------------------------------

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
void norm_one_tcal(int lowchan,int highchan,int *badchannels,float* tcalxx,float* tcalyy)
{
	int i;
	float avgxx=0.0,avgyy=0.0;
	int counter = 0;
        for (i=lowchan;i<highchan; i++)
        {
		if(!badchannels[i])
		{
			avgxx+=tcalxx[i];
			avgyy+=tcalyy[i];
			counter++;
		}
        }
	
	avgxx /= counter;
	avgyy /= counter;

        for (i=lowchan;i<highchan; i++)
        {
		if(!badchannels[i])
		{
			tcalxx[i]/=avgxx;
			tcalyy[i]/=avgyy;
		}
        }
        for (i=0;i<MAX_CHANNELS; i++)
        {
		if(i < lowchan)
		{
			tcalxx[i] = tcalxx[lowchan];
			tcalyy[i] = tcalyy[lowchan];
		}
		else if(badchannels[i])
		{
			tcalxx[i] = tcalxx[i-1];
			tcalyy[i] = tcalyy[i-1];
		}
		else if(i >= highchan)
		{
			tcalxx[i] = tcalxx[highchan-1];
			tcalyy[i] = tcalyy[highchan-1];
		}
		else
			continue;
	}
}
//--------------------------------------------------------------------------------------------------------
void write_tcal(float* tcalxx,float* tcalyy, int r)
{
	char filename[64];
	int i;
	sprintf(filename,"Tcalw_%04i.dat",r);
        FILE * tcalfile = fopen(filename,"w");
        if(tcalfile == NULL)
        {
                printf("Can't open Tcal.dat\n");
                exit(1);
        }
        for (i=0;i<MAX_CHANNELS;i++)
        {
			fprintf(tcalfile,"%2.6f %2.6f\n",tcalxx[i],tcalyy[i]);
			//printf("%2.6f %2.6f\n",tcalxx[i],tcalyy[i]);
	}

        fclose(tcalfile);

}
//--------------------------------------------------------------------------------------------------------
void compute_tcal(SpecRecord dataset[], int size, int lowchan, int highchan, float hif, float hiband, float freq[], int *badchannels, float* tcalxx,float* tcalyy, int rec, int window)
{
        int n, i;
        PolAvg avgcal[MAX_CHANNELS],avgdata[MAX_CHANNELS];
        static int count[MAX_CHANNELS] = {0};
	float minf=hif - hiband, maxf=hif + hiband;

        float mean[4], sigma[4];

	int min,max;
        for (i=0;i<MAX_CHANNELS; i++)
        {
		if(freq[i] > minf)
			max = i;
		if(freq[i] > maxf)
			min = i;
	}

	if(rec==0)
	{
        	memset(&avgcal, 0,MAX_CHANNELS*sizeof(PolAvg));
	        memset(&avgdata, 0, MAX_CHANNELS*sizeof(PolAvg));
	        for(n=0; n<window; n++)
       	        {
	       	 	if(dataset[n].flagBAD)
			{
			continue;
			}
        	        for(i=lowchan; i<highchan; i++)
                        {

                              if(dataset[n].flagRFI[i] == RFI_NONE)
                              {
                                        avgcal[i].xx += dataset[n].cal.xx[i];
                                        avgcal[i].yy += dataset[n].cal.yy[i];
                                        avgdata[i].xx += dataset[n].caloff.xx[i];
                                        avgdata[i].yy += dataset[n].caloff.yy[i];
						count[i]++;	
                              }
                        }
                }
	

	        for (i=lowchan;i<MAX_CHANNELS; i++)
	        {
			if(badchannels[i] || i >= highchan)
			{
				tcalxx[i] = tcalxx[i-1];
				tcalyy[i] = tcalyy[i-1];
			}
			else
			{
				if(i > min && i < max)
				{       
					if(count[i])
					{
						tcalxx[i] = 0.5*((avgcal[min].xx/avgdata[min].xx)+(avgcal[max].xx/avgdata[max].xx));
						tcalyy[i] = 0.5*((avgcal[min].yy/avgdata[min].yy)+(avgcal[max].yy/avgdata[max].yy));
					}
				}
				else
				{
					if(count[i])
					{
						tcalxx[i] = avgcal[i].xx/avgdata[i].xx;
						tcalyy[i] = avgcal[i].yy/avgdata[i].yy;
					}
				}
			}
		}
	}
	else if(rec < window/2 && rec > 0)
	{
		return;
	}
	else if(rec >= window/2 && rec < size - window/2)
	{
		int a = rec - window/2;
		int b = rec + window/2;
	       	if(!dataset[a].flagBAD)
		{
		        for (i=lowchan;i<MAX_CHANNELS; i++)
		        {
				if(badchannels[i] || i >= highchan)
				{
					tcalxx[i] = tcalxx[i-1];
					tcalyy[i] = tcalyy[i-1];
				}
				else
				{
					if(i > min && i < max)
					{       
						if(count[min])
						{
						tcalxx[i] -= (0.5/count[min])*((dataset[a].cal.xx[min]/dataset[a].caloff.xx[min])+(dataset[a].cal.xx[max]/dataset[a].caloff.xx[max]));
						tcalyy[i] -= (0.5/count[min])*((dataset[a].cal.yy[min]/dataset[a].caloff.yy[min])+(dataset[a].cal.yy[max]/dataset[a].caloff.yy[max]));
						tcalxx[i]/=(count[min]-1);
						tcalyy[i]/=(count[min]-1);
						}
					}
					else
					{	
						if(count[i])
						{
						tcalxx[i] -= (1.0/count[i])*(dataset[a].cal.xx[i]/dataset[a].caloff.xx[i]);
						tcalyy[i] -= (1.0/count[i])*(dataset[a].cal.yy[i]/dataset[a].caloff.yy[i]);
						tcalxx[i]/=(count[i]-1);
						tcalyy[i]/=(count[i]-1);
						}
					}
					
				}
				if(count[i])
					count[i]--;
			}		
			//count--;
		}

	       	if(!dataset[b].flagBAD)
		{
			//count++;

		        for (i=lowchan;i<MAX_CHANNELS; i++)
		        {
				count[i]++;
				if(badchannels[i] || i >= highchan)
				{
					tcalxx[i] = tcalxx[i-1];
					tcalyy[i] = tcalyy[i-1];
				}
				else
				{
					if(i > min && i < max)
					{       
						if(count[min])
						{
						tcalxx[i] = tcalxx[i]*((count[min]-1)/count[min]);
						tcalyy[i] = tcalyy[i]*((count[min]-1)/count[min]);
						tcalxx[i] += (0.5/count[min])*((dataset[b].cal.xx[min]/dataset[b].caloff.xx[min])+(dataset[b].cal.xx[max]/dataset[b].caloff.xx[max]));
						tcalyy[i] += (0.5/count[min])*((dataset[b].cal.yy[min]/dataset[b].caloff.yy[min])+(dataset[b].cal.yy[max]/dataset[b].caloff.yy[max]));
						}
					}
					else
					{
						if(count[i])
						{
						tcalxx[i] = tcalxx[i]*((count[i]-1)/count[i]);
						tcalyy[i] = tcalyy[i]*((count[i]-1)/count[i]);
						tcalxx[i] += (1.0/count[i])*(dataset[b].cal.xx[i]/dataset[b].caloff.xx[i]);
						tcalyy[i] += (1.0/count[i])*(dataset[b].cal.yy[i]/dataset[b].caloff.yy[i]);
						}
					}
				}
			}
		}
	}
	else
	{
		return;
	}
}
//--------------------------------------------------------------------------------------------------------

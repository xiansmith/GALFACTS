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
		if(!badchannels[i] && isfinite(tcalxx[i]) && isfinite(tcalyy[i]) && !isnan(tcalxx[i]) && !isnan(tcalyy[i]))
		{
			avgxx+=tcalxx[i];
			avgyy+=tcalyy[i];
			counter++;
		}
        }

	if(counter)
	{	
		avgxx /= counter;
		avgyy /= counter;
	}
	else
	{
		avgxx = 1.0;
		avgyy = 1.0;
	}
//	printf("avg %f %f\n",avgxx,avgyy);

	if(!isfinite(avgxx) || !isfinite(avgyy) || isnan(avgxx) || isnan(avgyy))
		printf("averages messed up \n");
        //for (i=lowchan;i<highchan; i++)
        for (i=0;i<MAX_CHANNELS; i++)
        {
		if(!badchannels[i] && isfinite(tcalxx[i]) && isfinite(tcalyy[i]) && !isnan(tcalxx[i]) && !isnan(tcalyy[i]))
		{
			tcalxx[i]/=avgxx;
			tcalyy[i]/=avgyy;
		}
		else
		{
			tcalxx[i] = 1.0;
			tcalyy[i] = 1.0;
		}
        }
        for (i=0;i<MAX_CHANNELS; i++)
        {
		/*if(i < lowchan)
		{
			tcalxx[i] = tcalxx[lowchan];
			tcalyy[i] = tcalyy[lowchan];
		}
		else if(i >= highchan)
		{
			tcalxx[i] = tcalxx[highchan-1];
			tcalyy[i] = tcalyy[highchan-1];
		}*/
		//else 
		if(badchannels[i] || !isfinite(tcalxx[i]) || !isfinite(tcalyy[i]) || isnan(tcalxx[i]) || isnan(tcalyy[i]))
		{
			if(i == 0)
			{
				if(!isfinite(tcalxx[i]) || isnan(tcalxx[i]))
					tcalxx[i] = 1.0;
				if(!isfinite(tcalyy[i]) || isnan(tcalyy[i]))
					tcalyy[i] = 1.0;
			}
			else
			{
				if(badchannels[i-1] || !isfinite(tcalxx[i-1]) || isnan(tcalxx[i-1]) || !isfinite(tcalyy[i-1]) || isnan(tcalyy[i-1]))
				{
					tcalxx[i] = 1.0;
					tcalyy[i] = 1.0;
				}
				else
				{
					tcalxx[i] = tcalxx[i-1];
					tcalyy[i] = tcalyy[i-1];
				}
				if(!isfinite(tcalxx[i]) || !isfinite(tcalyy[i]) || isnan(tcalxx[i]) || isnan(tcalyy[i]))
         	                	printf("Still NaNs found  :-( %f %f\n",tcalxx[i],tcalyy[i]);
			}
		}
		else
		{
			continue;
		}
	}
/*        for (i=0;i<MAX_CHANNELS; i++)
        {
		if(!isfinite(tcalxx[i]) || !isfinite(tcalyy[i]) || !isnan(tcalxx[i]) || !isnan(tcalyy[i]))
			printf("NaNs found\n");
	}
*/	return;
}
//--------------------------------------------------------------------------------------------------------
void write_tcal(float* tcalxx,float* tcalyy, int r, int low, int high)
{
	char filename[64];
	int i;
	sprintf(filename,"Tcalw_%04i.new",r);
        FILE * tcalfile = fopen(filename,"w");
        if(tcalfile == NULL)
        {
                printf("Can't open Tcal.dat\n");
                exit(1);
        }
        //for (i=0;i<MAX_CHANNELS;i++)
        for (i=low;i<high;i++)
        {
			fprintf(tcalfile,"%2.6f %2.6f\n",tcalxx[i],tcalyy[i]);
			//printf("%2.6f %2.6f\n",tcalxx[i],tcalyy[i]);
	}

        fclose(tcalfile);
	return;

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
        		for (i=0;i<MAX_CHANNELS; i++)
                        {

                              if(dataset[n].flagRFI[i] == RFI_NONE && !badchannels[i])
                              {
                                        avgcal[i].xx += dataset[n].cal.xx[i];
                                        avgcal[i].yy += dataset[n].cal.yy[i];
                                        avgdata[i].xx += dataset[n].caloff.xx[i];
                                        avgdata[i].yy += dataset[n].caloff.yy[i];
					count[i]++;	
                              }
                        }
                }

                //for(i=lowchan; i<highchan; i++)
        	for (i=0;i<MAX_CHANNELS; i++)
                {
			//if(i >= 3400 && i < 3800)
			//printf("Chan %d count %d dataxx %f\n",i,count[i],avgcal[i].xx/avgdata[i].xx);

                        if(dataset[n].flagRFI[i] == RFI_NONE && !badchannels[i])
                        {
				if(count[i])
				{
                                avgcal[i].xx /=count[i];
                                avgcal[i].yy /=count[i];
                                avgdata[i].xx /=count[i];
                                avgdata[i].yy /=count[i];
				}
                        }
                }

	        //for (i=lowchan;i<highchan; i++)
        	for (i=0;i<MAX_CHANNELS; i++)
	        {
			//if(badchannels[i] || i >= highchan)
			if(!badchannels[i])
			{
				if(i > min && i < max)
				{       
					if(count[min])
					{
						tcalxx[i] = 0.5*((avgcal[min].xx/avgdata[min].xx)+(avgcal[max].xx/avgdata[max].xx));
						tcalyy[i] = 0.5*((avgcal[min].yy/avgdata[min].yy)+(avgcal[max].yy/avgdata[max].yy));
					}
					else
					{
						tcalxx[i] = tcalxx[i-1];
		                                tcalyy[i] = tcalyy[i-1];
					}
				}
				else
				{
					//if(count[min])
					if(count[i])
					{
						tcalxx[i] = avgcal[i].xx/avgdata[i].xx;
						tcalyy[i] = avgcal[i].yy/avgdata[i].yy;
					}
					else
					{
						tcalxx[i] = tcalxx[i-1];
		                                tcalyy[i] = tcalyy[i-1];
					}
				}
			}
			else
			{
				tcalxx[i] = tcalxx[i-1];
				tcalyy[i] = tcalyy[i-1];
			}
		}
		//write_tcal(tcalxx,tcalyy, 0, lowchan,highchan);
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
	        	//for (i=lowchan;i<highchan; i++)
        		for (i=0;i<MAX_CHANNELS; i++)
		        {
				if(badchannels[i])
				{
					tcalxx[i] = tcalxx[i-1];
					tcalyy[i] = tcalyy[i-1];
				}
				else if(dataset[a].flagRFI[i] == RFI_NONE)
				{
					if(i > min && i < max)
					{       
						if(count[min+1])
						{
						tcalxx[i] = tcalxx[i]*(count[min+1]) - (0.5)*((dataset[a].cal.xx[min+1]/dataset[a].caloff.xx[min+1])+(dataset[a].cal.xx[max]/dataset[a].caloff.xx[max]));
						tcalyy[i] = tcalyy[i]*(count[min+1]) - (0.5)*((dataset[a].cal.yy[min+1]/dataset[a].caloff.yy[min+1])+(dataset[a].cal.yy[max]/dataset[a].caloff.yy[max]));
						tcalxx[i] = tcalxx[i]/(count[min+1]-1.0);
						tcalyy[i] = tcalyy[i]/(count[min+1]-1.0);
						
						}
					}
					else
					{	
						if(count[i] >= 2)
						{
						tcalxx[i] -= (1.0/count[i])*(dataset[a].cal.xx[i]/dataset[a].caloff.xx[i]);
						tcalyy[i] -= (1.0/count[i])*(dataset[a].cal.yy[i]/dataset[a].caloff.yy[i]);
						tcalxx[i] = tcalxx[i] *(count[i]/(count[i]-1.0));
						tcalyy[i] = tcalyy[i] *(count[i]/(count[i]-1.0));
						count[i]--;
						}
					}
					
				}
			}		
		}

	       	if(!dataset[b].flagBAD)
		{
        		for (i=0;i<MAX_CHANNELS; i++)
		        {
				if(badchannels[i])
				{
					tcalxx[i] = tcalxx[i-1];
					tcalyy[i] = tcalyy[i-1];
				}
				else if(dataset[b].flagRFI[i] == RFI_NONE)
				{
					count[i]++;
					if(i > min && i < max)
					{       
						if(count[min+1])
						{
						tcalxx[i] = tcalxx[i]*((count[min+1]-1.0)/count[min+1]);
						tcalyy[i] = tcalyy[i]*((count[min+1]-1.0)/count[min+1]);
						tcalxx[i] += (0.5/count[min+1])*((dataset[b].cal.xx[min+1]/dataset[b].caloff.xx[min+1])+(dataset[b].cal.xx[max]/dataset[b].caloff.xx[max]));
						tcalyy[i] += (0.5/count[min+1])*((dataset[b].cal.yy[min+1]/dataset[b].caloff.yy[min+1])+(dataset[b].cal.yy[max]/dataset[b].caloff.yy[max]));
						}
					}
					else
					{
						if(count[i])
						{
				//printf("Chan %d count %d befdataxx %f ",i,count[i],tcalxx[i]);
						tcalxx[i] = tcalxx[i]*((count[i]-1.0)/count[i]);
						tcalyy[i] = tcalyy[i]*((count[i]-1.0)/count[i]);
						tcalxx[i] = tcalxx[i] + (1.0/count[i])*(dataset[b].cal.xx[i]/dataset[b].caloff.xx[i]);
						tcalyy[i] = tcalyy[i] + (1.0/count[i])*(dataset[b].cal.yy[i]/dataset[b].caloff.yy[i]);
				//printf("aftdataxx %f bad %d\n",tcalxx[i],badchannels[i]);
						}
					}
				}
			}
		}

		/*	if(rec == 2000)	
			{
			printf("========REC 2000=============");
        		for (i=3400;i<3800; i++)
				printf("Chan %d count %d dataxx %f bad %d\n",i,count[i],tcalxx[i],badchannels[i]);*/
        		/*for (i=3400;i<3800; i++)
				//printf("Chan %d count %d dataxx %f\n",i,count[i],(1.0/count[i])*(dataset[b].cal.xx[i]/dataset[b].caloff.xx[i]));
				printf("Chan %d count %d dataxx %f add %f subt %f \n",i,count[i],tcalxx[i],(1.0/count[i])*(dataset[b].cal.xx[i]/dataset[b].caloff.xx[i]),(1.0/count[i])*(dataset[a].cal.xx[i]/dataset[a].caloff.xx[i]));
			}*/
	}
	else
	{
		return;
	}
	return;
}
//--------------------------------------------------------------------------------------------------------

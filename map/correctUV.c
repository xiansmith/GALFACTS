#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "correctUV.h"

void correct_UV(FluxWappData * wappdata, int chan, MapMetaData *md)
{
	//currently only works for full channel and not channel averages 
	//this is because the table contains MAX_CHANNEL values
	//For channel averages should we average the table values ??

	int d,i;
	for(d=0; d<wappdata->numDays; d++)
        {
		FluxDayData * daydata = &wappdata->daydata[d];

		//read the table
		//float epsilon[MAX_CHANNELS], phi[MAX_CHANNELS];
		float Uleak[MAX_CHANNELS], Vleak[MAX_CHANNELS];
		FILE *epsphi;
		char filename[64];
		sprintf(filename,"UVleakage%d.dat",d%7);
		epsphi = fopen(filename,"r");
		if(epsphi == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
			//fscanf(epsphi,"%f %f",&epsilon[i],&phi[i]);
			fscanf(epsphi,"%f %f",&Uleak[i],&Vleak[i]);
			}
			fclose(epsphi);

			if(md->avg)
			{	
				int j;
				sprintf(filename,"UVleakage%d_avg%04i.dat",d%7,chan);
				epsphi = fopen(filename,"w");
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					//int chn = j;
					for(k = j+1;k<j+md->avg;k++)
					{
						Uleak[j]+=Uleak[k];
						Vleak[j]+=Vleak[k];
					}
					Uleak[j]/=md->avg;
					Vleak[j]/=md->avg;
				}
				fprintf(epsphi,"%d %f %f\n",j,Uleak[j],Vleak[j]);
				fclose(epsphi);
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
					Uleak[0]+=Uleak[k];
					Vleak[0]+=Vleak[k];
				}
				Uleak[0]/=(md->avg_highchan - md->avg_lowchan);
				Vleak[0]/=(md->avg_highchan - md->avg_lowchan);

			}

			printf("INFO: read file %s\n", filename);
		}

		//apply the corrections
		int r = daydata->numRecords;
		for(i=0;i<r;i++)
		{
			//daydata->records[i].stokes.U -= 2*epsilon[chan]*daydata->records[i].stokes.I*cos(phi[chan]);
			//daydata->records[i].stokes.V -= 2*epsilon[chan]*daydata->records[i].stokes.I*sin(phi[chan]);
			daydata->records[i].stokes.U -= daydata->records[i].stokes.I*Uleak[chan];
			daydata->records[i].stokes.V -= daydata->records[i].stokes.I*Vleak[chan];
		}
        }
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "spec_corr.h"

void spec_corr(FluxWappData * wappdata, int chan, MapMetaData *md)
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
		float g[7][MAX_CHANNELS];
		int c;
		FILE *Icorr;
		char filename[64];
		sprintf(filename,"Icorr.dat");

		Icorr = fopen(filename,"r");
		if(Icorr == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
				fscanf(Icorr ,"%d %f",&c,&g[0][i]);
				fscanf(Icorr ,"%f",&g[1][i]);
				fscanf(Icorr ,"%f",&g[2][i]);
				fscanf(Icorr ,"%f",&g[3][i]);
				fscanf(Icorr ,"%f",&g[4][i]);
				fscanf(Icorr ,"%f",&g[5][i]);
				fscanf(Icorr ,"%f",&g[6][i]);
			}
			fclose(Icorr);

			if(md->avg)
			{	
				int j;
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					for(k = j+1;k<j+md->avg;k++)
					{
						g[0][j]+=g[0][k];
						g[1][j]+=g[1][k];
						g[2][j]+=g[2][k];
						g[3][j]+=g[3][k];
						g[4][j]+=g[4][k];
						g[5][j]+=g[5][k];
						g[6][j]+=g[6][k];
					}
					g[0][j]/=md->avg;
					g[1][j]/=md->avg;
					g[2][j]/=md->avg;
					g[3][j]/=md->avg;
					g[4][j]/=md->avg;
					g[5][j]/=md->avg;
					g[6][j]/=md->avg;
				}
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
						g[0][0]+=g[0][k];
						g[1][0]+=g[1][k];
						g[2][0]+=g[2][k];
						g[3][0]+=g[3][k];
						g[4][0]+=g[4][k];
						g[5][0]+=g[5][k];
						g[6][0]+=g[6][k];
				}
				g[0][0]/=(md->avg_highchan - md->avg_lowchan);
				g[1][0]/=(md->avg_highchan - md->avg_lowchan);
				g[2][0]/=(md->avg_highchan - md->avg_lowchan);
				g[3][0]/=(md->avg_highchan - md->avg_lowchan);
				g[4][0]/=(md->avg_highchan - md->avg_lowchan);
				g[5][0]/=(md->avg_highchan - md->avg_lowchan);
				g[6][0]/=(md->avg_highchan - md->avg_lowchan);
			}

		}

		//apply the corrections
		int r = daydata->numRecords;
		int beam;
                if(!strcmp(wappdata->wapp,"multibeam"))
                {
                        beam = d%7;
                }
                else
                {
                        beam = atoi(&wappdata->wapp[4]);
                }

		for(i=0;i<r;i++)
		{
			daydata->records[i].stokes.I *= g[beam][chan];
		}
        }
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pacorr.h"

void pa_corr(FluxWappData * wappdata, int chan, MapMetaData *md)
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
		float g[7][MAX_CHANNELS] ,p[7][MAX_CHANNELS];
		int c;
		FILE *pcorr;
		char filename[64];
		sprintf(filename,"Pcorr.dat");

		pcorr = fopen(filename,"r");
		if(pcorr == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
				fscanf(pcorr ,"%d %f %f",&c,&g[0][i],&p[0][i]);
				fscanf(pcorr ,"%f %f",&g[1][i],&p[1][i]);
				fscanf(pcorr ,"%f %f",&g[2][i],&p[2][i]);
				fscanf(pcorr ,"%f %f",&g[3][i],&p[3][i]);
				fscanf(pcorr ,"%f %f",&g[4][i],&p[4][i]);
				fscanf(pcorr ,"%f %f",&g[5][i],&p[5][i]);
				fscanf(pcorr ,"%f %f",&g[6][i],&p[6][i]);
			}
			fclose(pcorr);

			if(md->avg)
			{	
				int j;
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					for(k = j+1;k<j+md->avg;k++)
					{
						g[0][j]+=g[0][k];
						p[0][j]+=p[0][k];
						g[1][j]+=g[1][k];
						p[1][j]+=p[1][k];
						g[2][j]+=g[2][k];
						p[2][j]+=p[2][k];
						g[3][j]+=g[3][k];
						p[3][j]+=p[3][k];
						g[4][j]+=g[4][k];
						p[4][j]+=p[4][k];
						g[5][j]+=g[5][k];
						p[5][j]+=p[5][k];
						g[6][j]+=g[6][k];
						p[6][j]+=p[6][k];
					}
					g[0][j]/=md->avg;
					p[0][j]/=md->avg;
					g[1][j]/=md->avg;
					p[1][j]/=md->avg;
					g[2][j]/=md->avg;
					p[2][j]/=md->avg;
					g[3][j]/=md->avg;
					p[3][j]/=md->avg;
					g[4][j]/=md->avg;
					p[4][j]/=md->avg;
					g[5][j]/=md->avg;
					p[5][j]/=md->avg;
					g[6][j]/=md->avg;
					p[6][j]/=md->avg;
				}
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
						g[0][0]+=g[0][k];
						p[0][0]+=p[0][k];
						g[1][0]+=g[1][k];
						p[1][0]+=p[1][k];
						g[2][0]+=g[2][k];
						p[2][0]+=p[2][k];
						g[3][0]+=g[3][k];
						p[3][0]+=p[3][k];
						g[4][0]+=g[4][k];
						p[4][0]+=p[4][k];
						g[5][0]+=g[5][k];
						p[5][0]+=p[5][k];
						g[6][0]+=g[6][k];
						p[6][0]+=p[6][k];
				}
				g[0][0]/=(md->avg_highchan - md->avg_lowchan);
				p[0][0]/=(md->avg_highchan - md->avg_lowchan);
				g[1][0]/=(md->avg_highchan - md->avg_lowchan);
				p[1][0]/=(md->avg_highchan - md->avg_lowchan);
				g[2][0]/=(md->avg_highchan - md->avg_lowchan);
				p[2][0]/=(md->avg_highchan - md->avg_lowchan);
				g[3][0]/=(md->avg_highchan - md->avg_lowchan);
				p[3][0]/=(md->avg_highchan - md->avg_lowchan);
				g[4][0]/=(md->avg_highchan - md->avg_lowchan);
				p[4][0]/=(md->avg_highchan - md->avg_lowchan);
				g[5][0]/=(md->avg_highchan - md->avg_lowchan);
				p[5][0]/=(md->avg_highchan - md->avg_lowchan);
				g[6][0]/=(md->avg_highchan - md->avg_lowchan);
				p[6][0]/=(md->avg_highchan - md->avg_lowchan);
				//printf("averaged Uleak is %f Vleak is %f\n", Uleak[0], Vleak[0] );
			}

			//printf("INFO: read file %s\n", filename);
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

		p[beam][chan] = p[beam][chan]*M_PI/180.0;
		for(i=0;i<r;i++)
		{
			float Q_old = daydata->records[i].stokes.Q;
			float U_old = daydata->records[i].stokes.U;
			daydata->records[i].stokes.Q = g[beam][chan]*(Q_old*cos(2.0*p[beam][chan])-U_old*sin(2.0*p[beam][chan]));
			daydata->records[i].stokes.U = g[beam][chan]*(Q_old*sin(2.0*p[beam][chan])+U_old*cos(2.0*p[beam][chan]));
		}
        }
}

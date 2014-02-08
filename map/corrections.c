#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "corrections.h"

//contains all the corrections derived from the CALIBRATOR sources. Since the multiple steps in the correction all rely on original data its better to deal with them all at a same place. Otherwise it is messy to pass the same data back and forth between different functions

void corrections(FluxWappData * wappdata, int chan, MapMetaData *md)
{
	//currently only works for full channel and not channel averages 
	//this is because the table contains MAX_CHANNEL values
	//For channel averages should we average the table values ??

	int d,i,b;
	for(d=0; d<wappdata->numDays; d++)
        {
		FluxDayData * daydata = &wappdata->daydata[d];

		//read the tables
		float qc[7][MAX_CHANNELS];
		float uc[7][MAX_CHANNELS];
		float vc[7][MAX_CHANNELS];
		float Pc[7][MAX_CHANNELS];
		float spec_I[7][MAX_CHANNELS];
		float pa[7][MAX_CHANNELS];

		int c;
		FILE *table;
		char filename[64];
		sprintf(filename,"Qcorr.dat");

		table = fopen(filename,"r");
		if(table == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
				fscanf(table ,"%d",&c);
				//printf("%d ",c);
				for(b=0;b<7;b++)
				{
					fscanf(table ,"%f",&qc[b][i]);
					//printf("%f ",qc[b][i]);
				}
				//printf("\n");
			}
			fclose(table);

			if(md->avg)
			{	
				int j;
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					for(k = j+1;k<j+md->avg;k++)
					{
						for(b=0;b<7;b++)
						{
							qc[b][j]+=qc[b][k];
						}
					}
					for(b=0;b<7;b++)
					{
						qc[b][j]/=md->avg;
					}
				}
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
					for(b=0;b<7;b++)
					{
						qc[b][0]+=qc[b][k];
					}
				}
				for(b=0;b<7;b++)
				{
					qc[b][0]/=(md->avg_highchan - md->avg_lowchan);
				}
			}

		}

		sprintf(filename,"Ucorr.dat");
		table = fopen(filename,"r");
		if(table == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
				fscanf(table ,"%d",&c);
				//printf("%d ",c);
				for(b=0;b<7;b++)
				{
					fscanf(table ,"%f",&uc[b][i]);
					//printf("%f ",uc[b][i]);
				}
				//printf("\n");
			}
			fclose(table);

			if(md->avg)
			{	
				int j;
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					for(k = j+1;k<j+md->avg;k++)
					{
						for(b=0;b<7;b++)
						{
							uc[b][j]+=uc[b][k];
						}
					}
					for(b=0;b<7;b++)
					{
						uc[b][j]/=md->avg;
					}
				}
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
					for(b=0;b<7;b++)
					{
						uc[b][0]+=uc[b][k];
					}
				}
				for(b=0;b<7;b++)
				{
					uc[b][0]/=(md->avg_highchan - md->avg_lowchan);
				}
			}

		}

		sprintf(filename,"Vcorr.dat");
		table = fopen(filename,"r");
		if(table == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
				fscanf(table ,"%d",&c);
				//printf("%d ",c);
				for(b=0;b<7;b++)
				{
					fscanf(table ,"%f",&vc[b][i]);
					//printf("%f ",vc[b][i]);
				}
				//printf("\n");
			}
			fclose(table);

			if(md->avg)
			{	
				int j;
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					for(k = j+1;k<j+md->avg;k++)
					{
						for(b=0;b<7;b++)
						{
							vc[b][j]+=vc[b][k];
						}
					}
					for(b=0;b<7;b++)
					{
						vc[b][j]/=md->avg;
					}
				}
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
					for(b=0;b<7;b++)
					{
						vc[b][0]+=vc[b][k];
					}
				}
				for(b=0;b<7;b++)
				{
					vc[b][0]/=(md->avg_highchan - md->avg_lowchan);
				}
			}

		}

		sprintf(filename,"Icorr.dat");
		table = fopen(filename,"r");
		if(table == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
				fscanf(table ,"%d",&c);
				//printf("%d ",c);
				for(b=0;b<7;b++)
				{
					fscanf(table ,"%f",&spec_I[b][i]);
					//printf("%f ",spec_I[b][i]);
				}
				//printf("\n");
			}
			fclose(table);

			if(md->avg)
			{	
				int j;
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					for(k = j+1;k<j+md->avg;k++)
					{
						for(b=0;b<7;b++)
						{
							spec_I[b][j]+=spec_I[b][k];
						}
					}
					for(b=0;b<7;b++)
					{
						spec_I[b][j]/=md->avg;
					}
				}
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
					for(b=0;b<7;b++)
					{
						spec_I[b][0]+=spec_I[b][k];
					}
				}
				for(b=0;b<7;b++)
				{
					spec_I[b][0]/=(md->avg_highchan - md->avg_lowchan);
				}
			}

		}

		/* not needed for now can be enabled later
		sprintf(filename,"Pcorr.dat");
		table = fopen(filename,"r");
		if(table == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
				fscanf(table ,"%d",&c);
				for(b=0;b<7;b++)
				{
					fscanf(table ,"%f",&Pc[b][i]);
				}
			}
			fclose(table);

			if(md->avg)
			{	
				int j;
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					for(k = j+1;k<j+md->avg;k++)
					{
						for(b=0;b<7;b++)
						{
							Pc[b][j]+=Pc[b][k];
						}
					}
					for(b=0;b<7;b++)
					{
						Pc[b][j]/=md->avg;
					}
				}
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
					for(b=0;b<7;b++)
					{
						Pc[b][0]+=Pc[b][k];
					}
				}
				for(b=0;b<7;b++)
				{
					Pc[b][0]/=(md->avg_highchan - md->avg_lowchan);
				}
			}

		}

		sprintf(filename,"PAcorr.dat");
		table = fopen(filename,"r");
		if(table == NULL)
		{
			printf("ERROR: unable to open file %s\n", filename);
			return;
		}
		else
		{
			for(i = 0; i < MAX_CHANNELS; i++)
			{
				fscanf(table ,"%d",&c);
				for(b=0;b<7;b++)
				{
					fscanf(table ,"%f",&pa[b][i]);
				}
			}
			fclose(table);

			if(md->avg)
			{	
				int j;
				for(j = md->avg_lowchan;j < md->avg_highchan;j+=md->avg)
				{
					int k;
					for(k = j+1;k<j+md->avg;k++)
					{
						for(b=0;b<7;b++)
						{
							pa[b][j]+=pa[b][k];
						}
					}
					for(b=0;b<7;b++)
					{
						pa[b][j]/=md->avg;
					}
				}
			}
			if( chan == 0 ) {
				int k = 0;
				for(k = md->avg_lowchan;k< md->avg_highchan;k++)
				{
					for(b=0;b<7;b++)
					{
						pa[b][0]+=pa[b][k];
					}
				}
				for(b=0;b<7;b++)
				{
					pa[b][0]/=(md->avg_highchan - md->avg_lowchan);
				}
			}

		}
		*/
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
		
		//pa[beam][chan] = pa[beam][chan]*M_PI/180.0;

		//plugging in these for now
		pa[beam][chan] = -M_PI/4.0;
		Pc[beam][chan] = 1.0;
		for(i=0;i<r;i++)
		{
			daydata->records[i].stokes.Q = daydata->records[i].stokes.Q - qc[beam][chan]*daydata->records[i].stokes.I;
			daydata->records[i].stokes.U = daydata->records[i].stokes.U - uc[beam][chan]*daydata->records[i].stokes.I;
			daydata->records[i].stokes.V = daydata->records[i].stokes.V - vc[beam][chan]*daydata->records[i].stokes.I;
			daydata->records[i].stokes.I *= spec_I[beam][chan];
			float Q = daydata->records[i].stokes.Q;
			float U = daydata->records[i].stokes.U;
			daydata->records[i].stokes.Q = Pc[beam][chan]*(Q*cos(2*pa[beam][chan])-U*sin(2*pa[beam][chan]));
			daydata->records[i].stokes.U = Pc[beam][chan]*(Q*sin(2*pa[beam][chan])+U*cos(2*pa[beam][chan]));
		}
        }
}

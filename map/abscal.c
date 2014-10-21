#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "abscal.h"
#include "chebyshev.h"
//----------------------------------------------------------------------------------------------------------------------------------------
static void absolute_cal_day(FluxWappData * wappdata, int day, int beam, int chan, MapMetaData *md)
{
	FILE *table;
	char filename[64];
	float spec_I[7][MAX_CHANNELS];
	int c,i,b;
/*	sprintf(filename,"Icorr.dat");
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
*/


	FluxDayData * daydata = &wappdata->daydata[day];
	int R = daydata->numRecords;
	int r;
	for (r = 0; r < R; r++)
	{
		if (isfinite(daydata->records[r].stokes.I))
		{
			float RA = daydata->records[r].RA;
			//S1
			//daydata->records[r].stokes.I -= (0.0000479748017*RA*RA-0.004185405*RA-3.65278283);
			//S2
			//daydata->records[r].stokes.I -= (0.0000330741985*RA*RA-0.0116655548*RA-2.41301002);
			//N4
			daydata->records[r].stokes.I -= (0.0000000957942322*RA*RA+0.00644514212*RA-5.84658068);
		}
	}
	return;
}

//-------------------------------------------------------------------------------
void abs_cal(FluxWappData * wappdata, int chan, MapMetaData *md)
{
int d,beam;
for(d=0; d<wappdata->numDays; d++) 
	{
        if(!strcmp(wappdata->wapp,"multibeam"))
	{
		beam = d%7;
	}
	else
		beam = atoi(&wappdata->wapp[4]);

	absolute_cal_day(wappdata,d,beam,chan,md);
	}
}
//-------------------------------------------------------------------------------

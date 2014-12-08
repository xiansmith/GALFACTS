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
	int i,b;
	float c[3];
        sprintf(filename,"%s.abscal.txt",md->field);

	table = fopen(filename,"r");
	if(table == NULL)
	{
		printf("ERROR: unable to open file %s\n", filename);
		return;
	}
	else
	{
		for(i = 0; i < 3; i++)
		{
			fscanf(table ,"%f",&c[i]);
		}
		fclose(table);


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
				//daydata->records[r].stokes.I -= (0.0000000957942322*RA*RA+0.00644514212*RA-5.84658068);
				daydata->records[r].stokes.I -= (c[2]*RA*RA+c[1]*RA+c[0]);
			}
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "correctUV.h"

void correct_UV(FluxWappData * wappdata, int chan, MapMetaData *md)
{
	//currently only works for full channel and not channel averages 
	//this is because the table contains MAX_CHANNEL values
	//For channel averages should we average the table values ??

	int d;
	int i = 0, j = 0, k = 0;

	float Ulk, Vlk;
//	float Qlk[7];
	int cnt;

	for(d=0; d<wappdata->numDays; d++)
        {
		FluxDayData * daydata = &wappdata->daydata[d];

		int r = daydata->numRecords;
		Ulk = Vlk = 0.0;
		cnt = 0;

		for(i=0;i<r;i++)
		{
			if( daydata->records[i].RA < 202.70 || daydata->records[i].RA > 202.80 || daydata->records[i].DEC < 30.42 ||  daydata->records[i].DEC > 30.58)
			{
			Ulk =  Ulk + daydata->records[i].stokes.U/daydata->records[i].stokes.I;
			Vlk =  Vlk + daydata->records[i].stokes.V/daydata->records[i].stokes.I;
			//Qlk[b] =  Qlk[b] + daydata->records[i].stokes.Q/daydata->records[i].stokes.I;
			cnt++;
			}
                }

		if(cnt)
		{
			Ulk/=cnt;
			Vlk/=cnt;
		}

		for(i=0;i<r;i++)
		{
			daydata->records[i].stokes.U -=  Ulk*daydata->records[i].stokes.I;
			daydata->records[i].stokes.V -=  Vlk*daydata->records[i].stokes.I;
                }
	}

}

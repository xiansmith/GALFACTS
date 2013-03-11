#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "correctUV.h"

void correct_UV(FluxWappData * wappdata, int chan)
{
	//currently only works for full channel and not channel averages 
	//this is because the table contains MAX_CHANNEL values
	//For channel averages should we average the table values ??

	int d,i;
	for(d=0; d<wappdata->numDays; d++)
        {
		FluxDayData * daydata = &wappdata->daydata[d];

		//read the table
		float epsilon[MAX_CHANNELS], phi[MAX_CHANNELS];
		FILE *epsphi;
		char filename[64];
		sprintf(filename,"epsilon_phi%d.dat",d%7);
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
				fscanf(epsphi,"%f %f",&epsilon[i],&phi[i]);
			}
			fclose(epsphi);
			printf("INFO: read file %s\n", filename);
		}

		//apply the corrections
		int r = daydata->numRecords;
		for(i=0;i<r;i++)
		{
			daydata->records[i].stokes.U -= 2*epsilon[chan]*daydata->records[i].stokes.I*cos(phi[chan]);
			daydata->records[i].stokes.V -= 2*epsilon[chan]*daydata->records[i].stokes.I*sin(phi[chan]);
		}
        }
}


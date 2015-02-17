#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_pointing.h"
#define ARECIBO_ZENITH 18.34424

//-------------------------------------------------------------------------------
void fix_pointing(const char *field, FluxWappData *wappdata)
{
	int r, day;
//	printf("Field %s\n",field);
	
	FILE *fits;
	char filename[100];
	sprintf(filename,"%s.pointing_fits.dat",field);
	fits = fopen(filename,"r");

	if(NULL == fits)
		printf("ERR: Could not open file %s\n");

	float low_zenith[7][3];
	float high_zenith[7][3];
	float raerr[7];
	int i,j;
	for(i=1;i<7;i++)
	{
		//skip over beam 0 since the corrections are for the outer beams.
		//the beam 0 record in the array sits empty.
		for(j=0;j<3;j++)
			fscanf(fits,"%f",&low_zenith[i][j]);
		for(j=0;j<3;j++)
			fscanf(fits,"%f",&high_zenith[i][j]);
		fscanf(fits,"%f",&raerr[i]);
		raerr[i]/=3600.0;

		//debug print messages
/*		for(j=0;j<3;j++)
			printf("%f ",low_zenith[i][j]);
		printf("\n");
		for(j=0;j<3;j++)
			printf("%f ",high_zenith[i][j]);
		printf("\n");
		printf("%f",raerr[i]*3600.0);
		printf("\n");
*/	}	

	for(day=0; day<wappdata->numDays; day++)
	{
		//printf("Day %d\n",day);
		//if(day > 7)
		//	exit(1);
		FluxDayData * daydata = &wappdata->daydata[day];
		for(r=0; r<daydata->numRecords; r++)
		{
			float decerr;
			float ZEN = ARECIBO_ZENITH - daydata->records[r].DEC;
			int beam = day%7;

			switch(beam)
			{
				case 1:
					if(ZEN < 14.0)
						decerr = (low_zenith[beam][0] + low_zenith[beam][1]*ZEN + low_zenith[beam][2]*ZEN*ZEN)/3600.0;
					else
						decerr = (high_zenith[beam][0] + high_zenith[beam][1]*ZEN + high_zenith[beam][2]*ZEN*ZEN)/3600.0;
					daydata->records[r].DEC -= decerr;
					daydata->records[r].RA -= raerr[beam];
					//printf("%d RA %f %f DEC %f %f\n",day%7,daydata->records[r].RA,daydata->records[r].RA-raerr,daydata->records[r].DEC,daydata->records[r].DEC-decerr);
					break;
				case 2:
					if(ZEN < 14.0)
						decerr = (low_zenith[beam][0] + low_zenith[beam][1]*ZEN + low_zenith[beam][2]*ZEN*ZEN)/3600.0;
					else
						decerr = (high_zenith[beam][0] + high_zenith[beam][1]*ZEN + high_zenith[beam][2]*ZEN*ZEN)/3600.0;
					daydata->records[r].DEC -= decerr;
					daydata->records[r].RA -= raerr[beam];
					//printf("%d RA %f %f DEC %f %f\n",day%7,daydata->records[r].RA,daydata->records[r].RA-raerr,daydata->records[r].DEC,daydata->records[r].DEC-decerr);
					break;
				case 3:
					if(ZEN < 14.0)
						decerr = (low_zenith[beam][0] + low_zenith[beam][1]*ZEN + low_zenith[beam][2]*ZEN*ZEN)/3600.0;
					else
						decerr = (high_zenith[beam][0] + high_zenith[beam][1]*ZEN + high_zenith[beam][2]*ZEN*ZEN)/3600.0;
					daydata->records[r].DEC -= decerr;
					daydata->records[r].RA -= raerr[beam];
					//printf("%d RA %f %f DEC %f %f\n",day%7,daydata->records[r].RA,daydata->records[r].RA-raerr,daydata->records[r].DEC,daydata->records[r].DEC-decerr);
					break;
				case 4:
					if(ZEN < 14.0)
						decerr = (low_zenith[beam][0] + low_zenith[beam][1]*ZEN + low_zenith[beam][2]*ZEN*ZEN)/3600.0;
					else
						decerr = (high_zenith[beam][0] + high_zenith[beam][1]*ZEN + high_zenith[beam][2]*ZEN*ZEN)/3600.0;
					daydata->records[r].DEC -= decerr;
					daydata->records[r].RA -= raerr[beam];
					//printf("%d RA %f %f DEC %f %f\n",day%7,daydata->records[r].RA,daydata->records[r].RA-raerr,daydata->records[r].DEC,daydata->records[r].DEC-decerr);
					break;
				case 5:
					if(ZEN < 14.0)
						decerr = (low_zenith[beam][0] + low_zenith[beam][1]*ZEN + low_zenith[beam][2]*ZEN*ZEN)/3600.0;
					else
						decerr = (high_zenith[beam][0] + high_zenith[beam][1]*ZEN + high_zenith[beam][2]*ZEN*ZEN)/3600.0;
					daydata->records[r].DEC -= decerr;
					daydata->records[r].RA -= raerr[beam];
					//printf("%d RA %f %f DEC %f %f\n",day%7,daydata->records[r].RA,daydata->records[r].RA-raerr,daydata->records[r].DEC,daydata->records[r].DEC-decerr);
					break;
				case 6:
					if(ZEN < 14.0)
						decerr = (low_zenith[beam][0] + low_zenith[beam][1]*ZEN + low_zenith[beam][2]*ZEN*ZEN)/3600.0;
					else
						decerr = (high_zenith[beam][0] + high_zenith[beam][1]*ZEN + high_zenith[beam][2]*ZEN*ZEN)/3600.0;
					daydata->records[r].DEC -= decerr;
					daydata->records[r].RA -= raerr[beam];
					//printf("%d RA %f %f DEC %f %f\n",day%7,daydata->records[r].RA,daydata->records[r].RA-raerr,daydata->records[r].DEC,daydata->records[r].DEC-decerr);
					break;
				default:
					break;
			}
		}
	}
}

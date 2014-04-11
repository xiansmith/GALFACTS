#include <stdio.h>
#include <stdlib.h>
#include "denoising.h"
#include <math.h>
#include "programs/fitsio.h"
#include "common.h"

//#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )

static void print_usage(char * argv[])
{
        printf("\n");
        printf("Usage: %s <filename> <lowchan> <highchan>\n", argv[0]);
}

int main(int argc, char *argv[])
{
	FILE *f;
	if(argc !=4)
	{
		print_usage(argv);
		return(EXIT_FAILURE);
	}
	char *filename = argv[1];
	int lowchan = atoi(argv[2]);
	int highchan = atoi(argv[3]);

	f = fopen(filename,"r");
	if(NULL == f)
        {
                printf("ERROR: file %s not found.\n",filename);
                exit(1);
        }

	float b0[MAX_CHANNELS],b1[MAX_CHANNELS],b2[MAX_CHANNELS],b3[MAX_CHANNELS],b4[MAX_CHANNELS],b5[MAX_CHANNELS],b6[MAX_CHANNELS];

	int i;
	//for(i=0;i<MAX_CHANNELS;i++)
	for(i=0;i<highchan-lowchan;i++)
	{
		fscanf(f,"%f %f %f %f %f %f %f",&b0[i],&b1[i],&b2[i],&b3[i],&b4[i],&b5[i],&b6[i]);
		if(isnan(b0[i]))
		{
			printf("found nan\n");
			b0[i] = b0[i-1];
		}
		//if(isnan(y[i]))
		//{
	//		printf("found nan\n");
//			y[i] = y[i-1];
//		}
	}
	fclose(f);
	moving_average_filter(b0,highchan-lowchan,100);
	moving_average_filter(b1,highchan-lowchan,100);
	moving_average_filter(b2,highchan-lowchan,100);
	moving_average_filter(b3,highchan-lowchan,100);
	moving_average_filter(b4,highchan-lowchan,100);
	moving_average_filter(b5,highchan-lowchan,100);
	moving_average_filter(b6,highchan-lowchan,100);

	char outfile[30];
	sprintf(outfile,"GALFACTS_calibration.smooth");
	f = fopen(outfile,"w");
	//for(i=0;i<MAX_CHANNELS;i++)
	for(i=0;i<highchan-lowchan;i++)
	{
		//fprintf(f,"%f %f\n",x[i],y[i]);
		fprintf(f,"%f %f %f %f %f %f %f\n",b0[i],b1[i],b2[i],b3[i],b4[i],b5[i],b6[i]);
	}
	fclose(f);
}

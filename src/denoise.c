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
        printf("Usage: %s <filename> <beam>\n", argv[0]);
}

int main(int argc, char *argv[])
{
	FILE *f;
	if(argc !=3)
	{
		print_usage(argv);
		return(EXIT_FAILURE);
	}
	char *filename = argv[1];
	int beam = atoi(argv[2]);
	
	f = fopen(filename,"r");
	if(NULL == f)
        {
                printf("ERROR: file %s not found.\n",filename);
                exit(1);
        }

	float x[MAX_CHANNELS],y[MAX_CHANNELS];

	int i;
	for(i=0;i<MAX_CHANNELS;i++)
	{
		fscanf(f,"%f %f",&x[i],&y[i]);
		if(isnan(x[i]))
		{
			printf("found nan\n");
			x[i] = x[i-1];
		}
		if(isnan(y[i]))
		{
			printf("found nan\n");
			y[i] = y[i-1];
		}
	}
	fclose(f);
	moving_average_filter(x,MAX_CHANNELS,100);
	moving_average_filter(y,MAX_CHANNELS,100);

	char outfile[30];
	sprintf(outfile,"epsilon_phi%d.smooth",beam);
	f = fopen(outfile,"w");
	for(i=0;i<MAX_CHANNELS;i++)
	{
		fprintf(f,"%f %f\n",x[i],y[i]);
	}
	fclose(f);
}

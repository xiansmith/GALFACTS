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
	//printf("Usage: %s <cube prefix> <beam> <maxp> <minp> <lowchan> <highchan>\n", argv[0]);
	printf("Usage: %s <Icube> <Vcube> <frac>\n", argv[0]);
}

int main(int argc, char *argv[]) 
{

	int j;
	char *fname1, *fname2;
	FILE *Icube,*Vcube;
	header_param_list Ihpar,Vhpar;
	float *Iplane_data,*Vplane_data;
	int plane_size;
	float frac;
	if (argc != 4) {
		print_usage(argv);
		return EXIT_FAILURE;
	}
	fname1 = argv[1];
	fname2 = argv[2];
	frac = atof(argv[3]);
	Icube = fopen(fname1,"r");
        if(Icube == NULL)
        {
                printf("ERROR: cube file %s not found.\n",fname1);
		exit(1);
	}
	readfits_header(Icube,&Ihpar);

	Vcube = fopen(fname2,"r");
        if(Vcube == NULL)
        {
                printf("ERROR: cube file %s not found.\n",fname2);
		exit(1);
	}
	readfits_header(Vcube,&Vhpar);

	plane_size = Ihpar.naxis[0] * Ihpar.naxis[1];
	Iplane_data = malloc(sizeof(float)*plane_size);
	Vplane_data = malloc(sizeof(float)*plane_size);
	
	readfits_plane(Icube,Iplane_data,&Ihpar);
	readfits_plane(Vcube,Vplane_data,&Vhpar);
	float RA,DEC;
	for(j=0;j<plane_size;j++)
	{
		float pct;
		if(Iplane_data[j]!=0.0)
			pct = fabs(Vplane_data[j]/Iplane_data[j]);
		if(pct > frac)
		{
			RA = (j%Ihpar.naxis[0] - Ihpar.crpix[0])*Ihpar.cdelt[0] + Ihpar.crval[0];
			DEC = (j/Ihpar.naxis[0] - Ihpar.crpix[1])*Ihpar.cdelt[1]+ Ihpar.crval[1];
			printf("CIRCLE %2.4f %2.4f %2.4f\n",RA,DEC,fabs(Vplane_data[j]/Iplane_data[j])*1);
		}	
	}

	fclose(Icube);
	fclose(Vcube);

	return EXIT_SUCCESS;
}

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
	printf("Usage: %s <cube prefix> <lowchan> <highchan>\n", argv[0]);
}

int main(int argc, char *argv[]) 
{

	int i, n;
	int maxx,maxy,rad,minp;
	char *prefix, outfilename[130], cubename[130];
	FILE *Icube,*Ucube,*Vcube,*outfile;
	header_param_list Ihpar,Uhpar,Vhpar;
	float *plane_data;
	int plane_size;
	int num_planes;
	int num_comp;
	int lowchan,highchan;
	int beam;

	if (argc != 4) {
		print_usage(argv);
		return EXIT_FAILURE;
	}
	prefix = argv[1];
	//maxx = atoi(argv[3]);
	//maxy = atoi(argv[4]);
	//minp = atoi(argv[5]);
	//rad = atoi(argv[5]);
	lowchan = atoi(argv[2]);
	highchan = atoi(argv[3]);

	sprintf(cubename,"%s_%d_%04i_10chanavg_V.fits",prefix,lowchan,highchan);
	Icube = fopen(cubename,"r");
        if(Icube == NULL)
        {
                printf("ERROR: cube file %s not found.\n",cubename);
		exit(1);
	}
	readfits_header(Icube,&Ihpar);

	printf("---------------------Here 1-----------------\n");

	plane_size = Ihpar.naxis[0] * Ihpar.naxis[1];
	num_planes = Ihpar.naxis[2];
	plane_data = malloc(sizeof(float)*plane_size);
	
	if((highchan-lowchan)/10!=num_planes)
	{
		printf("ERROR: highchan - lowchan is not equal to the number of planes in this cube.\n");
		exit(1);
	}

	float *sI,*sU,*sV;
	sI = (float *)malloc(num_planes*sizeof(float));
	sU = (float *)malloc(num_planes*sizeof(float));
	sV = (float *)malloc(num_planes*sizeof(float));

	printf("num_planes: %i\n", num_planes);
	printf("plane_size: %i\n", plane_size);

	int *coords[2];
	coords[0] = malloc(sizeof(int)* num_comp);
	coords[1] = malloc(sizeof(int)* num_comp);


	float back;
	int j,count;
	double rms;
	double mean;
	for(i=0;i<num_planes;i++)
	{	
		mean = 0.0;
		rms = 0.0;
		count = 0;
		readfits_plane(Icube,plane_data,&Ihpar);
		for(j=0;j<plane_size;j++)
		{
			if(plane_data[j] < 100.0 && plane_data[j] > -100.0)
			{
			//rms+=(plane_data[j]*plane_data[j]);	
			mean+=(plane_data[j]);
			count++;
			}	
		}
		if(count)	
		mean/=count;
		count++;
		for(j=0;j<plane_size;j++)
		{
			if(plane_data[j] < 100.0 && plane_data[j] > -100.0)
			{
			rms+=(plane_data[j]-mean)*(plane_data[j]-mean);	
			count++;
			}	
		}
		if(count)
		rms/=count;
		printf("%d %2.6f\n",lowchan+i*10,sqrt(rms));
	}
	exit(1);
	printf("---------------------Here 2-----------------\n");
	
	//sprintf(cubename,"%s_BEAM%d_%04i_%04i_Ucube.fits",prefix,beam,lowchan,highchan-1);
	sprintf(cubename,"%s_%d_%04i_10chanavg_V.fits",prefix,lowchan,highchan);
	Ucube = fopen(cubename,"r");
        if(Ucube == NULL)
        {
                printf("ERROR: U cube file not found.\n");
		exit(1);
	}
	//sprintf(cubename,"%s_BEAM%d_%04i_%04i_Vcube.fits",prefix,beam,lowchan,highchan-1);
	sprintf(cubename,"%s_%d_%04i_10chanavg_V.fits",prefix,lowchan,highchan);
	Vcube = fopen(cubename,"r");
        if(Vcube == NULL)
        {
                printf("ERROR: V cube file not found.\n");
		exit(1);
	}
	readfits_header(Ucube,&Uhpar);
	readfits_header(Vcube,&Vhpar);
	float *epsilon,*phi;

	epsilon = (float *)malloc(num_planes*sizeof(float));
	phi = (float *)malloc(num_planes*sizeof(float));

	//sprintf(outfilename,"epsilon_phi%d.dat",beam);

	//outfile = fopen(outfilename,"w");
	printf("---------------------Here 3-----------------\n");

	for(i=0;i<num_planes;i++)
	{
		rms = 0.0;
		readfits_plane(Ucube,plane_data,&Uhpar);
		for(j=0;j<plane_size;j++)
		{
			if(plane_data[j] < 100.0 && plane_data[j] > -100.0)
			rms+=(plane_data[j]*plane_data[j]);	
		}
		rms/=plane_size;
		printf("Stokes U CChannel: %d RMS:%2.6f K\n",lowchan+i,sqrt(rms));
		rms = 0.0;
		readfits_plane(Vcube,plane_data,&Vhpar);
		for(j=0;j<plane_size;j++)
		{
			if(plane_data[j] < 100.0 && plane_data[j] > -100.0)
			rms+=(plane_data[j]*plane_data[j]);	
		}
		rms/=plane_size;
		printf("Stokes I Channel: %d RMS:%2.6f K\n",lowchan+i,sqrt(rms));
	}

	fclose(Icube);
	fclose(Ucube);
	fclose(Vcube);
	//fclose(outfile);

	return EXIT_SUCCESS;
}

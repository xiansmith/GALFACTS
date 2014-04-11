#include <stdio.h>
#include <stdlib.h>
//#include "denoising.h"
#include <math.h>
#include "programs/fitsio.h"
#include "common.h"

//#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )

static void print_usage(char * argv[])
{
	printf("\n");
	//printf("Usage: %s <cube prefix> <beam> <maxp> <minp> <lowchan> <highchan>\n", argv[0]);
	printf("Usage: %s <cube prefix> <beam>  <lowchan> <highchan>\n", argv[0]);
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
//	beam = atoi(argv[2]);
	//maxx = atoi(argv[2]);
	//maxy = atoi(argv[3]);
	//minp = atoi(argv[5]);
	//rad = atoi(argv[4]);
	lowchan = atoi(argv[2]);
	highchan = atoi(argv[3]);

	float *Ibeam[7];
	float Iavg[7];
	for(beam=0;beam<7;beam++)
	{
		Ibeam[beam] = (float*)malloc(sizeof(float)*MAX_CHANNELS);
		for(i=0;i<MAX_CHANNELS;i++)
			Ibeam[beam][i] = 1.0;
	}

	for(beam=0;beam<7;beam++)
	{
		Iavg[beam] = 0.0;
		sprintf(cubename,"%s_BEAM%d_%04i_%04i_I.fits",prefix,beam,lowchan,highchan-1);
		Icube = fopen(cubename,"r");
	        if(Icube == NULL)
	        {
	                printf("ERROR: cube file %s not found.\n",cubename);
			exit(1);
		}
		readfits_header(Icube,&Ihpar);
	

		plane_size = Ihpar.naxis[0] * Ihpar.naxis[1];
		num_planes = Ihpar.naxis[2];
		plane_data = malloc(sizeof(float)*plane_size);
		
		if((highchan-lowchan)!=num_planes)
		{
			printf("ERROR: highchan - lowchan is not equal to the number of planes in this cube.\n");
			exit(1);
		}
		
//		printf("num_planes: %i\n", num_planes);
//		printf("plane_size: %i\n", plane_size);
	
		int j,count =0;
		int max;
//		int max = 1830;
//		int min = 270;
		for(i=0;i<num_planes;i++)
		{
			readfits_plane(Icube,plane_data,&Ihpar);
			if(i==0)
				{
				float peak = 0.0;
				for(j=0;j<plane_size;j++)
				{
					if(peak < plane_data[j])
					{
						peak = plane_data[j];
						max = j;
					}
				}
			
				//printf("Beam %d Peak: %f Back: %f Pixel %d\n",beam,plane_data[max],plane_data[min],max);
				printf("Beam %d Peak: %f Pixel %d\n",beam,plane_data[max],max);
			}
			
			//printf("Beam %d Peak: %f Pixel %d\n",beam,plane_data[max],max);
			if(plane_data[max]< 100.0)
			{
				//Ibeam[beam][lowchan+i] = plane_data[max]-plane_data[min];
				//Iavg[beam] += (plane_data[max]-plane_data[min]);
				Ibeam[beam][lowchan+i] = plane_data[max];
				Iavg[beam] += plane_data[max];
				count++;
			}
		}
		
		
		/*for(i=lowchan;i<highchan;i++)
		{
			Iavg[beam]+=Ibeam[beam][i];
		}*/
		//printf("Iavg%d% 2.8f\n",beam,Iavg[beam]);
		Iavg[beam]/=count;
		fclose(Icube);
		free(plane_data);
	}
	sprintf(outfilename,"beamgains.dat");
	outfile = fopen(outfilename,"w");

	for(i=0;i<MAX_CHANNELS;i++)
	{
		for(beam=0;beam<7;beam++)
		{
			fprintf(outfile,"%2.8f ",Ibeam[beam][i]/Ibeam[0][i]);
		}
		fprintf(outfile,"\n");
	}

	for(beam=0;beam<7;beam++)
	{
		fprintf(outfile,"%2.8f ",Iavg[beam]/Iavg[0]);
		printf("%2.8f ",Iavg[beam]/Iavg[0]);
	}
	fprintf(outfile,"\n");
	printf("\n");
	fclose(outfile);

/*        for(n=0; n<num_planes; n++) x[n] = sI[n];
        diffusion_filter(x, num_planes, 1000);
        for(n=0; n<num_planes; n++) sI[n] = x[n];*/

	return EXIT_SUCCESS;
}


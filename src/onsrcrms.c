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
	printf("Usage: %s <cube prefix> <maxx> <maxy> <lowchan> <highchan>\n", argv[0]);
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

	if (argc != 6) {
		print_usage(argv);
		return EXIT_FAILURE;
	}
	prefix = argv[1];
	//beam = atoi(argv[2]);
	maxx = atoi(argv[2]);
	maxy = atoi(argv[3]);
	//minp = atoi(argv[5]);
	//rad = atoi(argv[5]);
	lowchan = atoi(argv[4]);
	highchan = atoi(argv[5]);

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
	double rms = 0.0;
	double mean = 0.0;
	for(i=0;i<num_planes;i++)
	{
		readfits_plane(Icube,plane_data,&Ihpar);
		//if(rad == -1)
		//{
			if(plane_data[maxy*Ihpar.naxis[0]+maxx] > -100.0 && plane_data[maxy*Ihpar.naxis[0]+maxx] < 100.0)
			{
			mean += plane_data[maxy*Ihpar.naxis[0]+maxx];
		//	
			sV[i] = plane_data[maxy*Ihpar.naxis[0]+maxx];
			count++;
			}
			else
			sV[i] = plane_data[maxy*Ihpar.naxis[0]+maxx];
			//	sV[i] = 0.0;

		//else
		//{
		/*	count = 0;
			back = 0.0;
			for(j=0;j<plane_size;j++)
			{
				if((j%Ihpar.naxis[0] < (maxx-rad) || j%Ihpar.naxis[0] > (maxx+rad))\
				||(j/Ihpar.naxis[0] < (maxy-rad) || j/Ihpar.naxis[0] > (maxy+rad)) )
				{
					back += plane_data[i];
					count++;
				}
			}
			back/=count;*/
			//sI[i] = plane_data[maxy*Ihpar.naxis[0]+maxx]-plane_data[rad];
			//sI[i] = plane_data[maxy*Ihpar.naxis[0]+maxx]-back;
		//}
		//if(i==0)
		//	printf("Peak: %f Background: %f\n",plane_data[maxy*Ihpar.naxis[0]+maxx],plane_data[minp]);
			//printf("Peak: %f Background: %f\n",plane_data[maxy*Ihpar.naxis[0]+maxx],back);
		//sI[i] = (plane_data[1770]+plane_data[1769]+plane_data[1771])/3.0-(plane_data[570]+plane_data[1750]+plane_data[1790]+plane_data[3030])/4.0;
		//sI[i] = plane_data[1770]-(plane_data[570]+plane_data[1750]+plane_data[1790]+plane_data[3030])/4.0;
	}

	if(count)
		mean/=count;
	count = 0;
	for(i=0;i<num_planes;i++)
	{
			//if(plane_data[maxy*Ihpar.naxis[0]+maxx] > -100.0 && plane_data[maxy*Ihpar.naxis[0]+maxx] < 100.0)
			if(sV[i] > -100.0 && sV[i] < 100.0)
			{
			rms+=(sV[i]-mean)*(sV[i]-mean);
			count++;
			}
	}
	if(count)
	rms/=count;
	printf("Rms: %f\n",rms);
	exit(1);
	printf("---------------------Here 2-----------------\n");
	
	sprintf(cubename,"%s_BEAM%d_%04i_%04i_Ucube.fits",prefix,beam,lowchan,highchan-1);
	Ucube = fopen(cubename,"r");
        if(Ucube == NULL)
        {
                printf("ERROR: U cube file not found.\n");
		exit(1);
	}
	sprintf(cubename,"%s_BEAM%d_%04i_%04i_Vcube.fits",prefix,beam,lowchan,highchan-1);
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

	sprintf(outfilename,"epsilon_phi%d.dat",beam);

	outfile = fopen(outfilename,"w");
	printf("---------------------Here 3-----------------\n");

	for(i=0;i<num_planes;i++)
	{
		readfits_plane(Ucube,plane_data,&Uhpar);
		if(rad == -1)
		{
			sU[i] = plane_data[maxy*Ihpar.naxis[0]+maxx];
		}
		else
		{
                        /*count = 0;
                        back = 0.0;
                        for(j=0;j<plane_size;j++)
                        {
                                if((j%Uhpar.naxis[0] < (maxx-rad) || j%Ihpar.naxis[0] > (maxx+rad))\
                                || (j/Uhpar.naxis[0] < (maxy-rad) || j/Ihpar.naxis[0] > (maxy+rad)) )
                                {
                                        back += plane_data[i];
                                        count++;
                                }
                        }
                        back/=count;*/

			//sU[i] = plane_data[maxy*Uhpar.naxis[0]+maxx]-back;
			sU[i] = plane_data[maxy*Uhpar.naxis[0]+maxx]-plane_data[rad];
		}
		//sU[i] = (plane_data[1770]+plane_data[1769]+plane_data[1771])/3.0-(plane_data[570]+plane_data[1750]+plane_data[1790]+plane_data[3030])/4.0;
		//sU[i] = plane_data[1770]-(plane_data[570]+plane_data[1750]+plane_data[1790]+plane_data[3030])/4.0;
		readfits_plane(Vcube,plane_data,&Vhpar);
		if(rad == -1)
		{
			sV[i] = plane_data[maxy*Ihpar.naxis[0]+maxx];
		}
		else
		{
                        /*count = 0;
                        back = 0.0;
                        for(j=0;j<plane_size;j++)
                        {
                                if((j%Vhpar.naxis[0] < (maxx-rad) || j%Vhpar.naxis[0] > (maxx+rad))\
                                ||(j/Vhpar.naxis[0] < (maxy-rad) || j/Vhpar.naxis[0] > (maxy+rad)) )
                                {
                                        back += plane_data[i];
                                        count++;
                                }
                        }
                        back/=count;*/

			sV[i] = plane_data[maxy*Ihpar.naxis[0]+maxx]-plane_data[rad];
			//sV[i] = plane_data[maxy*Vhpar.naxis[0]+maxx]-back;
		}
		//sV[i] = (plane_data[1770]+plane_data[1769]+plane_data[1771])/3.0-(plane_data[570]+plane_data[1750]+plane_data[1790]+plane_data[3030])/4.0;
		//sV[i] = plane_data[1770]-(plane_data[570]+plane_data[1750]+plane_data[1790]+plane_data[3030])/4.0;
	}

        float *x; 
	x = (float*)malloc((num_planes)*sizeof(float));

        for(n=0; n<num_planes; n++) x[n] = sI[n];
        diffusion_filter(x, num_planes, 1000);
        for(n=0; n<num_planes; n++) sI[n] = x[n];

        for(n=0; n<num_planes; n++) x[n] = sU[n];
        diffusion_filter(x, num_planes, 1000);
        for(n=0; n<num_planes; n++) sU[n] = x[n];

        for(n=0; n<num_planes; n++) x[n] = sV[n];
        diffusion_filter(x, num_planes, 1000);
        for(n=0; n<num_planes; n++) sV[n] = x[n];

	for(i=0;i<num_planes;i++)
	{
		//phi[i] = atan2(sV[i],sU[i]);
		//if(fabs(sU[i]/(2*cos(phi[i])*sI[i])) < fabs(sV[i]/(2*sin(phi[i])*sI[i])))
		//	epsilon[i] = sU[i]/(2*cos(phi[i])*sI[i]);
		//else
		//	epsilon[i] = sV[i]/(2*sin(phi[i])*sI[i]);
		epsilon[i] = sU[i]/(sI[i]);
		phi[i] = sV[i]/(sI[i]);
	}
	printf("---------------------New 4------------------\n");

        for(n=0; n<num_planes; n++) x[n] = epsilon[n];
        diffusion_filter(x, num_planes, 1000);
        for(n=0; n<num_planes; n++) epsilon[n] = x[n];

        for(n=0; n<num_planes; n++) x[n] = phi[n];
        diffusion_filter(x, num_planes, 1000);
        for(n=0; n<num_planes; n++) phi[n] = x[n];

	//for(i=0;i<num_planes;i++)
	for(i=0;i<MAX_CHANNELS;i++)
	{
		if(i < lowchan || i >= highchan)
			fprintf(outfile,"%f %f\n",0.0,0.0);
		else
		{
			//if(epsilon[i] < 0.0 || epsilon[i] > 1.0 || phi[i] < -M_PI || phi[i] > M_PI)
			if(epsilon[i] < -1.0 || epsilon[i] > 1.0 || phi[i] < -1.0 || phi[i] > 1.0)
				fprintf(outfile,"%f %f\n",0.0,0.0);
			else
				fprintf(outfile,"%f %f\n",epsilon[i],phi[i]);
		}
		//fprintf(outfile,"%f %f %f %f %f\n",epsilon[i],phi[i],sI[i],sU[i],sV[i]);
	}
	fclose(Icube);
	fclose(Ucube);
	fclose(Vcube);
	fclose(outfile);

	return EXIT_SUCCESS;
}


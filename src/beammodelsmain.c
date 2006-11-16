#include "beammodels.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include "programs/fitsio.h"
#include "fluxdata.h"
#include "jsd/jsd_futil.h"

int main(int argc, char *argv[])
{
	char outputfilename[64],cubefilename[64];
	int i;
	float max_response[MAX_CHANNELS],min_response[MAX_CHANNELS],avgQ[MAX_CHANNELS],avgU[MAX_CHANNELS],avgV[MAX_CHANNELS];
	float RA[MAX_CHANNELS],DEC[MAX_CHANNELS],radius;
	int startchannel,endchannel,beamno;
	FILE * beammodelfile;
	char subdirname[15];
	if(argc != 5)
	{
		printf("Usage: beammodels <radius> <startchannel> <endchannel> <beamno>\n");
		return EXIT_FAILURE;
	}
	else
	{
		radius = (float)atof(argv[1]);
		startchannel = atoi(argv[2]);
		endchannel = atoi(argv[3]);
		beamno = atoi(argv[4]);
	}
//	sprintf(avgfilename, "beam%d_Iavg.fits",beamno);
	sprintf(cubefilename, "beam%d_Icube.fits",beamno);
	get_peak_power_coord(cubefilename,RA,DEC,max_response,min_response,startchannel,endchannel);
	sprintf(cubefilename, "beam%d_Qcube.fits",beamno);
	get_avg(cubefilename,avgQ,startchannel,endchannel);
	sprintf(cubefilename, "beam%d_Ucube.fits",beamno);
	get_avg(cubefilename,avgU,startchannel,endchannel);
	sprintf(cubefilename, "beam%d_Vcube.fits",beamno);
	get_avg(cubefilename,avgV,startchannel,endchannel);
//	get_maxmin_power_response(cubefilename,max_response,min_response,startchannel,endchannel);

	sprintf(subdirname,"beam%d_model",beamno);
	mode_t mode = S_IRWXU|S_IRWXG|S_IRWXO;
	mkdir(subdirname,mode);
	for(i = startchannel;i < endchannel;i++)
	{
//		chdir(subdirname);
		sprintf(outputfilename, "%s/beam%d_model%03i.dat",subdirname,beamno,i);
		beammodelfile = fopen(outputfilename, "w");
		if(beammodelfile == NULL)
		{	
			printf("Error: Cannot open beam model file.\n");
			exit(EXIT_FAILURE);
		}		
//		chdir("..");
		make_beam_model(beammodelfile,i,beamno,RA[i],DEC[i],radius,max_response[i],min_response[i],avgQ[i],avgU[i],avgV[i]);
		fclose(beammodelfile);
		double maxI = 0.0;
		int linecount,j;
		FluxRecord *rec;
		beammodelfile = fopen(outputfilename, "r");
		if(beammodelfile == NULL)
		{	
			printf("Error: Cannot open beam model file.\n");
			exit(EXIT_FAILURE);
		}
		linecount = jsd_line_count(beammodelfile);
		rec = (FluxRecord *) malloc (sizeof (FluxRecord) * (linecount));
		for(j = 0;j < linecount;j++)
		{
			fscanf(beammodelfile,"%f %f %f %lf %lf %lf %lf",&rec[j].RA,&rec[j].DEC,&rec[j].AST,&rec[j].stokes.I,&rec[j].stokes.Q,\
			&rec[j].stokes.U, &rec[j].stokes.V);
			if(rec[j].stokes.I > maxI)
				maxI = rec[j].stokes.I;
		}
		fclose(beammodelfile);
		beammodelfile = fopen(outputfilename, "w");
		if(beammodelfile == NULL)
		{	
			printf("Error: Cannot open beam model file.\n");
			exit(EXIT_FAILURE);
		}		
		fprintf(beammodelfile,"dRA dDEC AST I Q U V\n");
		for(j = 0;j < linecount;j++)
		{
			fprintf(beammodelfile,"%f %f %f %lf %lf %lf %lf\n",rec[j].RA,rec[j].DEC,rec[j].AST,rec[j].stokes.I/maxI,rec[j].stokes.Q,\
			rec[j].stokes.U, rec[j].stokes.V);
		}
		fclose(beammodelfile);
	}
	printf("Beam model extracted successfully.\n");
	return EXIT_SUCCESS;
}

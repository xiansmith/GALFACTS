#include "beammodels.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include "programs/fitsio.h"
#include "fluxdata.h"
#include "jsd/jsd_futil.h"
//#include "chebyshev.h"

int multibeam;
int main(int argc, char *argv[])
{
	char outputfilename[64],cubefilename[64];
	int i;
	float max_response[MAX_CHANNELS],min_response[MAX_CHANNELS],avgQ[MAX_CHANNELS],avgU[MAX_CHANNELS],avgV[MAX_CHANNELS];
	float avgI[MAX_CHANNELS];
	float RApp[MAX_CHANNELS],DECpp[MAX_CHANNELS],radius;
	float RA,DEC,I;
	int startchannel,endchannel,beamno;
	FILE * beammodelfile;
	char subdirname[15];
	char * prefix;
	char *ccfilename;
	//if(argc != 6)
	if(argc != 6)
	{
		//printf("Usage: beammodels <radius> <startchannel> <endchannel> <beamno> <spectro_source name> <RA> <DEC>\n");
		printf("Usage: beammodels <radius> <startchannel> <endchannel> <beamno> <ccfile> \n");
		return EXIT_FAILURE;
	}
	else
	{
		radius = (float)atof(argv[1]);
		startchannel = atoi(argv[2]);
		endchannel = atoi(argv[3]);
		beamno = atoi(argv[4]);
		ccfilename = argv[5];
//		prefix = argv[5];
//		RA = (float)atof(argv[6]);
//		DEC = (float)atof(argv[7]);
	}

	multibeam = 0;
	//FluxWappData wappdata[7];
	FluxWappData * wappdata;
	char wapp[5];
	char **files;
/*	int numDays = get_date_dirs("./",&files);

        numDays = get_date_dirs("./", &files);
        if (numDays <= 0) {
                printf("ERROR: could not find any date dirs\n");
                return EXIT_FAILURE;
        }
*/
//	for (i = 0;i < 7;i++)
//	{
		//sprintf(wapp,"beam%d",beamno);
		//wappdata[i] = fluxwappdata_alloc(wapp,files,numDays);
		//wappdata = fluxwappdata_alloc(wapp,files,numDays);
//	}

//	sprintf(avgfilename, "beam%d_Iavg.fits",beamno);
	//sprintf(cubefilename, "%s_BEAM%d_%04i_%04i_Icube.fits",prefix,beamno,startchannel,endchannel-1);


/*	sprintf(cubefilename, "%s_BEAM%d_average_image_I.fits",prefix,beamno,startchannel,endchannel-1);
	printf("Cube:%s\n",cubefilename);
	get_peak_power_coord(cubefilename,RApp,DECpp,max_response,min_response,startchannel,endchannel);
*/
	FILE *ccfile = fopen(ccfilename,"r");
	if(!ccfile)
	{
		printf("Can't open %s\n",ccfilename);
	}
	else
	{
		fscanf(ccfile,"%f %f %f",&RA,&DEC,&I);
		fclose(ccfile);
	}
	
/*	get_avg(cubefilename,avgI,startchannel,endchannel);
	sprintf(cubefilename, "%s_BEAM%d_%04i_%04i_Qcube.fits",prefix,beamno,startchannel,endchannel-1);
	printf("Cube:%s\n",cubefilename);
	get_avg(cubefilename,avgQ,startchannel,endchannel);
	sprintf(cubefilename, "%s_BEAM%d_%04i_%04i_Ucube.fits",prefix,beamno,startchannel,endchannel-1);
	printf("Cube:%s\n",cubefilename);
	get_avg(cubefilename,avgU,startchannel,endchannel);
	sprintf(cubefilename, "%s_BEAM%d_%04i_%04i_Vcube.fits",prefix,beamno,startchannel,endchannel-1);
	printf("Cube:%s\n",cubefilename);
	get_avg(cubefilename,avgV,startchannel,endchannel);
*/
//	get_maxmin_power_response(cubefilename,max_response,min_response,startchannel,endchannel);

	//sprintf(subdirname,"beam%d_model",beamno);
	sprintf(subdirname,"beam_models");
	mode_t mode = S_IRWXU|S_IRWXG|S_IRWXO;
	mkdir(subdirname,mode);
	for(i = startchannel;i < endchannel;i++)
	{
//		chdir(subdirname);
//		printf("Avg:%f\n",avgI[i]);
		sprintf(outputfilename, "%s/beam%d_model%04i.dat",subdirname,beamno,i);
		//sprintf(outputfilename, "%s/beam%d_modelxxyy%04i.dat",subdirname,beamno,i);
	        //printf("output:%s\n",outputfilename);
//		beammodelfile = fopen(outputfilename, "w");
//		if(beammodelfile == NULL)
//		{	
//			printf("Error: Cannot open beam model file.\n");
//			exit(EXIT_FAILURE);
//		}		


//		chdir("..");

		//make_beam_model(beammodelfile,i,beamno,RA[i],DEC[i],radius,max_response[i],min_response[i],avgQ[i],avgU[i],avgV[i]);
		//make_beam_model(beammodelfile,i,beamno,RA[i],DEC[i],radius,max_response[i],avgI[i],avgQ[i],avgU[i],avgV[i]);
		//printf("RA %f DEC %f\n",RA[i],DEC[i]);
		//fluxwappdata_readchan(wappdata,i,CLEAN);
		make_beam_model(i,beamno,RA,DEC,radius);
		//make_beam_model(i,beamno,RApp[i],DECpp[i],radius);
//		fclose(beammodelfile);


		double maxI = 0.0;
                float RApeak=0.0, DECpeak =0.0;
		double maxQ = 0.0;
		int linecount,j;
		FluxRecord *rec;
		beammodelfile = fopen(outputfilename, "r");
		if(beammodelfile == NULL)
		{	
			printf("Error: Cannot open beam model file in 1.\n");
			exit(EXIT_FAILURE);
		}
		linecount = jsd_line_count(beammodelfile);
		float *V = (float*)malloc(sizeof(float)*linecount);
		rec = (FluxRecord *) malloc (sizeof (FluxRecord) * (linecount));
		for(j = 0;j < linecount;j++)
		{
			fscanf(beammodelfile,"%f %f %f %lf %lf %lf %lf",&rec[j].RA,&rec[j].DEC,&rec[j].AST,&rec[j].stokes.I,&rec[j].stokes.Q,\
			&rec[j].stokes.U, &rec[j].stokes.V);
			V[j] = rec[j].stokes.V;
			if(rec[j].stokes.I > maxI)
			{
                                RApeak = rec[j].RA;
                                DECpeak = rec[j].DEC;
				maxI = rec[j].stokes.I;
			}
			if(rec[j].stokes.Q > maxQ)
				maxQ = rec[j].stokes.Q;
		}
		//printf("MaxI %2.6f\n",maxI);
		fclose(beammodelfile);

		//float min,max;
                //chebyshev_minmax(V, linecount, &min, &max);
                //chebyshev_normalize(V, linecount, min, max);		

                float dist = 60.0*sqrt(RApeak*RApeak+DECpeak*DECpeak);
                maxI /= exp(-2.772589*(dist/3.4)*(dist/3.4));
		printf("Distance %f, Ratio %f\n",dist,exp(-2.772589*(dist/3.4)*(dist/3.4)));

		beammodelfile = fopen(outputfilename, "w");
		if(beammodelfile == NULL)
		{	
			printf("Error: Cannot open beam model file in 2.\n");
			exit(EXIT_FAILURE);
		}		
		//fprintf(beammodelfile,"dRA dDEC AST I Q U V\n");


		for(j = 0;j < linecount;j++)
		{
			fprintf(beammodelfile,"%f %f %f %lf %lf %lf %lf\n",rec[j].RA,rec[j].DEC,rec[j].AST,rec[j].stokes.I/maxI,rec[j].stokes.Q/maxI,\
			rec[j].stokes.U/maxI, rec[j].stokes.V/maxI);
			//fprintf(beammodelfile,"%f %f %f %lf %lf %lf %lf\n",rec[j].RA,rec[j].DEC,rec[j].AST,rec[j].stokes.I,rec[j].stokes.Q,\
			rec[j].stokes.U, rec[j].stokes.V);
			//fprintf(beammodelfile,"%f %f %f %lf %lf %lf %lf\n",rec[j].RA,rec[j].DEC,rec[j].AST,rec[j].stokes.I*0.5/maxI,rec[j].stokes.Q*0.5/maxQ,\
			rec[j].stokes.U/maxI, rec[j].stokes.V/maxI);
		}
		fclose(beammodelfile);
	}
//	printf("Beam model extracted successfully.\n");
	return EXIT_SUCCESS;
}

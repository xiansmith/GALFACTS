#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include "programs/fitsio.h"
//#include "fluxdata.h"
//#include "jsd/jsd_futil.h"
/*#include "beammodels.h"

void get_avg(char * filename,float avg[MAX_CHANNELS],int startchannel,int endchannel)
{
        FILE * beam_file;
        header_param_list beam_hpar;
        int data_len;
        int i,j;
//      float temp;
        int count;
        double sum;
        float * plane_data;
        beam_file = fopen(filename, "r");
        readfits_header(beam_file, &beam_hpar);
        //allocate memory for plane
        data_len = beam_hpar.naxis[0] * beam_hpar.naxis[1];
        plane_data = (float*) calloc(data_len, sizeof (float));
        //read the plane
        for (i=0; i< (endchannel-startchannel); i++)
        {
                //read a plane
                readfits_plane(beam_file, plane_data, &beam_hpar);
                sum = 0.0;
                count = 0;
                for (j=0; j<data_len; j++)
                {
                        if (!IS_BLANK_PIXEL(plane_data[j]) && isfinite(plane_data[j]))
                        {
                                sum += plane_data[j];
                                count++;
                        }
                }
                avg[startchannel + i] = sum/count;
//              printf("Calculated RA: %f, DEC: %f\n",RA[startchannel + i],DEC[startchannel + i]);
//              printf("Max response: %f, Min response: %f\n",max_response[startchannel+i],min_response[startchannel + i]);
        }
        free(plane_data);
        fclose(beam_file);
}
*/
int get_peak_power_coord(char * filename)
{
        FILE * beam_file;
        header_param_list beam_hpar;
        int data_len;
        int j;
//      float temp;
        int indx;
        float * plane_data;
        beam_file = fopen(filename, "r");
        readfits_header(beam_file, &beam_hpar);
        //allocate memory for plane
        data_len = beam_hpar.naxis[0] * beam_hpar.naxis[1];
        plane_data = (float*) calloc(data_len, sizeof (float));
        //read the plane
//        for (i=0; i< (endchannel-startchannel); i++)
//        {
                //read a plane
                float max_response = 0;
  //              min_response[startchannel+i] = 100;
                indx = 0;
                readfits_plane(beam_file, plane_data, &beam_hpar);
                for (j=0; j<data_len; j++)
                {
                        if (!IS_BLANK_PIXEL(plane_data[j]) && isfinite(plane_data[j]))
                        {
                                if(max_response < plane_data[j])
                                {
                                        max_response = plane_data[j];
                                        indx = j;
                                }
//                                if(min_response[startchannel+i] > plane_data[j])
//                                        min_response[startchannel+i] = plane_data[j];
                        }
                }
//                l = (int)(indx % beam_hpar.naxis[0]);
//                m = (int)((indx - l)/ beam_hpar.naxis[0]);
//                RA[startchannel + i] = (l - beam_hpar.crpix[0] + 1)*beam_hpar.cdelt[0] + beam_hpar.crval[0];
//                DEC[startchannel + i] = (m - beam_hpar.crpix[1] + 1)*beam_hpar.cdelt[1] + beam_hpar.crval[1];
//                printf("Calculated RA: %f, DEC: %f\n",RA[startchannel + i],DEC[startchannel + i]);
//                printf("Max response: %f, Min response: %f\n",max_response[startchannel+i],min_response[startchannel + i]);
//        }
        free(plane_data);
        fclose(beam_file);
	return indx;
}

void main(int argc, char* argv[])
{
	FILE *modelI,*modelQ,*modelU;
	FILE *dataI,*dataQ,*dataU;
        header_param_list model_hpar,data_hpar;
	char * subdirname;
	int beam,modelpeak,datapeak;
        if(argc != 3)
        {
                //printf("Usage: beammodels <radius> <startchannel> <endchannel> <beamno> <spectro_source name>\n");
                printf("Usage: rmtest <subdir> <beam>\n");
                exit(1);
        }
	else
	{
		subdirname = argv[1];
		beam = atoi(argv[2]);
	}

	float maxImodel[1800],maxIdata[1800];
	float minImodel[1800],minIdata[1800];
	float ratioI[1800],totI[1800],bandshapeI[1800],I2[1800];
	int i;

	char modelfile[64],datafile[64];
	//sprintf(modelfile,"S0206+330/S0206+330_BEAM%d_0150_1949_Icube.fits",beam);
	sprintf(modelfile,"S0226+343/S0226+343_BEAM%d_0150_1949_Icube.fits",beam);
	sprintf(datafile,"%s/%s_BEAM%d_0150_1949_Icube.fits",subdirname,subdirname,beam);
	//sprintf(datafile,"full/CORR_TEST_0150_1949_Icube.fits");
	modelpeak = get_peak_power_coord(modelfile);
	datapeak = get_peak_power_coord(datafile);
	printf("Data peak: %d Model peak %d\n",datapeak,modelpeak);
	modelI = fopen(modelfile,"r");
	if(modelI == NULL)
	{
		printf("Cannot open model file\n");
		exit(1);
	}
	readfits_header(modelI,&model_hpar);
	dataI = fopen(datafile,"r");
	if(dataI == NULL)
	{
		printf("Cannot open data file\n");
		exit(1);
	}
	readfits_header(dataI,&data_hpar);
	int modeldatalen = model_hpar.naxis[0] * model_hpar.naxis[1];
	int datalen = data_hpar.naxis[0] * data_hpar.naxis[1];
	//printf("datalen = %d\n",datalen);
	//printf("start = %f\n",data_hpar.crval[2]);
	//printf("delta = %f\n",data_hpar.cdelt[2]);
	//exit(1);
	float *modelplane = (float*) calloc(modeldatalen, sizeof (float));
	float *dataplane = (float*) calloc(datalen, sizeof (float));
	float avgmodelI = 0,avgdataI = 0;
	float avgmodelQ = 0,avgdataQ = 0;
	float avgmodelU = 0,avgdataU = 0;

	FILE * file;
	file = fopen("totI.dat","w");

	for(i = 0;i < 1800;i++)
	{
//		printf("i = %d\n",i);
                readfits_plane(modelI, modelplane, &model_hpar);
                readfits_plane(dataI, dataplane, &data_hpar);
		//maxImodel[i] = modelplane[1951];
		maxImodel[i] = modelplane[modelpeak];
		maxIdata[i] = dataplane[datapeak];
		//minImodel[i] = modelplane[331];
		minImodel[i] = modelplane[329];
		minIdata[i] = dataplane[329];
		avgmodelI += (maxImodel[i]-minImodel[i]);
		avgdataI += (maxIdata[i]-minIdata[i]);
	}
	float avgratioI;
	avgmodelI/=1800;
	avgdataI/=1800;
	printf("%f\n",avgdataI/avgmodelI);
	avgratioI = avgdataI/avgmodelI;
	//avgratioI = 0.0;///////////////////CAreful !!!!!!
	//exit(1);
	for(i = 0;i < 1800;i++)
	{
		//ratioI[i] = ((maxIdata[i]-minIdata[i])*maxImodel[i])/((maxImodel[i]-minImodel[i])*maxIdata[i]);
		bandshapeI[i] = (maxImodel[i]-minImodel[i])/avgmodelI;
		ratioI[i] = (maxIdata[i]-minIdata[i])/(maxImodel[i]-minImodel[i]);
		//ratioI[i] = (maxImodel[i]-minImodel[i])/avgmodelI;
		//ratioI[i] = (maxIdata[i]-minIdata[i])/avgmodelI;
		//printf("%f %f %f\n",(maxImodel[i]-minImodel[i])/maxImodel[i],(maxIdata[i]-minIdata[i])/maxIdata[i],ratioI[i]);
		//printf("%d %f\n",i,ratioI[i]);
		//printf("%f\n",maxIdata[i]-minIdata[i]);
		//totI[i] = maxIdata[i]-minIdata[i];
		totI[i] = (maxIdata[i]-minIdata[i])/bandshapeI[i];
		I2[i] = (maxIdata[i]-minIdata[i])*(maxIdata[i]-minIdata[i]);
		//fprintf(file,"%f\n",totI[i]);
		fprintf(file,"%d %f %f %f %f\n",i,ratioI[i],totI[i],(maxImodel[i]-minImodel[i]),bandshapeI[i]);
	}
	fclose(file);
	fclose(modelI);
	fclose(dataI);


	//sprintf(modelfile,"S0206+330/S0206+330_BEAM%d_0150_1949_Qcube.fits",beam);
	sprintf(modelfile,"S0226+343/S0226+343_BEAM%d_0150_1949_Qcube.fits",beam);
	sprintf(datafile,"%s/%s_BEAM%d_0150_1949_Qcube.fits",subdirname,subdirname,beam);
	//sprintf(datafile,"full/CORR_TEST_0150_1949_Qcube.fits");
	modelQ = fopen(modelfile,"r");
	if(modelQ == NULL)
	{
		printf("Cannot open model file\n");
		exit(1);
	}
	readfits_header(modelQ,&model_hpar);
	dataQ = fopen(datafile,"r");
	if(dataQ == NULL)
	{
		printf("Cannot open data file\n");
		exit(1);
	}
	readfits_header(dataQ,&data_hpar);
	int j;
	float resQ[1800];
	float maxQmodel[1800],maxQdata[1800];
	float minQmodel[1800],minQdata[1800];
	float ratioQ[1800],totQ[1800],bandshapeQ[1800],Q2[1800];
	float avgratioQ,ratioQavg=0;
	
	file = fopen("resQ.dat","w");
	for(i = 0;i < 1800;i++)
	{
                readfits_plane(modelQ, modelplane, &model_hpar);
                readfits_plane(dataQ, dataplane, &data_hpar);
		//maxQmodel[i] = modelplane[1951];
		maxQmodel[i] = modelplane[modelpeak];
		maxQdata[i] = dataplane[datapeak];
		//minQmodel[i] = modelplane[331];
		minQmodel[i] = modelplane[329];
		minQdata[i] = dataplane[329];
		totQ[i] = maxQdata[i]-minQdata[i];
		Q2[i] = totQ[i]*totQ[i];
		avgmodelQ += (maxQmodel[i]-minQmodel[i]);
		avgdataQ += (maxQdata[i]-minQdata[i]);
		ratioQ[i] = (maxQdata[i]-minQdata[i])/(maxQmodel[i]-minQmodel[i]);
		ratioQavg+=ratioQ[i];
		//if(ratioQ[i] > 1.0 || ratioQ[i] < -1.0)
		//	ratioQ[i] = 0.0;
		//resQ[i] = maxQdata[i] - (0.310778+0.0000115146*i)*maxQmodel[i];
		//resQ[i] = maxQdata[i] - (0.489273+0.0000220088*i)*maxQmodel[i];
		//resQ[i] = maxQdata[i]-minQdata[i];
		//resQ[i] = maxQdata[i] -minQdata[i] - (0.489273+0.0000220088*i)*(maxQmodel[i]-minQmodel[i]); //S0340/S0226
		//resQ[i] = maxQdata[i] -minQdata[i] - (0.507826)*(maxQmodel[i]-minQmodel[i]); //S0340/S0226
		//resQ[i] = maxQdata[i] -minQdata[i] - (0.464256+0.0000213166*i)*(maxQmodel[i]-minQmodel[i]); //S0340/S0226 b 1
		//resQ[i] = maxQdata[i] -minQdata[i] - (0.184031+0.00000840703*i)*(maxQmodel[i]-minQmodel[i]); //S0352/S0226
		//resQ[i] = maxQdata[i] -minQdata[i] - (0.310778+0.0000115146*i)*(maxQmodel[i]-minQmodel[i]); //S0217/S0226
		//resQ[i] = maxQdata[i] -minQdata[i] - (2.01364+0.00000123033*i)*(maxQmodel[i]-minQmodel[i]); //S0340/S0206
		//resQ[i] = maxQdata[i] -minQdata[i] - (0.757459+0.000000887724*i)*(maxQmodel[i]-minQmodel[i]); //S0352/S0206
		//resQ[i] = maxQdata[i];
		//resQ[i] = maxQmodel[i]-minQmodel[i];
		//printf("%f %f\n",(maxQmodel[i]-minQmodel[i])/maxQmodel[i],(maxQdata[i]-minQdata[i])/maxQdata[i]);
		//printf("%f %f\n",(maxQmodel[i]-minQmodel[i]),(maxQdata[i]-minQdata[i]));
		//printf("%f %f\n",(maxQmodel[i]),(maxQdata[i]));
		//printf("%d %f\n",i,resQ[i]);
		//fprintf(file,"%d %f %f\n",i,resQ[i]/totI[i],resQ[i]);
		//fprintf(file,"%d %f\n",i,ratioQ[i]);
	}
	fclose(modelQ);
	fclose(dataQ);
	avgmodelQ/=1800;
	avgdataQ/=1800;
	ratioQavg/=1800;
	avgratioQ = avgdataQ/avgmodelQ;
	printf("%f %f\n",avgdataQ/avgmodelQ,ratioQavg);

	for(i = 0;i < 1800;i++)
	{
		bandshapeQ[i] = (maxQmodel[i]-minQmodel[i])/avgmodelQ;
		resQ[i] = maxQdata[i] - minQdata[i] - avgratioI*(maxQmodel[i]-minQmodel[i]); //S0340/S0226
		//fprintf(file,"%d %f %f %f %f %f %f\n",i,resQ[i]/totI[i],resQ[i]/bandshapeQ[i],resQ[i],bandshapeQ[i],maxQmodel[i]-minQmodel[i],Q2[i]);
		//fprintf(file,"%d %f %f %f %f %f %f\n",i,resQ[i]/totI[i],resQ[i]/bandshapeQ[i],resQ[i],bandshapeQ[i],totQ[i],ratioQ[i]);
		fprintf(file,"%d %f %f %f\n",i,maxQdata[i]-minQdata[i],maxQmodel[i]-minQmodel[i],(maxQmodel[i]-minQmodel[i])*avgratioQ);
	}
	fclose(file);
	//sprintf(modelfile,"S0206+330/S0206+330_BEAM%d_0150_1949_Ucube.fits",beam);
	sprintf(modelfile,"S0226+343/S0226+343_BEAM%d_0150_1949_Ucube.fits",beam);
	sprintf(datafile,"%s/%s_BEAM%d_0150_1949_Ucube.fits",subdirname,subdirname,beam);
	//sprintf(datafile,"full/CORR_TEST_0150_1949_Ucube.fits");
	modelU = fopen(modelfile,"r");
	if(modelU == NULL)
	{
		printf("Cannot open model file\n");
		exit(1);
	}
	readfits_header(modelU,&model_hpar);
	dataU = fopen(datafile,"r");
	if(dataU == NULL)
	{
		printf("Cannot open data file\n");
		exit(1);
	}
	readfits_header(dataU,&data_hpar);
	float resU[1800];
	float maxUmodel[1800],maxUdata[1800];
	float minUmodel[1800],minUdata[1800];
	float ratioU[1800],bandshapeU[1800],totU[1800],U2[1800];
	float avgratioU,ratioUavg=0;
	file = fopen("resU.dat","w");
	
	for(i = 0;i < 1800;i++)
	{
                readfits_plane(modelU, modelplane, &model_hpar);
                readfits_plane(dataU, dataplane, &data_hpar);
		//maxUmodel[i] = modelplane[1951];
		maxUmodel[i] = modelplane[modelpeak];
		maxUdata[i] = dataplane[datapeak];
		//minUmodel[i] = modelplane[331];
		minUmodel[i] = modelplane[329];
		minUdata[i] = dataplane[329];
		avgmodelU += (maxUmodel[i]-minUmodel[i]);
		avgdataU += (maxUdata[i]-minUdata[i]);
		totU[i] = maxUdata[i]-minUdata[i];
		U2[i] = totU[i]*totU[i];
		ratioU[i] = (maxUdata[i]-minUdata[i])/(maxUmodel[i]-minUmodel[i]);
		ratioUavg+=ratioU[i];
		//if(ratioQ[i] > 1.0 || ratioQ[i] < -1.0)
		//	ratioQ[i] = 0.0;
		//resU[i] = maxUdata[i] - (0.310778+0.0000115146*i)*maxUmodel[i];
		//resU[i] = maxUdata[i] - (0.489723+0.0000220088*i)*maxUmodel[i];
		//resU[i] = maxUdata[i]-minUdata[i];
		//resU[i] = maxUdata[i] - minUdata[i] - (0.489273+0.0000220088*i)*(maxUmodel[i]-minUmodel[i]); //S0340/S0226
		//resU[i] = maxUdata[i] -minUdata[i] - (0.507826)*(maxUmodel[i]-minUmodel[i]); //S0340/S0226
		//resU[i] = maxUdata[i] - minUdata[i] - (0.464256+0.0000213166*i)*(maxUmodel[i]-minUmodel[i]); //S0340/S0226 b1
		//resU[i] = maxUdata[i] - minUdata[i] - (0.184031+0.00000840703*i)*(maxUmodel[i]-minUmodel[i]); //S0352/S0226
		//resU[i] = maxUdata[i] -minUdata[i] - (0.310778+0.0000115146*i)*(maxUmodel[i]-minUmodel[i]); //S0217/S0226
		//resU[i] = maxUdata[i] - minUdata[i] - (0.757459+0.000000887724*i)*(maxUmodel[i]-minUmodel[i]); //S0352/S0206
		//resU[i] = maxUdata[i];
		//resU[i] = maxUmodel[i]-minUmodel[i];
		//printf("%f %f\n",(maxQmodel[i]-minQmodel[i])/maxQmodel[i],(maxQdata[i]-minQdata[i])/maxQdata[i]);
		//printf("%f %f\n",(maxQmodel[i]-minQmodel[i]),(maxQdata[i]-minQdata[i]));
		//printf("%f %f\n",(maxQmodel[i]),(maxQdata[i]));
		//printf("%d %f\n",i,resU[i]);
		//fprintf(file,"%d %f %f\n",i,resU[i]/totI[i],resU[i]);
		//fprintf(file,"%d %f\n",i,ratioU[i]);
	}
	fclose(modelU);
	fclose(dataU);
	avgmodelU/=1800;
	avgdataU/=1800;
	ratioUavg/=1800;
	avgratioU = avgdataU/avgmodelU;
	printf("%f %f\n",avgdataU/avgmodelU,ratioUavg);
	for(i = 0;i < 1800;i++)
	{
		bandshapeU[i] = (maxUmodel[i]-minUmodel[i])/avgmodelU;
		resU[i] = maxUdata[i] - minUdata[i] - avgratioI*(maxUmodel[i]-minUmodel[i]); //S0340/S0226
		//fprintf(file,"%d %f %f %f %f %f %f\n",i,resU[i]/totI[i],resU[i]/bandshapeU[i],resU[i],bandshapeU[i],maxUmodel[i]-minUmodel[i],U2[i]);
		//fprintf(file,"%d %f %f %f %f %f %f\n",i,resU[i]/totI[i],resU[i]/bandshapeU[i],resU[i],bandshapeU[i],totU[i],ratioU[i]);
		fprintf(file,"%d %f %f %f\n",i,maxUdata[i]-minUdata[i],maxUmodel[i]-minUmodel[i],(maxUmodel[i]-minUmodel[i])*avgratioU);
	}
	fclose(file);

	float pI[1800],pA[1800];
	file = fopen("pIA.dat","w");
	float lambda2,lambda,nu;
	float c = 299792458.0;
	for(i = 0;i < 1800;i++)
	{
		pI[i] = (sqrt(resQ[i]*resQ[i]+resU[i]*resU[i]))/totI[i];
		pA[i] = 0.5*atan2(resU[i],resQ[i]);
		pA[i] = 0.5*atan2(maxUdata[i]-minUdata[i],maxQdata[i]-minQdata[i]);
		if(pA[i] < 0)
			pA[i] += M_PI;
		//printf("%f \n",pI[i]);
		//nu = (1450000000.0+(data_hpar.crval[2]-1450000000.0)*2+data_hpar.cdelt[2]*i);
		nu = (data_hpar.crval[2]+data_hpar.cdelt[2]*i);
		lambda = c/nu;
		//printf("%f %f %f\n",nu,lambda,data_hpar.crval[2]);
		lambda2 = lambda*lambda;

		fprintf(file,"%d %f %f %f %f %f\n",i,pI[i],pA[i],lambda2,sqrt(Q2[i]+U2[i]),I2[i]);
	}
	fclose(file);
	
}

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include "programs/fitsio.h"
#include "fluxdata.h"
#include "jsd/jsd_futil.h"
#include "beammodels.h"
//#include "common.h"

/*void get_maxmin_power_response(char * filename,float max_response[MAX_CHANNELS],float min_response[MAX_CHANNELS],int startchannel,int endchannel)
{
	FILE * beam_file;
	header_param_list beam_hpar;
	int data_len;
	int i,j;
	float * plane_data;
	beam_file = fopen(filename, "r");
	readfits_header(beam_file, &beam_hpar);
//	printf("In function call\n");
	//allocate memory for plane
	data_len = beam_hpar.naxis[0] * beam_hpar.naxis[1];
	plane_data = (float*) calloc(data_len, sizeof (float));
	//read the plane
//	readfits_plane(beam_file, plane_data, &beam_hpar);
	for (i=0; i<(endchannel - startchannel); i++) 
	{
		//read a plane
		max_response[startchannel+i] = 0;
		min_response[startchannel+i] = 100;
		readfits_plane(beam_file, plane_data, &beam_hpar);
		for (j=0; j<data_len; j++) 
		{
			if (!IS_BLANK_PIXEL(plane_data[j]) && isfinite(plane_data[j])) 
			{
				if(max_response[startchannel+i] < plane_data[j])
					max_response[startchannel+i] = plane_data[j];
				if(min_response[startchannel+i] > plane_data[j])
					min_response[startchannel+i] = plane_data[j];
			}
		}
//		printf("Max Value : %f\n",peak_response[startchannel+i]);
	}
	free(plane_data);
	fclose(beam_file);
}*/

void get_avg(char * filename,float avg[MAX_CHANNELS],int startchannel,int endchannel)
{
	FILE * beam_file;
	header_param_list beam_hpar;
	int data_len;
	int i,j;
//	float temp;
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
//		printf("Calculated RA: %f, DEC: %f\n",RA[startchannel + i],DEC[startchannel + i]);
//		printf("Max response: %f, Min response: %f\n",max_response[startchannel+i],min_response[startchannel + i]);
	}
	free(plane_data);
	fclose(beam_file);
}

void get_peak_power_coord(char * filename,float RA[MAX_CHANNELS],float DEC[MAX_CHANNELS],float max_response[MAX_CHANNELS],float min_response[MAX_CHANNELS],int startchannel,int endchannel)
{
	FILE * beam_file;
	header_param_list beam_hpar;
	int data_len;
	int i,j;
//	float temp;
	int indx,l,m;
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
		max_response[startchannel+i] = 0;
		min_response[startchannel+i] = 100;
		indx = 0;
		readfits_plane(beam_file, plane_data, &beam_hpar);
		for (j=0; j<data_len; j++) 
		{
			if (!IS_BLANK_PIXEL(plane_data[j]) && isfinite(plane_data[j])) 
			{
				if(max_response[startchannel+i] < plane_data[j])
				{
					max_response[startchannel+i] = plane_data[j];		
					indx = j;
				}
				if(min_response[startchannel+i] > plane_data[j])
					min_response[startchannel+i] = plane_data[j];
			}
		}	
		l = (int)(indx % beam_hpar.naxis[0]);
		//m = (int)((indx - l)/ beam_hpar.naxis[0]);
		m = (int)((indx)/ beam_hpar.naxis[0]);
		RA[startchannel + i] = (l - beam_hpar.crpix[0] + 1)*beam_hpar.cdelt[0] + beam_hpar.crval[0];
		DEC[startchannel + i] = (m - beam_hpar.crpix[1] + 1)*beam_hpar.cdelt[1] + beam_hpar.crval[1];
		printf("Calculated RA: %f, DEC: %f\n",RA[startchannel + i],DEC[startchannel + i]);
		printf("Max response: %f, Min response: %f\n",max_response[startchannel+i],min_response[startchannel + i]);
	}
	free(plane_data);
	fclose(beam_file);
}

//void make_beam_model(FILE * beammodelfile,int channel,int beamno,float RA,float DEC,float radius,float max_response,float min_response,float avgQ,float avgU,float avgV)
void make_beam_model(int channel,int beamno,float RA,float DEC,float radius)
{
	int i,j,numDays,count;
	char **files;
	char infilename[64],outfilename[80];
	FILE * infile,* beammodelfile;
        char header[80+1];
	FluxRecord pRec;
	numDays = get_date_dirs("./", &files);
	float maxI = 0.0;
	float maxV = -10.0;
	float minV = 10.0;
	//float maxV = 10.0;
	float maxU = 10.0;
	float maxQ = 10.0;
/*	for(j = 0;j < numDays;j++)
	{
		sprintf(infilename,"%s/beam%d/balance%03i.dat",files[j],beamno,channel);
		infile = fopen(infilename, "r");
	        if (infile == NULL) {
        		printf("ERROR: can't open input file %s\n", infilename);
		        continue;
	        }

		count = jsd_line_count(infile);
	        fgets(header, 80, infile);

		for(i = 0;i < count;i++)
		{
			fscanf(infile,"%f %f %f %lf %lf %lf %lf",&pRec.RA,&pRec.DEC,&pRec.AST,&pRec.stokes.I,&pRec.stokes.Q,&pRec.stokes.U, &pRec.stokes.V);

			if(maxI < pRec.stokes.I && !isnan(pRec.stokes.I)) 
			{
				maxI = pRec.stokes.I;
			}
		}
	        fclose(infile);
	}
*/
//	float RA0[10000],DEC0[10000],cRA=0.0,cDEC=0.0;

//	for(beamno =  0;beamno <7; beamno++)
//	{
		int cnt = 0;
		int day = 0;
		int rec = 0;
		//float maxI,minI,avgQ,avgU,avgV;
		for(j = 0;j < numDays;j++)
		{
			sprintf(infilename,"%s/beam%d/balance%04i.dat",files[j],beamno,channel);
			//sprintf(infilename,"%s/beam%d/balanceraw%04i.dat",files[j],beamno,channel);
//			sprintf(infilename,"%s/beam%d/fluxtime%04i.dat",files[j],beamno,channel);
			infile = fopen(infilename, "r");
		        if (infile == NULL) {
        			printf("ERROR: can't open input file %s\n", infilename);
			        continue;
		        }
			//printf("Reading %s\n",infilename);
			//sprintf(outfilename,"beam%d_model/beam%d_model%04i.dat",beamno,beamno,channel);
			sprintf(outfilename,"beam_models/beam%d_model%04i.dat",beamno,channel);
			//sprintf(outfilename,"beam%d_model/beam%d_modelxxyy%04i.dat",beamno,beamno,channel);
			beammodelfile = fopen(outfilename, "a");
		        if (beammodelfile == NULL) {
        			printf("ERROR: can't open output file %s\n", outfilename);
			        continue;
		        }
	
			count = jsd_line_count(infile);
	        	fgets(header, 80, infile);

			for(i = 0;i < count;i++)
			{
				fscanf(infile,"%f %f %f %lf %lf %lf %lf",&pRec.RA,&pRec.DEC,&pRec.AST,&pRec.stokes.I,&pRec.stokes.Q,&pRec.stokes.U, &pRec.stokes.V);
				if(maxI < pRec.stokes.I && !isnan(pRec.stokes.I)) 
				{
					maxI = pRec.stokes.I;
					day = j;
					rec = i;
					//maxV = pRec.stokes.V;
					//maxQ = pRec.stokes.Q;
					//maxU = pRec.stokes.U;
				}
				//if(fabs(maxV) < fabs(pRec.stokes.V) && !isnan(pRec.stokes.V)) 
				/*if(maxQ < (pRec.stokes.Q) && !isnan(pRec.stokes.Q)) 
				{
					maxQ = pRec.stokes.Q;
					//maxV = pRec.stokes.V;
				}
				if(maxU > (pRec.stokes.U) && !isnan(pRec.stokes.U)) 
				{
					maxU = pRec.stokes.U;
					//maxV = pRec.stokes.V;
				}*/
				//if(fabs(maxV) < fabs(pRec.stokes.V) && !isnan(pRec.stokes.V)) 
				if(maxV < pRec.stokes.V && !isnan(pRec.stokes.V)) 
				{
					maxV = pRec.stokes.V;
					//maxV = pRec.stokes.V;
				}
				if(minV > pRec.stokes.V && !isnan(pRec.stokes.V)) 
				{
					minV = pRec.stokes.V;
					//maxV = pRec.stokes.V;
				}
/*				if(i)
				{
					cRA = (RA0[cnt] - RA0[cnt-1])*0.5;
					cDEC = (DEC0[cnt] - DEC0[cnt-1])*0.5;
				}
				else
				{
					cRA = 0.0;
					cDEC = 0.0;
				}

				//for MOCK pointing errors
				switch(beamno)
				{
	        	            case 0:
					RA0[cnt] = pRec.RA;
					DEC0[cnt] = pRec.DEC;
        	        	        break; 
		                    case 1:
        	                	pRec.RA = RA0[cnt] + 2.7417/(60*cos(DEC0[cnt]*M_PI/180)) + cRA;
					pRec.DEC = DEC0[cnt] - 5.5426/60 + cDEC;
        		                break; 
	                	    case 2:
        	                	pRec.RA = RA0[cnt] + 5.4833/(60*cos(DEC0[cnt]*M_PI/180)) + cRA;
					pRec.DEC = DEC0[cnt] + cDEC;
	        	                break; 
		                    case 3:
        	                	pRec.RA = RA0[cnt] + 2.7417/(60*cos(DEC0[cnt]*M_PI/180)) + cRA;
					pRec.DEC = DEC0[cnt] + 5.5426/60 + cDEC;
        		                break; 
	        	            case 4:
        	                	pRec.RA = RA0[cnt] - 2.7417/(60*cos(DEC0[cnt]*M_PI/180)) + cRA;
					pRec.DEC = DEC0[cnt] + 5.5426/60 + cDEC;
        		                break; 
	        	            case 5:
        	                	pRec.RA = RA0[cnt] + 5.4833/(60*cos(DEC0[cnt]*M_PI/180)) + cRA;
					pRec.DEC = DEC0[cnt] + cDEC;
        		                break; 
	        	            case 6:
        	                	pRec.RA = RA0[cnt] - 2.7417/(60*cos(DEC0[cnt]*M_PI/180)) + cRA;
					pRec.DEC = DEC0[cnt] - 5.5426/60 + cDEC;
        		                break; 
        		            default:
	                	        printf("WARN: invalid beam (%i) specified.\n", beamno);
				}
				cnt++;
*/

				if(pRec.DEC <= (DEC + radius) && pRec.DEC >= (DEC - radius) && pRec.RA <= (RA + radius) && pRec.RA >= (RA - radius) \
				&& !isnan(pRec.stokes.I)) 
				{
					//fprintf(beammodelfile,"%f %f %f %lf %lf %lf %lf\n",(pRec.RA - RA),(pRec.DEC - DEC),pRec.AST,(pRec.stokes.I - min_response) \
					/(max_response- min_response),	(pRec.stokes.Q - avgQ)/(max_response- min_response),  \
					(pRec.stokes.U - avgU)/(max_response- min_response), \
					(pRec.stokes.V - avgV)/(max_response- min_response));
					//fprintf(beammodelfile,"%f %f %f %lf %lf %lf %lf\n",(pRec.RA - RA),(pRec.DEC - DEC),\
					pRec.AST, pRec.stokes.I,(pRec.stokes.Q - avgQ)/(max_response- min_response),  \
					(pRec.stokes.U - avgU)/(max_response- min_response), \
					(pRec.stokes.V - avgV)/(max_response- min_response));
					fprintf(beammodelfile,"%f %f %f %lf %lf %lf %lf\n",(pRec.RA - RA),(pRec.DEC - DEC),pRec.AST, \
					pRec.stokes.I,pRec.stokes.Q,pRec.stokes.U,pRec.stokes.V);
				}
			}
		        fclose(infile);
		        fclose(beammodelfile);
		}
		//printf("MaxI : %f\n",maxI);
		
		//printf("Beam %d V/I at maxI : %f\n",beamno,maxV/maxI);
		//printf("Beam %d Q/I at maxI : %f\n",beamno,maxQ/maxI);
		//printf("Beam %d U/I at maxI : %f\n",beamno,maxU/maxI);
		printf("Beam %d +V/I: %2.6f -V/I: %2.6f MaxI %2.6f Day %d Rec %d\n",beamno,maxV/maxI,minV/maxI,maxI,day,rec);
		//printf("Beam %d maxV %f /maxI %f : %f\n",beamno,maxV,maxI,maxV/maxI);
//	}
}

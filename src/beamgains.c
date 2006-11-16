#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "common.h"
#include "programs/fitsio.h"
#include "beamgains.h"

void get_peak_power_response(char * filename,float peak_response[MAX_CHANNELS],int startchannel,int endchannel)
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
	for (i=0; i<beam_hpar.naxis[2]; i++) 
	{
		//read a plane
		peak_response[startchannel+i] = 0;
		readfits_plane(beam_file, plane_data, &beam_hpar);
		for (j=0; j<data_len; j++) 
		{
			if (!IS_BLANK_PIXEL(plane_data[j]) && isfinite(plane_data[j])) 
			{
				if(peak_response[startchannel+i] < plane_data[j])
					peak_response[startchannel+i] = plane_data[j];
			}
		}
		printf("Max Value : %f\n",peak_response[startchannel+i]);
	}
	free(plane_data);
	fclose(beam_file);
}

int main(int argc, char *argv[])
{
	char filename[64];
	char outputfilename[64];
	float peak_response[7][MAX_CHANNELS];
	float temp[MAX_CHANNELS];
	float temp_gain;
	float sums[NUM_BEAMS];
	FILE * beamgainfile,* testfile;
	int i,j;
	int startchannel,endchannel;
//	char * pathname;
	if(argc != 4)
	{
		printf("Usage: beamgains <startchannel> <endchannel> <region>\n");
		return EXIT_FAILURE;
	}
	else
	{
		startchannel = atoi(argv[1]);
		endchannel = atoi(argv[2]);
	}
//	chdir(pathname);
	sprintf(outputfilename, "beamgains%s.dat",argv[3]);
	beamgainfile = fopen(outputfilename, "w");
	testfile = fopen("test.dat", "w");
	fprintf(beamgainfile, "#CHAN BEAM0   BEAM1      BEAM2      BEAM3      BEAM4      BEAM5      BEAM6\n");
	printf("Calculating beam gains...\n");
	for(i = 0;i < NUM_BEAMS;i++)
	{
//		printf("In loop 1\n");
	        sprintf(filename, "beam%d_Icube.fits",i);
		get_peak_power_response(filename,temp,startchannel,endchannel);
		fprintf(testfile,"beam%d\n",i);
		for(j = startchannel;j < endchannel;j++)
		{
//			printf("In nested loop\n");
			peak_response[i][j] = temp[j];
			fprintf(testfile,"%d %1.8f\n",j,temp[j]);
//			printf("%f\n",temp[j]);
		}
		sums[i] = 0;
	}
	for(i = startchannel;i < endchannel;i++)
	{
//		printf("In loop 2\n");
		fprintf(beamgainfile,"%d",i);
		for(j = 0;j < NUM_BEAMS;j++)
		{
			temp_gain = peak_response[j][i]/peak_response[0][i];
			fprintf(beamgainfile," %1.8f",temp_gain);
			sums[j] += temp_gain;
		}
		fprintf(beamgainfile,"\n");

	}	
	j = 0;
	fprintf(beamgainfile,"%d",j);//avgs
	for(j = 0;j < NUM_BEAMS;j++)
	{
			sums[j] = sums[j]/(endchannel-startchannel);
			fprintf(beamgainfile," %1.8f",sums[j]);
	}

	printf("Done !\n");
	fclose(beamgainfile);
	fclose(testfile);
	return EXIT_SUCCESS;
}

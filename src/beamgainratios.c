#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "common.h"
#include "programs/fitsio.h"
#include "beamgains.h"


int main(int argc, char *argv[])
{
	char header[81];
	float temp1,temp2;
	char filename1[64],filename2[64],outputfilename[64];
	FILE * beamgainfile1,*beamgainfile2,* outputfile;
	int i,j,k;
	int startchannel,endchannel;
	if(argc != 5)
	{
		printf("Usage: beamgains <startchannel> <endchannel> <region1> <region2>\n");
		return EXIT_FAILURE;
	}
	else
	{
		startchannel = atoi(argv[1]);
		endchannel = atoi(argv[2]);
	}
	sprintf(filename1,"beamgains%s.dat",argv[3]);
	sprintf(filename2,"beamgains%s.dat",argv[4]);
	beamgainfile1 = fopen(filename1, "r");
	beamgainfile2 = fopen(filename2, "r");
	sprintf(outputfilename,"beamgainratios_%s_%s.dat",argv[3],argv[4]);
	outputfile = fopen(outputfilename, "w");
	if(beamgainfile1 == NULL || beamgainfile2 == NULL)
		printf("ERROR: Cannot open input file.\n");
	fprintf(outputfile, "#CHAN BEAM0   BEAM1      BEAM2      BEAM3      BEAM4      BEAM5      BEAM6\n");
	printf("Calculating beam gains...\n");
	fgets(header,80,beamgainfile1);
	fgets(header,80,beamgainfile2);

	for(i = startchannel;i < endchannel+1;i++)
	{
		fprintf(outputfile,"%d",i);
		fscanf(beamgainfile1,"%d \n",&k);
		fscanf(beamgainfile2,"%d \n",&k);

		for(j = 0;j < NUM_BEAMS;j++)
		{
			fscanf(beamgainfile1,"%f",&temp1);
			fscanf(beamgainfile2,"%f",&temp2);
			fprintf(outputfile," %1.8f",(temp1/temp2));
		}
		fprintf(outputfile,"\n");

	}	
	fclose(beamgainfile1);
	fclose(beamgainfile2);
	fclose(outputfile);
	printf("Done!\n");
	return EXIT_SUCCESS;
}

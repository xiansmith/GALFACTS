#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"cimafits2.h"
#include"programs/fitsio.h"

int main(int argc,char* argv[])
{
	int beam;
	char *proj_code;
	char *date;
	char *datadir;
	int band;
	int trailer;
	char buf[LINELEN+1];
	if(argc != 7)
	{
		printf("usage: read_cimafits_header <beam> <proj_code> <date> <band> <datadir> <trailing nnnnn>\n");
		return 0;
	}
	else
	{
		beam = atoi(argv[1]);
		proj_code = argv[2];
		date = argv[3];
		band = atoi(argv[4]);
		datadir = argv[5];
		trailer = atoi(argv[6]);
	}	

	FILE *datafile;
	char datafilename[100+1];
	sprintf(datafilename,"%s/%s.%s.b%1ds%1dg0.%5d.fits",datadir,proj_code,date,beam,band,trailer);
	if ( (datafile = fopen(datafilename, "r") ) == NULL )
	{ 
		printf("ERROR: can't open data file for reading '%s'\n", datafilename);
		return 0;
	}
	printf("Opened the datafile:%s\n",datafilename);
	
	int pcount,pcounts,theap,naxis1,naxis2;
	
	int found = FALSE;
	do
	{
		fread(buf,sizeof(char),LINELEN,datafile);
		buf[LINELEN] = '\0';
		if(strncmp(buf,"NAXIS1  ",8)==0)
		{
			sscanf(&buf[10],"%d",&naxis1);
			printf("NAXIS1: %d\n",naxis1);
			found = TRUE;
		}
	}while(!found);

	found = FALSE;
	do
	{
		fread(buf,sizeof(char),LINELEN,datafile);
		buf[LINELEN] = '\0';
		if(strncmp(buf,"NAXIS2  ",8)==0)
		{
			sscanf(&buf[10],"%d",&naxis2);
			printf("NAXIS2: %d\n",naxis2);
			found = TRUE;
		}
	}while(!found);
	found = FALSE;
	do
	{
		fread(buf,sizeof(char),LINELEN,datafile);
		buf[LINELEN] = '\0';
		if(strncmp(buf,"PCOUNT  ",8)==0)
		{
				sscanf(&buf[10],"%d",&pcount);
				printf("PCOUNT: %d\n",pcount);
				found = TRUE;
			}
	}while(!found);

	found = FALSE;
	do
	{
		fread(buf,sizeof(char),LINELEN,datafile);
		buf[LINELEN] = '\0';
		if(strncmp(buf,"THEAP   ",8)==0)
		{
			sscanf(&buf[10],"%d",&theap);
			printf("THEAP: %d\n",theap);
			found = TRUE;
		}
	}while(!found);
	found = FALSE;
	do
	{
		fread(buf,sizeof(char),LINELEN,datafile);
		buf[LINELEN] = '\0';
		if(strncmp(buf,"PCOUNTS ",8)==0)
		{
			sscanf(&buf[10],"%d",&pcounts);
			printf("PCOUNTS: %d\n",pcounts);
			found = TRUE;
		}
	}while(!found);	

	long int offset = MAIN_HEADER + BINTABLE_HEADER;

	cimafits_row c1;
	fseek(datafile,offset+0*naxis1,SEEK_SET);
	fread(&c1,sizeof(cimafits_row),1,datafile);
	printf("Spec,N/A,N/A,num_stokes,num_dumps_row\n");
	printf("%s\n",c1.tdim1);
	fclose(datafile);
	return 1;
}

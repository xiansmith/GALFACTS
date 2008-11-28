#include<stdio.h>
#include"galfactsLib.h"
#include"spec.h"
#include"dtm2spec.h"

#define MAX_ROWS 109
static inline void cnvrt_end_sint(short int *x)
{
	short int result;
	((unsigned char*)&result)[0] = ((unsigned char*)x)[1];	
	((unsigned char*)&result)[1] = ((unsigned char*)x)[0];	
	*x = result;
}

static inline void cnvrt_end_int(int *x)
{
	int result;
	((unsigned char*)&result)[0] = ((unsigned char*)x)[3];	
	((unsigned char*)&result)[1] = ((unsigned char*)x)[2];	
	((unsigned char*)&result)[2] = ((unsigned char*)x)[1];	
	((unsigned char*)&result)[3] = ((unsigned char*)x)[0];	
	*x = result;
}

static inline void cnvrt_end_db(double *x)
{
	double result;
	((unsigned char*)&result)[0] = ((unsigned char*)x)[7];	
	((unsigned char*)&result)[1] = ((unsigned char*)x)[6];	
	((unsigned char*)&result)[2] = ((unsigned char*)x)[5];	
	((unsigned char*)&result)[3] = ((unsigned char*)x)[4];	
	((unsigned char*)&result)[4] = ((unsigned char*)x)[3];	
	((unsigned char*)&result)[5] = ((unsigned char*)x)[2];	
	((unsigned char*)&result)[6] = ((unsigned char*)x)[1];	
	((unsigned char*)&result)[7] = ((unsigned char*)x)[0];	
	*x = result;
}

int main(int argc,char* argv[])
{
	int num_files;
	int beam;
	char *proj_code;
	char *date;
	char *datadir;
	int band;
	int start_file;
	
	if(argc !=8)
	{
		printf("usage: dtm2spec <proj_code> <date> <band> <beam> <datadir> <start_file nnnnn> <num_files>\n");
		return 0;
	}
	else
	{
		proj_code = argv[1];
		date = argv[2];
		band = atoi(argv[3]);
		beam = atoi(argv[4]);
		datadir = argv[5];
                start_file = atoi(argv[6]);
		num_files = atoi(argv[7]);
	}

	FILE *datafile,*specfile;
	char datafilename[100+1],specfilename[40+1];
	
	sprintf(specfilename,"%s.%s.b%1ds%1d.spec",proj_code,date,beam,band);
	if ( (specfile = fopen(specfilename, "wb") ) == NULL )
	{ 
		printf("ERROR: can't open data file for writing '%s'\n", specfilename);
		return 0;
	}

	printf("Size of row struct: %d\n",sizeof(GFLIB_ROW));
	printf("Size of data struct: %d\n",sizeof(GFLIB_DATA));
	printf("Size of stat struct: %d\n",sizeof(GFLIB_PDEVSTAT));
	char buf[LINELEN+1];
	long int offset = MAIN_HEADER + BINTABLE_HEADER;
	int f,g,k,l,found = FALSE;
	for(f=0;f<num_files;f++)
	{
		int naxis1,naxis2;
		GFLIB_ROW row1,row2;
		SpecPointingBlock SPBlock;	
		float Aon,Aoff,Bon,Boff,Uon,Uoff,Von,Voff;

		sprintf(datafilename,"%s/%s_dtm.%s.b%1ds%1dg0.%.5d.fits",datadir,proj_code,date,beam,band,start_file+f);
		if ( (datafile = fopen(datafilename, "r") ) == NULL )
		{ 
			printf("ERROR: can't open data file for reading '%s'\n", datafilename);
			return 0;
		}
		printf("Opened the datafile:%s\n",datafilename);
	
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


		for(g=0;g < naxis2-1;g++)
		{
			printf("Reading row :%d\n",g+1);
			fseek(datafile,offset+g*naxis1,SEEK_SET);
			fread(&row1,sizeof(GFLIB_ROW),1,datafile);
			fread(&row2,sizeof(GFLIB_ROW),1,datafile);
			
			cnvrt_end_db(&row1.RA);
			cnvrt_end_db(&row2.RA);
			cnvrt_end_db(&row1.DEC);
			cnvrt_end_db(&row2.DEC);
			cnvrt_end_db(&row1.AST);
			cnvrt_end_db(&row2.AST);
			for(k=0;k<DUMPS_PER_ROW;k++)
			{
				SPBlock.centralBeam.raj_true_in_hours=(row1.RA + k*(row2.RA-row1.RA)/DUMPS_PER_ROW)/15;
				SPBlock.centralBeam.decj_true_in_degrees=(row1.DEC + k*(row2.DEC-row1.DEC)/DUMPS_PER_ROW);
				SPBlock.centralBeam.arecibo_local_mean_sidereal_time_in_sec=\
				(row1.AST+k*(row2.AST-row1.AST)/DUMPS_PER_ROW);

				cnvrt_end_sint(&row1.staton[k].fftAccum);
				cnvrt_end_sint(&row1.statoff[k].fftAccum);
				
				fwrite(&SPBlock,sizeof(SpecPointingBlock),1,specfile);
				for(l=0;l<RAW_CHANNELS;l++)
				{

					cnvrt_end_int(&row1.dataon[k].A[l]);
					cnvrt_end_int(&row1.dataon[k].B[l]);
					cnvrt_end_int(&row1.dataon[k].U[l]);
					cnvrt_end_int(&row1.dataon[k].V[l]);
					cnvrt_end_int(&row1.dataoff[k].A[l]);
					cnvrt_end_int(&row1.dataoff[k].B[l]);
					cnvrt_end_int(&row1.dataoff[k].U[l]);
					cnvrt_end_int(&row1.dataoff[k].V[l]);
				
					Aon = (float)(row1.dataon[k].A[l]/row1.staton[k].fftAccum);
					Bon = (float)(row1.dataon[k].A[l]/row1.staton[k].fftAccum);
					Uon = (float)(row1.dataon[k].A[l]/row1.staton[k].fftAccum);
					Von = (float)(row1.dataon[k].A[l]/row1.staton[k].fftAccum);
					Aoff = (float)(row1.dataoff[k].A[l]/row1.statoff[k].fftAccum);
					Boff = (float)(row1.dataoff[k].A[l]/row1.statoff[k].fftAccum);
					Uoff = (float)(row1.dataoff[k].A[l]/row1.statoff[k].fftAccum);
					Voff = (float)(row1.dataoff[k].A[l]/row1.statoff[k].fftAccum);

					fwrite(&Aon,sizeof(float),1,specfile);
					fwrite(&Bon,sizeof(float),1,specfile);
					fwrite(&Uon,sizeof(float),1,specfile);
					fwrite(&Von,sizeof(float),1,specfile);
					fwrite(&Aoff,sizeof(float),1,specfile);
					fwrite(&Boff,sizeof(float),1,specfile);
					fwrite(&Uoff,sizeof(float),1,specfile);
					fwrite(&Voff,sizeof(float),1,specfile);
				}//l loop for each channel
			}//k loop num dumps
			cnvrt_end_db(&row2.RA);
			cnvrt_end_db(&row2.DEC);
			cnvrt_end_db(&row2.AST);

		}//naxis2 loop g
		fclose(datafile);
	}//num files loop f
	fclose(specfile);
}

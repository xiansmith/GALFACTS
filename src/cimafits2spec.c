/*
	Conversion utility for converting Mock spectrometer output files
	to 2 streams (Low time high spectral resolution and vice versa).
	The output files are in .spec format compatible with WAPP output
	files.

	sguram 15 Nov 2008
*/

#include<stdio.h>
#include<string.h>
#include"cimafits2.h"
#include"spec.h"
#include"programs/fitsio.h"
#include<malloc.h>
#include<stdlib.h>
//#define IGNORE_ROWS 3
//#define STATUS_HEAP_SIZE 2610000
#define MAX_ROWS 261
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
	int lths_time_comp;
	int lths_spec_comp;
	int htls_time_comp;
	int htls_spec_comp;
	char *proj_code;
	char *date;
	char *datadir;
	int band;
	int start_file;
	int ignore_rows;
	float dec_min,dec_max;	
	if(argc !=10)
	{
		printf("usage: cimafits2spec <num_files> <beam> \n <lths_time_comp> <htls_spec_comp> \
		<proj_code> \n <date> <band> <datadir> <start_file nnnnn> \n");
		return 0;
	}
	else
	{
		num_files = atoi(argv[1]);
		beam = atoi(argv[2]);
		lths_time_comp = atoi(argv[3]);
		lths_spec_comp = 1;
//		htls_time_comp = atoi(argv[5]);		
		htls_time_comp = 1;
		htls_spec_comp = atoi(argv[4]);
		proj_code = argv[5];
		date = argv[6];
		band = atoi(argv[7]);
		datadir = argv[8];
        start_file = atoi(argv[9]);
//      ignore_rows = atoi(argv[10]);
		ignore_rows = 0;
//      dec_min = atof(argv[11]);
//      dec_max = atof(argv[12]);
	}	

	
	int config_lths_not_written=1,config_htls_not_written=1;

	FILE *datafile,*lths_file,*htls_file;
	char datafilename[100+1],lthsfilename[40+1],htlsfilename[40+1];
	
	sprintf(htlsfilename,"%s.%s.b%1ds%1d.htls.spec",proj_code,date,beam,band);
	if ( (htls_file = fopen(htlsfilename, "wb") ) == NULL )
	{ 
		printf("ERROR: can't open data file for writing '%s'\n", htlsfilename);
		return 0;
	}

	sprintf(lthsfilename,"%s.%s.b%1ds%1d.lths.spec",proj_code,date,beam,band);
	if ( (lths_file = fopen(lthsfilename, "wb") ) == NULL )
	{ 
		printf("ERROR: can't open data file for writing '%s'\n", lthsfilename);
		return 0;
	}

	char buf[81];
	cimafits_row c1,c2;
	long int offset = MAIN_HEADER + BINTABLE_HEADER;
	int size_pstat = sizeof(pdev_stat);
	int size_pdatum = sizeof(pdev_datum);
	int size_int = sizeof(int);
//	int size_long = sizeof(long int);
//	pdev_stat pstat[DUMPS_PER_ROW];
//	pdev_datum pdatum[DUMPS_PER_ROW];
	pdev_stat *pstat;
	pdev_datum *pdatum;

	pstat = (pdev_stat *)malloc(size_pstat*DUMPS_PER_ROW);
	pdatum = (pdev_datum *)malloc(size_pdatum*DUMPS_PER_ROW);

//	nSpecPointingBlock *spointing_htls,*spointing_lths;	
	SpecPointingBlock *spointing_htls,*spointing_lths;
//	spointing_htls = (nSpecPointingBlock *)malloc(sizeof(nSpecPointingBlock)*DUMPS_PER_ROW/(2*htls_time_comp));
//	spointing_lths = (nSpecPointingBlock *)malloc(sizeof(nSpecPointingBlock)*DUMPS_PER_ROW/(2*lths_time_comp));
	spointing_htls = (SpecPointingBlock *)malloc(sizeof(SpecPointingBlock)*DUMPS_PER_ROW/(2*htls_time_comp));
	spointing_lths = (SpecPointingBlock *)malloc(sizeof(SpecPointingBlock)*DUMPS_PER_ROW/(2*lths_time_comp));	

	Spec_PolSet *spolset_htls_on,*spolset_lths_on;
	Spec_PolSet *spolset_htls_on_final,*spolset_lths_on_final;
	Spec_PolSet *spolset_htls_off,*spolset_lths_off;
	Spec_PolSet *spolset_htls_off_final,*spolset_lths_off_final;
	spolset_htls_on = (Spec_PolSet *)malloc(sizeof(Spec_PolSet)*(DUMPS_PER_ROW/(2*htls_time_comp)));
	spolset_htls_off = (Spec_PolSet *)malloc(sizeof(Spec_PolSet)*(DUMPS_PER_ROW/(2*htls_time_comp)));
	spolset_lths_on = (Spec_PolSet *)malloc(sizeof(Spec_PolSet)*(DUMPS_PER_ROW/(2*lths_time_comp)));
	spolset_lths_off = (Spec_PolSet *)malloc(sizeof(Spec_PolSet)*(DUMPS_PER_ROW/(2*lths_time_comp)));

	spolset_htls_on_final = (Spec_PolSet *)malloc(sizeof(Spec_PolSet)*(DUMPS_PER_ROW/(2*htls_time_comp)));
	spolset_htls_off_final = (Spec_PolSet *)malloc(sizeof(Spec_PolSet)*(DUMPS_PER_ROW/(2*htls_time_comp)));
	spolset_lths_on_final = (Spec_PolSet *)malloc(sizeof(Spec_PolSet)*(DUMPS_PER_ROW/(2*lths_time_comp)));
	spolset_lths_off_final = (Spec_PolSet *)malloc(sizeof(Spec_PolSet)*(DUMPS_PER_ROW/(2*lths_time_comp)));


	float *XXon_htls,*XXoff_htls,*YYon_htls,*YYoff_htls;
	float *XYon_htls,*XYoff_htls,*YXon_htls,*YXoff_htls;
	float *XXon_lths,*XXoff_lths,*YYon_lths,*YYoff_lths;
	float *XYon_lths,*XYoff_lths,*YXon_lths,*YXoff_lths;


	XXon_htls = (float *)calloc((RAW_CHANNELS/htls_spec_comp),sizeof(float));
	XXoff_htls = (float *)calloc((RAW_CHANNELS/htls_spec_comp),sizeof(float));	
	XYon_htls = (float *)calloc((RAW_CHANNELS/htls_spec_comp),sizeof(float));
	XYoff_htls = (float *)calloc((RAW_CHANNELS/htls_spec_comp),sizeof(float));
	YXon_htls = (float *)calloc((RAW_CHANNELS/htls_spec_comp),sizeof(float));
	YXoff_htls = (float *)calloc((RAW_CHANNELS/htls_spec_comp),sizeof(float));	
	YYon_htls = (float *)calloc((RAW_CHANNELS/htls_spec_comp),sizeof(float));	
	YYoff_htls = (float *)calloc((RAW_CHANNELS/htls_spec_comp),sizeof(float));	

	XXon_lths = (float *)calloc((RAW_CHANNELS/lths_spec_comp),sizeof(float));
	XXoff_lths = (float *)calloc((RAW_CHANNELS/lths_spec_comp),sizeof(float));	
	XYon_lths = (float *)calloc((RAW_CHANNELS/lths_spec_comp),sizeof(float));
	XYoff_lths = (float *)calloc((RAW_CHANNELS/lths_spec_comp),sizeof(float));
	YXon_lths = (float *)calloc((RAW_CHANNELS/lths_spec_comp),sizeof(float));
	YXoff_lths = (float *)calloc((RAW_CHANNELS/lths_spec_comp),sizeof(float));	
	YYon_lths = (float *)calloc((RAW_CHANNELS/lths_spec_comp),sizeof(float));	
	YYoff_lths = (float *)calloc((RAW_CHANNELS/lths_spec_comp),sizeof(float));
			
	int h;
	for(h = 0;h < DUMPS_PER_ROW/(htls_time_comp*2);h++)			
	{
		spolset_htls_on[h].A = (unsigned int *)calloc((RAW_CHANNELS),size_int);
		spolset_htls_on[h].B = (unsigned int *)calloc((RAW_CHANNELS),size_int);
		spolset_htls_on[h].V = ( int *)calloc((RAW_CHANNELS),size_int);
		spolset_htls_on[h].U = ( int *)calloc((RAW_CHANNELS),size_int);
		spolset_htls_off[h].A = (unsigned  int *)calloc((RAW_CHANNELS),size_int);
		spolset_htls_off[h].B = (unsigned  int *)calloc((RAW_CHANNELS),size_int);
		spolset_htls_off[h].V = (int *)calloc((RAW_CHANNELS),size_int);
		spolset_htls_off[h].U = (int *)calloc((RAW_CHANNELS),size_int);
		spolset_htls_on_final[h].A = (unsigned int *)calloc((RAW_CHANNELS/htls_spec_comp),size_int);
		spolset_htls_on_final[h].B = (unsigned int *)calloc((RAW_CHANNELS/htls_spec_comp),size_int);
		spolset_htls_on_final[h].V = (int *)calloc((RAW_CHANNELS/htls_spec_comp),size_int);
		spolset_htls_on_final[h].U = (int *)calloc((RAW_CHANNELS/htls_spec_comp),size_int);
		spolset_htls_off_final[h].A = (unsigned int *)calloc((RAW_CHANNELS/htls_spec_comp),size_int);
		spolset_htls_off_final[h].B = (unsigned int *)calloc((RAW_CHANNELS/htls_spec_comp),size_int);
		spolset_htls_off_final[h].V = (int *)calloc((RAW_CHANNELS/htls_spec_comp),size_int);
		spolset_htls_off_final[h].U = (int *)calloc((RAW_CHANNELS/htls_spec_comp),size_int);
	}

	for(h = 0;h < DUMPS_PER_ROW/(lths_time_comp*2);h++)			
	{
		spolset_lths_on[h].A = (unsigned int *)calloc((RAW_CHANNELS),size_int);
		spolset_lths_on[h].B = (unsigned int *)calloc((RAW_CHANNELS),size_int);
		spolset_lths_on[h].V = (int *)calloc((RAW_CHANNELS),size_int);
		spolset_lths_on[h].U = (int *)calloc((RAW_CHANNELS),size_int);
		spolset_lths_off[h].A = (unsigned int *)calloc((RAW_CHANNELS),size_int);
		spolset_lths_off[h].B = (unsigned  int *)calloc((RAW_CHANNELS),size_int);
		spolset_lths_off[h].V = ( int *)calloc((RAW_CHANNELS),size_int);
		spolset_lths_off[h].U = ( int *)calloc((RAW_CHANNELS),size_int);
		spolset_lths_on_final[h].A = (unsigned int *)calloc((RAW_CHANNELS/lths_spec_comp),size_int);
		spolset_lths_on_final[h].B = (unsigned  int *)calloc((RAW_CHANNELS/lths_spec_comp),size_int);
		spolset_lths_on_final[h].V = (int *)calloc((RAW_CHANNELS/lths_spec_comp),size_int);
		spolset_lths_on_final[h].U = (int *)calloc((RAW_CHANNELS/lths_spec_comp),size_int);
		spolset_lths_off_final[h].A = (unsigned int *)calloc((RAW_CHANNELS/lths_spec_comp),size_int);
		spolset_lths_off_final[h].B = (unsigned int *)calloc((RAW_CHANNELS/lths_spec_comp),size_int);
		spolset_lths_off_final[h].V = (int *)calloc((RAW_CHANNELS/lths_spec_comp),size_int);
		spolset_lths_off_final[h].U = (int *)calloc((RAW_CHANNELS/lths_spec_comp),size_int);
	}

	int f,found = FALSE;
	for(f=0;f<num_files;f++)
	{
		int pcount,pcounts,theap,naxis1,naxis2;
			
		sprintf(datafilename,"%s/%s.%s.b%1ds%1dg0.%.5d.fits",datadir,proj_code,date,beam,band,start_file+f);
		if ( (datafile = fopen(datafilename, "r") ) == NULL )
		{ 
			printf("ERROR: can't open data file for reading '%s'\n", datafilename);
			return 0;
		}
		printf("#Opened the datafile:%s\n",datafilename);
	
		do
		{
			fread(buf,sizeof(char),LINELEN,datafile);
			buf[LINELEN] = '\0';
			if(strncmp(buf,"NAXIS1  ",8)==0)
			{
				sscanf(&buf[10],"%d",&naxis1);
				printf("#NAXIS1: %d\n",naxis1);
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
				printf("#NAXIS2: %d\n",naxis2);
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
				printf("#PCOUNT: %d\n",pcount);
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
				printf("#THEAP: %d\n",theap);
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
				printf("#PCOUNTS: %d\n",pcounts);
				found = TRUE;
			}
		}while(!found);		
	

		int g,i;
		int naxis2_fix;	
//		if(f == 0)
//			g = ignore_rows;
//		else
			g = 0;
			naxis2_fix = naxis2-1;
//		if(f == (num_files-1))
//			naxis2_fix = naxis2 - ignore_rows;
//		else
//			naxis2_fix = naxis2 - 1;
//		printf("Size of cimafits row %d\n",sizeof(cimafits_row));
//		int outside_bound;
		for(;g < naxis2_fix;g++)
		{
			printf("#Reading row %d\n",g+1);
			fseek(datafile,offset+g*naxis1,SEEK_SET);
			fread(&c1,sizeof(cimafits_row),1,datafile);
			fread(&c2,sizeof(cimafits_row),1,datafile);
			
			cnvrt_end_db(&c1.crval2);
			cnvrt_end_db(&c2.crval2);
			cnvrt_end_db(&c1.crval3);
			cnvrt_end_db(&c2.crval3);
			cnvrt_end_db(&c1.lst);
			cnvrt_end_db(&c2.lst);
			cnvrt_end_db(&c1.mjdxxobs);
			cnvrt_end_db(&c1.azimuth);
			cnvrt_end_db(&c1.elevatio);
			cnvrt_end_db(&c2.azimuth);
			cnvrt_end_db(&c2.elevatio);
			cnvrt_end_db(&c1.req_raj);
			cnvrt_end_db(&c1.req_decj);
			cnvrt_end_db(&c2.req_raj);
			cnvrt_end_db(&c2.req_decj);
			cnvrt_end_db(&c1.alfa_ang);
			
//			printf("RA:%f DEC:%f\n",c1.crval2,c1.crval3);
		
//			if(c1.crval3 > dec_max || c1.crval3 < dec_min) //may not always work keep an eye
//			{
//				outside_bound = TRUE;
//				printf("Row outside normal DEC range\n");
//			}
//			else
//				outside_bound = FALSE;

			fseek(datafile,offset+theap+g*DUMPS_PER_ROW*size_pstat,SEEK_SET);
			for(i = 0;i < DUMPS_PER_ROW;i++)
			{
				//fread(&pstat[i],size_pstat,DUMPS_PER_ROW,datafile);
				fread(&pstat[i],size_pstat,1,datafile);
			}
			fseek(datafile,offset+theap+MAX_ROWS*size_pstat*DUMPS_PER_ROW+g*DUMPS_PER_ROW*size_pdatum,SEEK_SET); 
			for(i = 0;i < DUMPS_PER_ROW;i++)
			{
				//fread(&pdatum[i],size_pdatum,DUMPS_PER_ROW,datafile);
				fread(&pdatum[i],size_pdatum,1,datafile);
			}
			FILE * htls_cfg_file;
			
			if(config_htls_not_written)
			{
				sprintf(htlsfilename,"%s.%s.b%ds%d.htls.spec_cfg",proj_code,date,beam,band);
				if ( (htls_cfg_file = fopen(htlsfilename, "w") ) == NULL )
				{ 
					printf("ERROR: can't open config file '%s'\n", htlsfilename);
					return 0;
				}

				cnvrt_end_db(&c1.cdelt5);
				fprintf(htls_cfg_file,"%f\n",c1.cdelt5*htls_time_comp*1000);
				fprintf(htls_cfg_file,"%i\n",RAW_CHANNELS*4/htls_spec_comp);
				cnvrt_end_db(&c1.crval1);
				fprintf(htls_cfg_file,"%f\n",c1.crval1/1000000);
				cnvrt_end_db(&c1.cdelt1);
				fprintf(htls_cfg_file,"%f\n",-1*c1.cdelt1*htls_spec_comp/1000);
				fprintf(htls_cfg_file,"1 %i 4 2 0\n",RAW_CHANNELS/htls_spec_comp);
				fprintf(htls_cfg_file,"%s\n",proj_code);
				fprintf(htls_cfg_file,"%f\n",c1.mjdxxobs);
				fprintf(htls_cfg_file,"AO\n");
				fprintf(htls_cfg_file,"Integration time (ms):%f\n",c1.cdelt5*htls_time_comp*1000);
				fprintf(htls_cfg_file,"MJD: %f\n",c1.mjdxxobs);
				fprintf(htls_cfg_file,"Center freq (MHz): %f\n",c1.crval1/1000000);
				fprintf(htls_cfg_file,"Channel band (kHz): %f\n",-1*c1.cdelt1*htls_spec_comp/1000);
				fprintf(htls_cfg_file,"Number of channels/record: %d\n",RAW_CHANNELS/htls_spec_comp);
				fprintf(htls_cfg_file,"RA at start (degrees): %f\n",c1.crval2);
				fprintf(htls_cfg_file,"DEC at start (degrees): %f\n",c1.crval3);
				fprintf(htls_cfg_file,"LST at start (seconds): %f\n",c1.lst);
				fprintf(htls_cfg_file,"ALFA angle (degrees) at start: %f\n",c1.alfa_ang);
				fprintf(htls_cfg_file,"Project ID: %s\n",proj_code);
				fclose(htls_cfg_file);
				config_htls_not_written = 0;
			}


			FILE * lths_cfg_file;
			
			if(config_lths_not_written)
			{
				sprintf(lthsfilename,"%s.%s.b%ds%d.lths.spec_cfg",proj_code,date,beam,band);
				if ( (lths_cfg_file = fopen(lthsfilename, "w") ) == NULL )
				{ 
					printf("ERROR: can't open config file '%s'\n", lthsfilename);
					return 0;
				}

				fprintf(lths_cfg_file,"%f\n",c1.cdelt5*lths_time_comp*1000);
				fprintf(lths_cfg_file,"%i\n",RAW_CHANNELS*4/lths_spec_comp);
				fprintf(lths_cfg_file,"%f\n",c1.crval1/1000000);
				fprintf(lths_cfg_file,"%f\n",-1*c1.cdelt1*lths_spec_comp/1000);
				fprintf(lths_cfg_file,"1 %i 4 2 0\n",RAW_CHANNELS/lths_spec_comp);
				fprintf(lths_cfg_file,"%s\n",proj_code);
				fprintf(lths_cfg_file,"%f\n",c1.mjdxxobs);
				fprintf(lths_cfg_file,"AO\n");
				fprintf(lths_cfg_file,"Integration time (ms):%f\n",c1.cdelt5*lths_time_comp*1000);
				fprintf(lths_cfg_file,"MJD: %f\n",c1.mjdxxobs);
				fprintf(lths_cfg_file,"Center freq (MHz): %f\n",c1.crval1/1000000);
				fprintf(lths_cfg_file,"Channel band (kHz): %f\n",-1*c1.cdelt1*lths_spec_comp/1000);
				fprintf(lths_cfg_file,"Number of channels/record: %d\n",RAW_CHANNELS/lths_spec_comp);
				fprintf(lths_cfg_file,"RA at start (degrees): %f\n",c1.crval2);
				fprintf(lths_cfg_file,"DEC at start (degrees): %f\n",c1.crval3);
				fprintf(lths_cfg_file,"LST at start (seconds): %f\n",c1.lst);
				fprintf(lths_cfg_file,"ALFA angle (degrees) at start: %f\n",c1.alfa_ang);
				fprintf(lths_cfg_file,"Project ID: %s\n",proj_code);
				fclose(lths_cfg_file);
				config_lths_not_written = 0;
			}

			//init to zero
			for(h = 0;h < DUMPS_PER_ROW/(htls_time_comp*2);h++)			
			{
				int l;
				for(l=0;l<RAW_CHANNELS;l++)
				{
					spolset_htls_on[h].A[l] =0;
					spolset_htls_on[h].B[l] =0;
					spolset_htls_on[h].U[l] =0;
					spolset_htls_on[h].V[l] =0;
					spolset_htls_off[h].A[l] =0;
					spolset_htls_off[h].B[l] =0;
					spolset_htls_off[h].U[l] =0;
					spolset_htls_off[h].V[l] =0;
				}
				spolset_htls_on[h].fft_weight = 0;	
				spolset_htls_off[h].fft_weight = 0;	
				for(l=0;l<RAW_CHANNELS/htls_spec_comp;l++)
				{
					spolset_htls_on_final[h].A[l] =0;
					spolset_htls_on_final[h].B[l] =0;
					spolset_htls_on_final[h].U[l] =0;
					spolset_htls_on_final[h].V[l] =0;
					spolset_htls_off_final[h].A[l] =0;
					spolset_htls_off_final[h].B[l] =0;
					spolset_htls_off_final[h].U[l] =0;
					spolset_htls_off_final[h].V[l] =0;
				}	
				spolset_htls_on_final[h].fft_weight = 0;	
				spolset_htls_off_final[h].fft_weight = 0;	
			}
			


			int h1 =0,h2 =0,flag=0;
			int num_on=0,num_off=0,onoff=0,onoff_count=0,oncount=0,offcount=0;
			for(h = 0;h < DUMPS_PER_ROW/(htls_time_comp*2);h++)			
			{
//				printf("h %d\n",h);
/*				spointing_htls[h].raj_true_in_degrees = c1.crval2 + h*(c2.crval2-c1.crval2)*(2*htls_time_comp)/(DUMPS_PER_ROW);
				spointing_htls[h].decj_true_in_degrees = c1.crval3 + h*(c2.crval3-c1.crval3)*(2*htls_time_comp)/(DUMPS_PER_ROW);
				spointing_htls[h].arecibo_local_mean_sidereal_time_in_sec = c1.lst + h*(c2.lst-c1.lst)*(2*htls_time_comp)/(DUMPS_PER_ROW);
				spointing_htls[h].az_cur_in_degrees = c1.azimuth + h*(c2.azimuth-c1.azimuth)*(2*htls_time_comp)/(DUMPS_PER_ROW);
				spointing_htls[h].za_cur_in_degrees = c1.elevatio + h*(c2.elevatio-c1.elevatio)*(2*htls_time_comp)/(DUMPS_PER_ROW);
				spointing_htls[h].raj_requested_in_degrees = c1.req_raj + h*(c2.req_raj-c1.req_raj)*(2*htls_time_comp)/(DUMPS_PER_ROW);
				spointing_htls[h].decj_requested_in_degrees = c1.req_decj + h*(c2.req_decj-c1.req_decj)*(2*htls_time_comp)/(DUMPS_PER_ROW);
				spointing_htls[h].mjd_last_five_digits = c1.mjdxxobs;
				spointing_htls[h].alfa_rotation_angle_in_degrees = c1.alfa_ang;
*/
				spointing_htls[h].centralBeam.raj_true_in_hours=c1.crval2 + h*(c2.crval2-c1.crval2)*(2*htls_time_comp)/(DUMPS_PER_ROW)/15;
				spointing_htls[h].centralBeam.decj_true_in_degrees=c1.crval3 + h*(c2.crval3-c1.crval3)*(2*htls_time_comp)/(DUMPS_PER_ROW);
				spointing_htls[h].centralBeam.atlantic_solar_time_now_in_sec=\
				c1.lst + h*(c2.lst-c1.lst)*(2*htls_time_comp)/(DUMPS_PER_ROW);			
//				printf("%f %f %f\n",spointing_htls[h].centralBeam.raj_true_in_hours,spointing_htls[h].centralBeam.decj_true_in_degrees, \
				spointing_htls[h].centralBeam.atlantic_solar_time_now_in_sec);

				int k,l;
				num_on=0,num_off=0;
				for(k = 0;k < htls_time_comp*2;k++)
				{
					cnvrt_end_sint(&pstat[h*htls_time_comp*2+k].calOn);
//					cnvrt_end_sint(&pstat[h*htls_time_comp*2+k+1].calOn);
//					printf("Cal:%d loopcount: %d\ h1: %d h2:%d\n",pstat[h*htls_time_comp*2+k].calOn,h*htls_time_comp*2+k,h1,h2);
/*					if(outside_bound)
					{
						//onoff is 0 if writing to num_off
						if(!onoff)
						{
							cnvrt_end_sint(&pstat[h*htls_time_comp*2+k].fftAccum);
//							if(pstat[h*htls_time_comp*2+k].fftAccum != 168)
//								printf("fft:%d\n",pstat[h*htls_time_comp*2+k].fftAccum);
							spolset_htls_off[h1].fft_weight += (int)pstat[h*htls_time_comp*2+k].fftAccum;
							for(l=0;l<RAW_CHANNELS;l++)
							{
								cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_A[l]);
								cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_B[l]);
								cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_U[l]);
								cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_V[l]);
								spolset_htls_off[h1].A[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_A[l]);
								spolset_htls_off[h1].B[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_B[l]);
								spolset_htls_off[h1].U[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_U[l]);
								spolset_htls_off[h1].V[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_V[l]);
							}//l loop htls	
							num_off++;
						}//onoff if
						else
						{
							cnvrt_end_sint(&pstat[h*htls_time_comp*2+k].fftAccum);
//							if(pstat[h*htls_time_comp*2+k].fftAccum != 168)
//								printf("fft:%d\n",pstat[h*htls_time_comp*2+k].fftAccum);
							spolset_htls_off[h1].fft_weight += (int)pstat[h*htls_time_comp*2+k].fftAccum;
							for(l=0;l<RAW_CHANNELS;l++)
							{
								cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_A[l]);
								cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_B[l]);
								cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_U[l]);
								cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_V[l]);
								spolset_htls_on[h2].A[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_A[l]);
								spolset_htls_on[h2].B[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_B[l]);
								spolset_htls_on[h2].U[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_U[l]);
								spolset_htls_on[h2].V[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_V[l]);
							}//l loop htls	
							num_on++;
						}//onoff else
						onoff_count++;
					}
					else*/
					if(!onoff)					
//					if(pstat[h*htls_time_comp*2+k].calOn == 0)
					{
//						flag = 0;
						//int l;
						if(h1==DUMPS_PER_ROW/(htls_time_comp*2))
							printf("Problematic h1 %d\n",h1);

						cnvrt_end_sint(&pstat[h*htls_time_comp*2+k].fftAccum);
//						if(pstat[h*htls_time_comp*2+k].fftAccum != 168)
//							printf("fft:%d\n",pstat[h*htls_time_comp*2+k].fftAccum);
						spolset_htls_off[h1].fft_weight += (int)pstat[h*htls_time_comp*2+k].fftAccum;
						for(l=0;l<RAW_CHANNELS;l++)
						{
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_A[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_B[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_U[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_V[l]);
							spolset_htls_off[h1].A[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_A[l]);
							spolset_htls_off[h1].B[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_B[l]);
							spolset_htls_off[h1].U[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_U[l]);
							spolset_htls_off[h1].V[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_V[l]);
						}//l loop htls	
						num_off++;
						offcount++;
						if(offcount == 10)
						{
							onoff = 1;
							offcount = 0;
						}
					}//fi
/*					else if(pstat[h*htls_time_comp*2+k].calOn == 1 && flag == 0)
					{
						flag = 1;
						cnvrt_end_sint(&pstat[h*htls_time_comp*2+k].fftAccum);
//						if(pstat[h*htls_time_comp*2+k].fftAccum != 168)	
//							printf("fft:%d\n",pstat[h*htls_time_comp*2+k].fftAccum);
						spolset_htls_off[h1].fft_weight += (int)pstat[h*htls_time_comp*2+k].fftAccum;
						for(l=0;l<RAW_CHANNELS;l++)
						{
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_A[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_B[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_U[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_V[l]);
							spolset_htls_off[h1].A[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_A[l]);
							spolset_htls_off[h1].B[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_B[l]);
							spolset_htls_off[h1].U[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_U[l]);
							spolset_htls_off[h1].V[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_V[l]);
						}
						num_off++;
					}
					else if(pstat[h*htls_time_comp*2+k].calOn == 1  && pstat[h*htls_time_comp*2+k+1].calOn == 0)
					{
						flag = 1;
						cnvrt_end_sint(&pstat[h*htls_time_comp*2+k].fftAccum);
//						if(pstat[h*htls_time_comp*2+k].fftAccum != 168)	
//							printf("fft:%d\n",pstat[h*htls_time_comp*2+k].fftAccum);
						spolset_htls_on[h2].fft_weight += (int)pstat[h*htls_time_comp*2+k].fftAccum;
						for(l=0;l<RAW_CHANNELS;l++)
						{
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_A[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_B[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_U[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_V[l]);
							spolset_htls_on[h2].A[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_A[l]);
							spolset_htls_on[h2].B[l] += (unsigned int)(pdatum[h*htls_time_comp*2+k].pol_B[l]);
							spolset_htls_on[h2].U[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_U[l]);
							spolset_htls_on[h2].V[l] += (int)(pdatum[h*htls_time_comp*2+k].stokes_V[l]);
						}
						num_on++;
					}
*/					else
					{
						if(h2==DUMPS_PER_ROW/(htls_time_comp*2))
							printf("Problematic h2 %d\n",h2);

						cnvrt_end_sint(&pstat[h*htls_time_comp*2+k].fftAccum);
//						if(pstat[h*htls_time_comp*2+k].fftAccum != 168)	
//							printf("fft:%d\n",pstat[h*htls_time_comp*2+k].fftAccum);
						spolset_htls_on[h2].fft_weight += (int)pstat[h*htls_time_comp*2+k].fftAccum;
						for(l=0;l<RAW_CHANNELS;l++)
						{
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_A[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].pol_B[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_U[l]);
							cnvrt_end_sint(&pdatum[h*htls_time_comp*2+k].stokes_V[l]);
							spolset_htls_on[h2].A[l] += (unsigned int)pdatum[h*htls_time_comp*2+k].pol_A[l];
							spolset_htls_on[h2].B[l] += (unsigned  int)pdatum[h*htls_time_comp*2+k].pol_B[l];
							spolset_htls_on[h2].U[l] += (int)pdatum[h*htls_time_comp*2+k].stokes_U[l];
							spolset_htls_on[h2].V[l] += (int)pdatum[h*htls_time_comp*2+k].stokes_V[l];
						}//l loop htls	
						num_on++;
						oncount++;
						if(oncount == 10)
						{
							oncount = 0;
							onoff = 0;
						}
//						flag = 1;
					}//final else	
					cnvrt_end_sint(&pstat[h*htls_time_comp*2+k+1].calOn);
/*					if(onoff_count == 20) //cal on off cycle keep an eye on this hard coded value
					{
						if(onoff)
							onoff = 0;
						else
							onoff = 1;
						onoff_count = 0;
					}
*/
/*					if(spolset_htls_on[h2].A[100])
					printf("CalON %d %d %d %d\n",spolset_htls_on[h2].A[100],spolset_htls_on[h2].B[100],spolset_htls_on[h2].U[100],spolset_htls_on[h2].V[100]);
					if(spolset_htls_off[h1].A[100])
					printf("CalOFF %d %d %d %d\n",spolset_htls_off[h1].A[100],spolset_htls_off[h1].B[100],spolset_htls_off[h1].U[100],spolset_htls_off[h1].V[100]);
//					printf(" numon : %d numoff %d h1 %d h2 %d\n",num_on,num_off,h1,h2);
*/					if(num_on == htls_time_comp)
					{
	//					printf("num on %d \n",num_on);
//						for(l=0;l<RAW_CHANNELS;l++)
//						{
//							spolset_htls_on[h2].A[l] /= spolset_htls_on[h2].fft_weight;
//							spolset_htls_on[h2].B[l] /= spolset_htls_on[h2].fft_weight;
//							spolset_htls_on[h2].U[l] /= spolset_htls_on[h2].fft_weight;
//							spolset_htls_on[h2].V[l] /= spolset_htls_on[h2].fft_weight;
//						}
						num_on = 0;
						h2++;
					}
					if(num_off == htls_time_comp)
					{
	//					printf("num off %d \n",num_off);
//						for(l=0;l<RAW_CHANNELS;l++)
//						{
//							spolset_htls_off[h1].A[l] /= spolset_htls_off[h1].fft_weight;
//							spolset_htls_off[h1].B[l] /= spolset_htls_off[h1].fft_weight;
//							spolset_htls_off[h1].U[l] /= spolset_htls_off[h1].fft_weight;
//							spolset_htls_off[h1].V[l] /= spolset_htls_off[h1].fft_weight;
//						}
						num_off = 0;
						h1++;
					}
				}//k loop htls
				

				//printf("DIAGNOSTIC:h1:%d,h2:%d\n",h1,h2);
			}//h loop for htls


			//init to zero
			for(h = 0;h < DUMPS_PER_ROW/(lths_time_comp*2);h++)			
			{
				int l;
				for(l=0;l<RAW_CHANNELS;l++)
				{
					spolset_lths_on[h].A[l] =0;
					spolset_lths_on[h].B[l] =0;
					spolset_lths_on[h].U[l] =0;
					spolset_lths_on[h].V[l] =0;
					spolset_lths_off[h].A[l] =0;
					spolset_lths_off[h].B[l] =0;
					spolset_lths_off[h].U[l] =0;
					spolset_lths_off[h].V[l] =0;
				}	
				spolset_lths_on[h].fft_weight = 0;	
				spolset_lths_off[h].fft_weight = 0;	
				for(l=0;l<RAW_CHANNELS/lths_spec_comp;l++)
				{
					spolset_lths_on_final[h].A[l] =0;
					spolset_lths_on_final[h].B[l] =0;
					spolset_lths_on_final[h].U[l] =0;
					spolset_lths_on_final[h].V[l] =0;
					spolset_lths_off_final[h].A[l] =0;
					spolset_lths_off_final[h].B[l] =0;
					spolset_lths_off_final[h].U[l] =0;
					spolset_lths_off_final[h].V[l] =0;
				}	
				spolset_lths_on_final[h].fft_weight = 0;	
				spolset_lths_off_final[h].fft_weight = 0;	
			}

			h1 =0;
			h2 =0;
//			flag=0;
			onoff = 0;
			oncount = 0;
			offcount = 0;
			for(h = 0;h < DUMPS_PER_ROW/(lths_time_comp*2);h++)			
			{
/*				spointing_lths[h].raj_true_in_degrees = c1.crval2 + h*(c2.crval2-c1.crval2)*(2*lths_time_comp)/(DUMPS_PER_ROW);
				spointing_lths[h].decj_true_in_degrees = c1.crval3 + h*(c2.crval3-c1.crval3)*(2*lths_time_comp)/(DUMPS_PER_ROW);
				spointing_lths[h].arecibo_local_mean_sidereal_time_in_sec = c1.lst + h*(c2.lst-c1.lst)*(2*lths_time_comp)/(DUMPS_PER_ROW);
				spointing_lths[h].az_cur_in_degrees = c1.azimuth + h*(c2.azimuth-c1.azimuth)*(2*lths_time_comp)/(DUMPS_PER_ROW);
				spointing_lths[h].za_cur_in_degrees = c1.elevatio + h*(c2.elevatio-c1.elevatio)*(2*lths_time_comp)/(DUMPS_PER_ROW);
				spointing_lths[h].raj_requested_in_degrees = c1.req_raj + h*(c2.req_raj-c1.req_raj)*(2*lths_time_comp)/(DUMPS_PER_ROW);
				spointing_lths[h].decj_requested_in_degrees = c1.req_decj + h*(c2.req_decj-c1.req_decj)*(2*lths_time_comp)/(DUMPS_PER_ROW);
				spointing_lths[h].mjd_last_five_digits = c1.mjdxxobs;
				spointing_lths[h].alfa_rotation_angle_in_degrees = c1.alfa_ang;
*/
				spointing_lths[h].centralBeam.raj_true_in_hours=c1.crval2 + h*(c2.crval2-c1.crval2)*(2*lths_time_comp)/(DUMPS_PER_ROW)/15;
				spointing_lths[h].centralBeam.decj_true_in_degrees=c1.crval3 + h*(c2.crval3-c1.crval3)*(2*lths_time_comp)/(DUMPS_PER_ROW);
				spointing_lths[h].centralBeam.atlantic_solar_time_now_in_sec=\
				c1.lst + h*(c2.lst-c1.lst)*(2*lths_time_comp)/(DUMPS_PER_ROW);		
				//printf("%f %f\n",spointing_lths[h].raj_true_in_degrees,spointing_lths[h].decj_true_in_degrees);
				//fwrite(&spointing_lths[h],sizeof(nSpecPointingBlock),1,lths_file);				

				int k,l;
				num_on=0;
				num_off=0;
//				oncount = 0;
//				offcount = 0;
//				onoff = 0;
				for(k = 0;k < lths_time_comp*2;k++)
				{

//					printf("Cal:%d loopcount: %d\ h1: %d h2:%d\n",pstat[h*lths_time_comp*2+k].calOn,h*lths_time_comp*2+k,h1,h2);
					//cnvrt_end_sint(&pstat[h*lths_time_comp+k].calOn)
/*					if(outside_bound)
					{
						//onoff is 0 if writing to num_off
						if(!onoff)
						{
//							if(pstat[h*htls_time_comp*2+k].fftAccum != 168)
//								printf("fft:%d\n",pstat[h*htls_time_comp*2+k].fftAccum);
							spolset_lths_off[h1].fft_weight += (int)pstat[h*lths_time_comp*2+k].fftAccum;
							for(l=0;l<RAW_CHANNELS;l++)
							{
								spolset_lths_off[h1].A[l] += (unsigned int)(pdatum[h*lths_time_comp*2+k].pol_A[l]);
								spolset_lths_off[h1].B[l] += (unsigned int)(pdatum[h*lths_time_comp*2+k].pol_B[l]);
								spolset_lths_off[h1].U[l] += (int)(pdatum[h*lths_time_comp*2+k].stokes_U[l]);
								spolset_lths_off[h1].V[l] += (int)(pdatum[h*lths_time_comp*2+k].stokes_V[l]);
							}//l loop htls	
							num_off++;
						}//onoff if
						else
						{
//							if(pstat[h*htls_time_comp*2+k].fftAccum != 168)
//								printf("fft:%d\n",pstat[h*htls_time_comp*2+k].fftAccum);
							spolset_lths_off[h1].fft_weight += (int)pstat[h*lths_time_comp*2+k].fftAccum;
							for(l=0;l<RAW_CHANNELS;l++)
							{
								spolset_lths_on[h2].A[l] += (unsigned int)(pdatum[h*lths_time_comp*2+k].pol_A[l]);
								spolset_lths_on[h2].B[l] += (unsigned int)(pdatum[h*lths_time_comp*2+k].pol_B[l]);
								spolset_lths_on[h2].U[l] += (int)(pdatum[h*lths_time_comp*2+k].stokes_U[l]);
								spolset_lths_on[h2].V[l] += (int)(pdatum[h*lths_time_comp*2+k].stokes_V[l]);
							}//l loop htls	
							num_on++;
						}//onoff else
					}
					else*/
//					if(pstat[h*lths_time_comp*2+k].calOn == 0)
					if(!onoff)					
					{
//						flag = 0;
//						int l;
						if(h1==DUMPS_PER_ROW/(lths_time_comp*2))
							printf("Problematic h1 %d\n",h1);
						if(offcount < 9)
						{
							spolset_lths_off[h1].fft_weight += (int)pstat[h*lths_time_comp*2+k].fftAccum;
							for(l=0;l<RAW_CHANNELS;l++)
							{
								spolset_lths_off[h1].A[l] += pdatum[h*lths_time_comp*2+k].pol_A[l];
								spolset_lths_off[h1].B[l] += pdatum[h*lths_time_comp*2+k].pol_B[l];
								spolset_lths_off[h1].U[l] += pdatum[h*lths_time_comp*2+k].stokes_U[l];
								spolset_lths_off[h1].V[l] += pdatum[h*lths_time_comp*2+k].stokes_V[l];
							}//l loop htls	
						}
						num_off++;
	                                        offcount++;     
                                                if(offcount == 10)
                                                {                                                                                                                                    onoff = 1;
                                                        offcount = 0;                                                                                                        }

					}//fi
/*					else if(pstat[h*lths_time_comp*2+k].calOn == 1 && flag == 0)
					{
						flag = 1;
						num_off++;
					}
					else if(pstat[h*lths_time_comp*2+k].calOn == 1  && pstat[h*lths_time_comp*2+k+1].calOn == 0)
					{
						flag = 1;
						num_on++;
					}
*/					else
					{
//						int l;
						if(h2==DUMPS_PER_ROW/(lths_time_comp*2))
							printf("Problematic h2 %d\n",h2);
						if(oncount < 9)
						{
							spolset_lths_on[h2].fft_weight += (int)pstat[h*lths_time_comp*2+k].fftAccum;
							for(l=0;l<RAW_CHANNELS;l++)
							{
								spolset_lths_on[h2].A[l] += pdatum[h*lths_time_comp*2+k].pol_A[l];
								spolset_lths_on[h2].B[l] += pdatum[h*lths_time_comp*2+k].pol_B[l];
								spolset_lths_on[h2].U[l] += pdatum[h*lths_time_comp*2+k].stokes_U[l];
								spolset_lths_on[h2].V[l] += pdatum[h*lths_time_comp*2+k].stokes_V[l];
							}//l loop htls	
						}
						num_on++;
                                                oncount++; 
                                                if(oncount == 10)
                                                {                 
                                                        onoff = 0;
                                                        oncount = 0;                        
                                                }

//						flag = 1;
					}//final else
//					int l;
					if(num_on == lths_time_comp)
					{
//						for(l=0;l<RAW_CHANNELS;l++)
//						{
//							spolset_lths_on[h2].A[l] /= spolset_lths_on[h2].fft_weight;
//							spolset_lths_on[h2].B[l] /= spolset_lths_on[h2].fft_weight;
//							spolset_lths_on[h2].U[l] /= spolset_lths_on[h2].fft_weight;
//							spolset_lths_on[h2].V[l] /= spolset_lths_on[h2].fft_weight;
//						}
						num_on = 0;
						h2++;
//						onoff = 0;
					}
					if(num_off == lths_time_comp)
					{
//						for(l=0;l<RAW_CHANNELS;l++)
//						{
//							spolset_lths_off[h1].A[l] /= spolset_lths_off[h1].fft_weight;
//							spolset_lths_off[h1].B[l] /= spolset_lths_off[h1].fft_weight;
//							spolset_lths_off[h1].U[l] /= spolset_lths_off[h1].fft_weight;
//							spolset_lths_off[h1].V[l] /= spolset_lths_off[h1].fft_weight;
//						}
						num_off = 0;
						h1++;
//						onoff = 1;
					}
				}//k loop lths
				
			}//h loop for lths

			//printf("DIAGNOSTIC:h1%d,h2%d",h1,h2);

			int m,n,p;
			for(p=0;p<DUMPS_PER_ROW/(2*htls_time_comp);p++)
			{
				for(m = 0;m < RAW_CHANNELS/htls_spec_comp;m++)
				{
					XXoff_htls[m] = 0;
					YYoff_htls[m] = 0;
					XYoff_htls[m] = 0;
					YXoff_htls[m] = 0;
					XXon_htls[m] = 0;
					YYon_htls[m] = 0;
					XYon_htls[m] = 0;
					YXon_htls[m] = 0;
				}					
				for(m = 0;m < RAW_CHANNELS/htls_spec_comp;m++)
				{
					for(n=0;n<htls_spec_comp;n++)
					{
						spolset_htls_off_final[p].A[m] += spolset_htls_off[p].A[m*htls_spec_comp+n];
						spolset_htls_off_final[p].B[m] += spolset_htls_off[p].B[m*htls_spec_comp+n];
						spolset_htls_off_final[p].U[m] += spolset_htls_off[p].U[m*htls_spec_comp+n];
						spolset_htls_off_final[p].V[m] += spolset_htls_off[p].V[m*htls_spec_comp+n];
						spolset_htls_on_final[p].A[m] += spolset_htls_on[p].A[m*htls_spec_comp+n];
						spolset_htls_on_final[p].B[m] += spolset_htls_on[p].B[m*htls_spec_comp+n];
						spolset_htls_on_final[p].U[m] += spolset_htls_on[p].U[m*htls_spec_comp+n];
						spolset_htls_on_final[p].V[m] += spolset_htls_on[p].V[m*htls_spec_comp+n];
					}//n loop
//					spolset_htls_off_final[p].A[m] /= htls_spec_comp;
//					spolset_htls_off_final[p].B[m] /= htls_spec_comp;
//					spolset_htls_off_final[p].U[m] /= htls_spec_comp;	
//					spolset_htls_off_final[p].V[m] /= htls_spec_comp;
//					spolset_htls_on_final[p].A[m] /= htls_spec_comp;
//					spolset_htls_on_final[p].B[m] /= htls_spec_comp;
//					spolset_htls_on_final[p].U[m] /= htls_spec_comp;	
//					spolset_htls_on_final[p].V[m] /= htls_spec_comp;
					XXoff_htls[m] = ((float)spolset_htls_off_final[p].A[m]/((float)spolset_htls_off[p].fft_weight*2));				
					YYoff_htls[m] = ((float)spolset_htls_off_final[p].B[m]/((float)spolset_htls_off[p].fft_weight*2));
					XYoff_htls[m] = (((float)spolset_htls_off_final[p].U[m]+(float)spolset_htls_off_final[p].V[m])/((float)spolset_htls_off[p].fft_weight*2));
					YXoff_htls[m] = (((float)spolset_htls_off_final[p].U[m]-(float)spolset_htls_off_final[p].V[m])/((float)spolset_htls_off[p].fft_weight*2));
					XXon_htls[m] = ((float)spolset_htls_on_final[p].A[m]/((float)spolset_htls_on[p].fft_weight*2));				
					YYon_htls[m] = ((float)spolset_htls_on_final[p].B[m]/((float)spolset_htls_on[p].fft_weight*2));
					XYon_htls[m] = (((float)spolset_htls_on_final[p].U[m]+(float)spolset_htls_on_final[p].V[m])/((float)spolset_htls_on[p].fft_weight*2));
					YXon_htls[m] = (((float)spolset_htls_on_final[p].U[m]-(float)spolset_htls_on_final[p].V[m])/((float)spolset_htls_on[p].fft_weight*2));										
//					printf("A: %u B %u U %d V %d W %d\n",spolset_htls_on_final[p].A[m],spolset_htls_on_final[p].B[m],spolset_htls_on_final[p].U[m],spolset_htls_on_final[p].V[m],spolset_htls_on[p].fft_weight);
//					printf("XX: %f YY %f XY %f YX %f\n",XXon_htls[m],YYon_htls[m],XYon_htls[m],YXon_htls[m]);
				}//m loop

//				fwrite(&spointing_htls[p],sizeof(nSpecPointingBlock),1,htls_file);
				fwrite(&spointing_htls[p],sizeof(SpecPointingBlock),1,htls_file);
				fwrite(XXon_htls,sizeof(float),RAW_CHANNELS/htls_spec_comp,htls_file);
				fwrite(YYon_htls,sizeof(float),RAW_CHANNELS/htls_spec_comp,htls_file);			
				fwrite(XYon_htls,sizeof(float),RAW_CHANNELS/htls_spec_comp,htls_file);			
				fwrite(YXon_htls,sizeof(float),RAW_CHANNELS/htls_spec_comp,htls_file);			
				fwrite(XXoff_htls,sizeof(float),RAW_CHANNELS/htls_spec_comp,htls_file);
				fwrite(YYoff_htls,sizeof(float),RAW_CHANNELS/htls_spec_comp,htls_file);			
				fwrite(XYoff_htls,sizeof(float),RAW_CHANNELS/htls_spec_comp,htls_file);			
				fwrite(YXoff_htls,sizeof(float),RAW_CHANNELS/htls_spec_comp,htls_file);
								
								
/*				fwrite(spolset_htls_on_final[p].A,size_int,RAW_CHANNELS/htls_spec_comp,htls_file);				
				fwrite(spolset_htls_on_final[p].B,size_int,RAW_CHANNELS/htls_spec_comp,htls_file);				
				fwrite(spolset_htls_on_final[p].U,size_int,RAW_CHANNELS/htls_spec_comp,htls_file);				
				fwrite(spolset_htls_on_final[p].V,size_int,RAW_CHANNELS/htls_spec_comp,htls_file);				
				fwrite(&spolset_htls_on_final[p].fft_weight,size_int,1,htls_file);				
				fwrite(spolset_htls_off_final[p].A,size_int,RAW_CHANNELS/htls_spec_comp,htls_file);				
				fwrite(spolset_htls_off_final[p].B,size_int,RAW_CHANNELS/htls_spec_comp,htls_file);				
				fwrite(spolset_htls_off_final[p].U,size_int,RAW_CHANNELS/htls_spec_comp,htls_file);				
				fwrite(spolset_htls_off_final[p].V,size_int,RAW_CHANNELS/htls_spec_comp,htls_file);				
				fwrite(&spolset_htls_off_final[p].fft_weight,size_int,1,htls_file);				
*/			}// p loop
		

			for(p=0;p<DUMPS_PER_ROW/(2*lths_time_comp);p++)
			{
				for(m = 0;m < RAW_CHANNELS/lths_spec_comp;m++)
				{
					XXoff_lths[m] = 0;
					YYoff_lths[m] = 0;
					XYoff_lths[m] = 0;
					YXoff_lths[m] = 0;
					XXon_lths[m] = 0;
					YYon_lths[m] = 0;
					XYon_lths[m] = 0;
					YXon_lths[m] = 0;
				}	
								
				for(m = 0;m < RAW_CHANNELS/lths_spec_comp;m++)
				{
					for(n=0;n<lths_spec_comp;n++)
					{
						spolset_lths_off_final[p].A[m] += spolset_lths_off[p].A[m*lths_spec_comp+n];
						spolset_lths_off_final[p].B[m] += spolset_lths_off[p].B[m*lths_spec_comp+n];
						spolset_lths_off_final[p].U[m] += spolset_lths_off[p].U[m*lths_spec_comp+n];
						spolset_lths_off_final[p].V[m] += spolset_lths_off[p].V[m*lths_spec_comp+n];
						spolset_lths_on_final[p].A[m] += spolset_lths_on[p].A[m*lths_spec_comp+n];
						spolset_lths_on_final[p].B[m] += spolset_lths_on[p].B[m*lths_spec_comp+n];
						spolset_lths_on_final[p].U[m] += spolset_lths_on[p].U[m*lths_spec_comp+n];
						spolset_lths_on_final[p].V[m] += spolset_lths_on[p].V[m*lths_spec_comp+n];
					}//n loop
				
//					spolset_lths_off_final[p].A[m] /= lths_spec_comp;
//					spolset_lths_off_final[p].B[m] /= lths_spec_comp;
//					spolset_lths_off_final[p].U[m] /= lths_spec_comp;	
//					spolset_lths_off_final[p].V[m] /= lths_spec_comp;
//					spolset_lths_on_final[p].A[m] /= lths_spec_comp;
//					spolset_lths_on_final[p].B[m] /= lths_spec_comp;
//					spolset_lths_on_final[p].U[m] /= lths_spec_comp;	
//					spolset_lths_on_final[p].V[m] /= lths_spec_comp;

					XXoff_lths[m] = ((float)spolset_lths_off_final[p].A[m]/((float)spolset_lths_off[p].fft_weight*2));				
					YYoff_lths[m] = ((float)spolset_lths_off_final[p].B[m]/((float)spolset_lths_off[p].fft_weight*2));
					XYoff_lths[m] = (((float)spolset_lths_off_final[p].U[m]+(float)spolset_lths_off_final[p].V[m])/((float)spolset_lths_off[p].fft_weight*2));
					YXoff_lths[m] = (((float)spolset_lths_off_final[p].U[m]-(float)spolset_lths_off_final[p].V[m])/((float)spolset_lths_off[p].fft_weight*2));
					XXon_lths[m] = ((float)spolset_lths_on_final[p].A[m]/((float)spolset_lths_on[p].fft_weight*2));				
					YYon_lths[m] = ((float)spolset_lths_on_final[p].B[m]/((float)spolset_lths_on[p].fft_weight*2));
					XYon_lths[m] = (((float)spolset_lths_on_final[p].U[m]+(float)spolset_lths_on_final[p].V[m])/((float)spolset_lths_on[p].fft_weight*2));
					YXon_lths[m] = (((float)spolset_lths_on_final[p].U[m]-(float)spolset_lths_on_final[p].V[m])/((float)spolset_lths_on[p].fft_weight*2));
				}//m loop



//				fwrite(&spointing_lths[p],sizeof(nSpecPointingBlock),1,lths_file);				
				fwrite(&spointing_lths[p],sizeof(SpecPointingBlock),1,lths_file);	
				fwrite(XXon_lths,sizeof(float),RAW_CHANNELS/lths_spec_comp,lths_file);
				fwrite(YYon_lths,sizeof(float),RAW_CHANNELS/lths_spec_comp,lths_file);			
				fwrite(XYon_lths,sizeof(float),RAW_CHANNELS/lths_spec_comp,lths_file);			
				fwrite(YXon_lths,sizeof(float),RAW_CHANNELS/lths_spec_comp,lths_file);			
				fwrite(XXoff_lths,sizeof(float),RAW_CHANNELS/lths_spec_comp,lths_file);
				fwrite(YYoff_lths,sizeof(float),RAW_CHANNELS/lths_spec_comp,lths_file);			
				fwrite(XYoff_lths,sizeof(float),RAW_CHANNELS/lths_spec_comp,lths_file);			
				fwrite(YXoff_lths,sizeof(float),RAW_CHANNELS/lths_spec_comp,lths_file);
/*				fwrite(spolset_lths_on_final[p].A,size_int,RAW_CHANNELS/lths_spec_comp,lths_file);				
				fwrite(spolset_lths_on_final[p].B,size_int,RAW_CHANNELS/lths_spec_comp,lths_file);				
				fwrite(spolset_lths_on_final[p].U,size_int,RAW_CHANNELS/lths_spec_comp,lths_file);				
				fwrite(spolset_lths_on_final[p].V,size_int,RAW_CHANNELS/lths_spec_comp,lths_file);				
				fwrite(&spolset_lths_on_final[p].fft_weight,size_int,1,lths_file);				
				fwrite(spolset_lths_off_final[p].A,size_int,RAW_CHANNELS/lths_spec_comp,lths_file);				
				fwrite(spolset_lths_off_final[p].B,size_int,RAW_CHANNELS/lths_spec_comp,lths_file);				
				fwrite(spolset_lths_off_final[p].U,size_int,RAW_CHANNELS/lths_spec_comp,lths_file);				
				fwrite(spolset_lths_off_final[p].V,size_int,RAW_CHANNELS/lths_spec_comp,lths_file);
				fwrite(&spolset_lths_off_final[p].fft_weight,size_int,1,lths_file);				
*/			}// p loop
	
		}//naxis2 loop

		fclose(datafile);
	}//num_files for loop


//	free(&spointing_htls);
//	free(&spointing_lths);

//	free(&spolset_htls_on);
//	free(&spolset_htls_off);
//	free(&spolset_lths_on);
//	free(&spolset_lths_off);
//	free(&spolset_htls_on_final);
//	free(&spolset_htls_off_final);
//	free(&spolset_lths_on_final);
//	free(&spolset_lths_off_final);
	fclose(htls_file);
	fclose(lths_file);
	return 1;
}

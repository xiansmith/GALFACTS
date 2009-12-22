#include <stdio.h>
#include <string.h>
//Phil's header file defining data-types for the decimated files
#include "galfactsLib.h"
//Data structures for spec file format (Defined by Desh)
#include "spec.h"
//Header file for this program
#include "dfrq2spec.h"
//Jeff Dever's file utility
#include "jsd/jsd_futil.h"

/*each dfrq file contains 109 rows (each row has data worth 600ms) except
 * for the last file on a given day.
*/
#define MAX_ROWS 109


//following are routines to convert endian-ness for various data-types
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
	char *proj_code; //a2130
	char *date; //date of observation
	char *datadir; //path where the data-files reside

	/* A File containing the list of decimated file names for a given day of
	 * observations that need to be combined and converted into a single spec
	 *  file. Typically there are hundred's of files for each day of
	 *  observation. I was originally using number of files argument
	 *  instead of filelist,but Roberto asked for filelist provision instead.
	 *  He can produce a filelist using ls command. This removes difficulty
	 *  with file name numbering starting from a different number each day due to
	 *  cima quirks. The filelist can be produces for a given day by issuing the
	 *  following command:
	 *  ls -1 *dfrq* > [filelist filename]
	 *  assuming a folder contains all files from a single band for a given
	 *  day of observations.
	 */
	char *filelistname;

	int band; //band number 0 or 1
	int config_not_written=1;

	if(argc !=7)
	{
		printf("usage: dfrq2spec <proj_code> <date> <band> <beam> <datadir> <filelist>\n");
		return 0;
	}
	else
	{
		proj_code = argv[1];
		date = argv[2];
		band = atoi(argv[3]);
		beam = atoi(argv[4]);
		datadir = argv[5];
 		filelistname = argv[6];
	}

	FILE *datafile,*specfile,*cfg_file,*listfile;
	char datafilename[100+1],specfilename[40+1],cfgfilename[40+1];
	listfile = fopen(filelistname,"r");
	num_files = jsd_line_count(listfile);
	printf("Number of files to be processed:%d\n",num_files);

	sprintf(specfilename,"%s.%s.b%1ds%1d.htls.spec",proj_code,date,beam,band);
	if ( (specfile = fopen(specfilename, "wb") ) == NULL )
	{
		printf("ERROR: can't open data file for writing '%s'\n", specfilename);
		return 0;
	}

	char buf[LINELEN+1];
	long int offset = MAIN_HEADER + BINTABLE_HEADER;
	int f,g,k,l,found = FALSE;
	//read in files one by one and process them
	for(f=0;f<num_files;f++)
	{
		int naxis1,naxis2;
		GFLIB_ROW row1,row2;
		SpecPointingBlock SPBlock;
		//Warning: Make sure MAX_CHANNELS is set to 128 in common.h
		float A[MAX_CHANNELS],B[MAX_CHANNELS];
		float XX[MAX_CHANNELS],YY[MAX_CHANNELS];
		float U[MAX_CHANNELS],V[MAX_CHANNELS];
		float XY[MAX_CHANNELS],YX[MAX_CHANNELS];

		//Use filelist to know which files to process
		fscanf(listfile,"%s",datafilename);
		if ( (datafile = fopen(datafilename, "r") ) == NULL )
		{
			printf("ERROR: can't open data file for reading '%s'\n", datafilename);
			return 0;
		}
		printf("Opened the datafile:%s\n",datafilename);

		//Read in NAXIS1 from fits header
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

		//Read in NAXIS2 from fits header
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

		//Process each file row by row
		for(g=0;g < naxis2-1;g++)
		{
			printf("Reading row :%d\n",g+1);
			fseek(datafile,offset+g*naxis1,SEEK_SET);

			/*Read 2 rows of data at a time, there is one time stamp for each row
			 * therefore need to read 2 rows so as to be able to interpolate time
			 * stamps. Once problem with this currently is that the last row in each file
			 * is not process and that data is lost. This is only 1/109th of the data though
			 * or 600ms per 65.4 seconds
			 */
			fread(&row1,sizeof(GFLIB_ROW),1,datafile);
			fread(&row2,sizeof(GFLIB_ROW),1,datafile);

			// convert endian-ness of data
			cnvrt_end_db(&row1.RA);
			cnvrt_end_db(&row2.RA);
			cnvrt_end_db(&row1.DEC);
			cnvrt_end_db(&row2.DEC);
			cnvrt_end_db(&row1.UTC);
			cnvrt_end_db(&row2.UTC);
			cnvrt_end_db(&row1.mjdxxobs);
			cnvrt_end_db(&row1.azimuth);
			cnvrt_end_db(&row1.elevatio);
			cnvrt_end_db(&row2.azimuth);
			cnvrt_end_db(&row2.elevatio);
			cnvrt_end_db(&row1.req_raj);
			cnvrt_end_db(&row1.req_decj);
			cnvrt_end_db(&row2.req_raj);
			cnvrt_end_db(&row2.req_decj);
			cnvrt_end_db(&row1.alfa_ang);

			//Write the spec_cfg files during the first pass
			if(config_not_written)
			{
				sprintf(cfgfilename,"%s.%s.b%ds%d.htls.spec_cfg",proj_code,date,beam,band);
				if ( (cfg_file = fopen(cfgfilename, "w") ) == NULL )
				{
					printf("ERROR: can't open config file '%s'\n", cfgfilename);
					return 0;
				}

				cnvrt_end_db(&row1.tdelt);
				fprintf(cfg_file,"%f\n",row1.tdelt*1000);
				fprintf(cfg_file,"%i\n",RAW_CHANNELS);
				cnvrt_end_db(&row1.cf);
				fprintf(cfg_file,"%f\n",row1.cf/1000000);
				cnvrt_end_db(&row1.fdelt);
				fprintf(cfg_file,"%f\n",-1*row1.fdelt*RAW_CHANNELS/1000);
				fprintf(cfg_file,"1 %i 4 2 0\n",RAW_CHANNELS);
				fprintf(cfg_file,"%s\n",proj_code);
				fprintf(cfg_file,"%f\n",row1.mjdxxobs);
				fprintf(cfg_file,"AO\n");
				fprintf(cfg_file,"Integration time (ms):%f\n",row1.tdelt*1000);
				fprintf(cfg_file,"MJD: %f\n",row1.mjdxxobs);
				fprintf(cfg_file,"Center freq (MHz): %f\n",row1.cf/1000000);
				fprintf(cfg_file,"Channel band (kHz): %f\n",-1*row1.fdelt/1000);
				fprintf(cfg_file,"Number of channels/record: %d\n",RAW_CHANNELS);
				fprintf(cfg_file,"RA at start (degrees): %f\n",row1.RA);
				fprintf(cfg_file,"DEC at start (degrees): %f\n",row1.DEC);
				fprintf(cfg_file,"UTC at start (seconds): %f\n",row1.UTC);
				fprintf(cfg_file,"ALFA angle (degrees) at start: %f\n",row1.alfa_ang);
				fprintf(cfg_file,"Project ID: %s\n",proj_code);
				fclose(cfg_file);
				config_not_written = 0;
			}

			for(k=0;k<DUMPS_PER_ROW;k++)
			{
				//Interpolate time stamps
				if(!beam)
				{
					SPBlock.centralBeam.raj_true_in_hours=row1.RA/15 + k*(row2.RA-row1.RA)/(DUMPS_PER_ROW)/15;
					SPBlock.centralBeam.decj_true_in_degrees=row1.DEC + k*(row2.DEC-row1.DEC)/(DUMPS_PER_ROW);
					SPBlock.centralBeam.atlantic_solar_time_now_in_sec=\
					row1.UTC + k*(row2.UTC-row1.UTC)/(DUMPS_PER_ROW);
				}
				else
				{
					SPBlock.outerBeams[beam-1].raj_true_in_hours=row1.RA/15 + k*(row2.RA-row1.RA)/(DUMPS_PER_ROW)/15;
					SPBlock.outerBeams[beam-1].decj_true_in_degrees=row1.DEC + k*(row2.DEC-row1.DEC)/(DUMPS_PER_ROW);
					SPBlock.centralBeam.atlantic_solar_time_now_in_sec=\
					row1.UTC + k*(row2.UTC-row1.UTC)/(DUMPS_PER_ROW);
				}
				cnvrt_end_sint(&row1.stat[k].fftAccum);
			//	cnvrt_end_sint(&row1.statoff[k].fftAccum);
				fwrite(&SPBlock,sizeof(SpecPointingBlock),1,specfile);

				/*Process data converting A,B,U,V into XX,XY,YY,YX values and normalizing with
				 * respect to fft accumulations
				 */
				for(l=0;l<RAW_CHANNELS;l++)
				{
					//byteswap
					cnvrt_end_int(&row1.data[k].A[l]);
					cnvrt_end_int(&row1.data[k].B[l]);
					cnvrt_end_int(&row1.data[k].U[l]);
					cnvrt_end_int(&row1.data[k].V[l]);
			//		cnvrt_end_int(&row1.dataoff[k].A[l]);
			//		cnvrt_end_int(&row1.dataoff[k].B[l]);
			//		cnvrt_end_int(&row1.dataoff[k].U[l]);
			//		cnvrt_end_int(&row1.dataoff[k].V[l]);
					//normalize
					A[l] = (float)(row1.data[k].A[l]/row1.stat[k].fftAccum);
					B[l] = (float)(row1.data[k].B[l]/row1.stat[k].fftAccum);
					U[l] = (float)((int)row1.data[k].U[l]/row1.stat[k].fftAccum);
					V[l] = (float)((int)row1.data[k].V[l]/row1.stat[k].fftAccum);
			//		Aoff[l] = (float)(row1.dataoff[k].A[l]/row1.statoff[k].fftAccum);
			//		Boff[l] = (float)(row1.dataoff[k].B[l]/row1.statoff[k].fftAccum);
			//		Uoff[l] = (float)((int)row1.dataoff[k].U[l]/row1.statoff[k].fftAccum);
			//		Voff[l] = (float)((int)row1.dataoff[k].V[l]/row1.statoff[k].fftAccum);
					//convert to xx yy xy yx
					XX[l] = A[l]/2;
			//		XXoff[l] = Aoff[l]/2;
					YY[l] = B[l]/2;
			//		YYoff[l] = Boff[l]/2;
					XY[l] = (U[l]+V[l])/2;
			//		XYoff[l] = (Uoff[l]+Voff[l])/2;
					YX[l] = (U[l]-V[l])/2;
			//		YXoff[l] = (Uoff[l]-Voff[l])/2;
//					XXon[l] = Aon[l];
//					XXoff[l] = Aoff[l];
//					YYon[l] = Bon[l];
//					YYoff[l] = Boff[l];
//					XYon[l] = Uon[l];
//					XYoff[l] = Uoff[l];
//					YXon[l] = Von[l];
//					YXoff[l] = Voff[l];
//					printf("%lf %lf %lf %lf\n",XX[l],YY[l],XY[l],YX[l]);

				}//l loop for each channel

				//Write out the data to spec file
				fwrite(&XX,sizeof(float),MAX_CHANNELS,specfile);
				fwrite(&YY,sizeof(float),MAX_CHANNELS,specfile);
				fwrite(&XY,sizeof(float),MAX_CHANNELS,specfile);
				fwrite(&YX,sizeof(float),MAX_CHANNELS,specfile);
				//fwrite(&XXoff,sizeof(float),MAX_CHANNELS,specfile);
				//fwrite(&YYoff,sizeof(float),MAX_CHANNELS,specfile);
				//fwrite(&XYoff,sizeof(float),MAX_CHANNELS,specfile);
				//fwrite(&YXoff,sizeof(float),MAX_CHANNELS,specfile);
			}//k loop num dumps
			cnvrt_end_db(&row2.RA);
			cnvrt_end_db(&row2.DEC);
			cnvrt_end_db(&row2.UTC);

		}//naxis2 loop g
		fclose(datafile);
	}//num files loop f
	fclose(specfile);
	return 1;
}

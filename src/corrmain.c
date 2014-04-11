#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "programs/fitsio.h"
#include "fluxdata.h"
#include "beammodels.h"
#include "jsd/jsd_futil.h"
#include "map.h"
#include "grid.h"
#include "cfit/chebyshev.h"
//#include "clean.h"
int multibeam;
static void print_usage(const char * prog)
{
	printf(
	"Usage: %s <fitsmap> <lowchan> <highchan> "
	"<grid patch radius> <clean patch radius> <clean loopgain> <stop rms> <title> <show progress>\n"
	"\n"
	"\tfitsmap - the FITS map to clean\n"
	"\tlowchan, highchan - the channel range to use\n"
	"\tgrid patch radius - area in pixels to calculate point spread response\n"
	"\tclean patch radius - area in pixels to subtract dirty beam\n"
	"\tclean loopgain - the gain of each clean iteration (should be between 0 and 1)\n"
	"\tstop rms - the rms threshold where the clean iteraton will stop\n"
	"\ttitle - the title for the clean map\n"
	"\tshow progress - if true will create a fits file showing the way the iterations progressed\n"
	"\n"
	"eg: %s mutlibeam_Icube.fits 30 230 8 .1333 0.2 35.5 \"GALFACTS Test Calibration region120 MULTIBEAM \
	CLEANED\" \n"
	"\n", prog, prog); 
}

int main(int argc, char * argv[])
{
	FluxWappData *wappdata;
	MapMetaData md;
	char ** files;
	//int numDays,numDays1,i,j,k,l,m,showprogress,progress_counter = 0,linecount,beammodellinecount;
	int numDays,i,j,k,l,m,beammodellinecount;
//	int beam;//SSG
	char * dirtymapname;
	//char beammodelfilename[64],maptimedomainfilename[64],tempheader[80],copyfilename[64],cleancomponentfilename[64];
	char beammodelfilename[64],cleancomponentfilename[64];
	float mappeakRA[MAX_CHANNELS],mappeakDEC[MAX_CHANNELS],mapmaxpower[MAX_CHANNELS],mapminpower[MAX_CHANNELS];
	//float loop_gain;
	//double stop_rms,rms;
	float clean_patch;
	char *wapp;
	//FILE * dirtymap_file,* progressfileI, *beammodelfile,*maptimedomainfile, *copyfile, *cleancomponentfile, *annfile;
	FILE * dirtymap_file, *beammodelfile, *cleancomponentfile;
//	FILE * progressfileQ ,* progressfileU ,* progressfileV;
	//header_param_list dirtymap_hpar,progressfile_hpar;
	header_param_list dirtymap_hpar;
	//multibeam = 1;
//	int data_len;
//	float * plane_data;
	
	/* Process command line arguments */ 
	/* Convert args into more usable units were required */

	if (argc != 8) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	} else {
		wapp = argv[1];
		dirtymapname = argv[2];
		md.lowchan = atoi(argv[3]);
		md.highchan = atoi(argv[4]);
		md.patch = atoi(argv[5]);
		clean_patch = (float )atof(argv[6]);
//		loop_gain = (float) atof(argv[6]);
//		stop_rms = (double) atof(argv[7]); 
		md.title = argv[7];
//		showprogress = atoi(argv[9]);
		dirtymap_file = fopen(dirtymapname, "r");
		readfits_header(dirtymap_file, &dirtymap_hpar);
		md.beamwidths = md.patch;
//		md.fcen = (float) atof(argv[2]);
		md.cellsize = fabs(dirtymap_hpar.cdelt[0]); 
		md.RAcen = dirtymap_hpar.crval[0];
		md.DECcen = dirtymap_hpar.crval[1];
		md.RArange =  fabs(dirtymap_hpar.naxis[0]*dirtymap_hpar.cdelt[0]);
		md.DECrange = fabs(dirtymap_hpar.naxis[1]*dirtymap_hpar.cdelt[1]);
		md.n1 = dirtymap_hpar.naxis[0];
		md.n2 = dirtymap_hpar.naxis[1];
		md.n3 = dirtymap_hpar.naxis[2];
		md.fwhm = 2.0;
		md.fstart = dirtymap_hpar.crval[2];
		md.ramin = md.RAcen - (md.n1*md.cellsize)/2;
		md.ramax = md.RAcen + (md.n1*md.cellsize)/2;
		md.decmin = md.DECcen - (md.n2*md.cellsize)/2;
		md.decmax = md.DECcen + (md.n2*md.cellsize)/2;
		fclose(dirtymap_file);
	printf("Dirty map name: %s\n",dirtymapname);
	printf("Map size: %d x %d\n", md.n1, md.n2);
	printf("DEC center: %f\n", md.DECcen);
	printf("DEC range: %f ... %f\n", md.decmin, md.decmax);
	printf("RA range: %f ... %f\n", md.ramin, md.ramax);
	printf("Map centre: %.2f %.2f\n", md.RAcen, md.DECcen); 
	printf("Cell size: %7.5f\n", md.cellsize);
	printf("Channel range: (%i, %i]\n", md.lowchan, md.highchan);
	printf("Channel count: %i\n", md.n3);
	printf("Grid Patch radius: %i\n", md.patch);
	printf("Clean Patch radius: %f\n", clean_patch);
//	printf("Center Frequency: %g\n", md->fcen); 
	printf("Start Frequency: %g\n", md.fstart);
// 	printf("Balance Order: %i\n", balorder);
// 	printf("Clean Loop Gain: %f\n", loop_gain);
 //	printf("Stop rms: %f\n", stop_rms);

	}

	//numDays = get_date_dirs("./", &files);
	numDays = get_date_dirs("/n/fox/processed/FIELD1/run6/band0", &files);
	printf("DIAGNOSTIC: No of date dirs %d\n",numDays);
	if (numDays <= 0) {
		printf("ERROR: could not find any date dirs\n");
		return EXIT_FAILURE;
	}
        if (!strcmp(wapp,"multibeam"))
        {
                numDays = numDays * 7;
                multibeam = 1;
        }
        else
                multibeam = 0;

//	FluxWappData *wappdata;
        wappdata = fluxwappdata_alloc(wapp, files, numDays);





        init_psf_map(md.fwhm, md.cellsize, md.beamwidths);
	init_psf_lookup_table(65537, 9.1);
	
/*	if(showprogress)
	{
		progressfileI = fopen("CLEAN_ProgressI.fits", "w");
		progressfileQ = fopen("CLEAN_ProgressQ.fits", "w");
		progressfileU = fopen("CLEAN_ProgressU.fits", "w");
		progressfileV = fopen("CLEAN_ProgressV.fits", "w");
		progressfile_hpar = dirtymap_hpar;
		progressfile_hpar.naxis[2] = 0;
		writefits_header(progressfileI, &progressfile_hpar);
		writefits_header(progressfileQ, &progressfile_hpar);
		writefits_header(progressfileU, &progressfile_hpar);
		writefits_header(progressfileV, &progressfile_hpar);
	}
*/	
//	annfile = fopen("CLEAN_ann.ann","w");
//	if(annfile == NULL)
//	{
//		printf("Error : Could not open clean annotation file.\n");
//		exit(EXIT_FAILURE);	
//	}
//	fprintf(annfile, "COLOUR RED\n");




	cleancomponentfile = fopen("N1cc.dat","r");

	if(cleancomponentfile == NULL)
	{
		printf("Error : Could not open clean component file.\n");
		exit(EXIT_FAILURE);
	}

	int cclinecount = jsd_line_count(cleancomponentfile);
	float *ccStokesI = (float*)malloc(cclinecount * sizeof(float));
	float *ccRA = (float*)malloc(cclinecount * sizeof(float));
	float *ccDEC = (float*)malloc(cclinecount * sizeof(float));
	for(i = 0;i< cclinecount;i++)
		fscanf(cleancomponentfile,"%f %f %f",&ccRA[i],&ccDEC[i],&ccStokesI[i]);

	fclose(cleancomponentfile);
	printf("DIAGNOSTIC: Found %d Clean Components\n",cclinecount);



			float *dataI, *dataQ, *dataU, *dataV, *weight;
			int numbytes; 
			numbytes = md.n1 * md.n2 * sizeof (float);
			dataI  = (float *) malloc (numbytes);
			dataQ  = (float *) malloc (numbytes);
			dataU  = (float *) malloc (numbytes);
			dataV  = (float *) malloc (numbytes);
			weight  = (float *) malloc (numbytes);
	
			if (!dataI || !dataQ || !dataU || !dataV || !weight) {
				printf("ERROR: memory allocation of %i bytes failed!\n", numbytes);
				return EXIT_FAILURE;
			}

	start_fits_cubes(wapp,&md);



	do
	{
		//get_peak_power_coord(dirtymapname,mappeakRA,mappeakDEC,mapmaxpower,mapminpower,md.lowchan,md.highchan);
//		printf("DIAGNOSTIC: Multibeam CLEAN iteration %d.\n",progress_counter+1);
//		rms = 0.0;
//		long int rms_counter = 0;
		for(i = md.lowchan;i < md.highchan;i++)
		{
 i = 0;
                        fluxwappdata_readchan(wappdata, i, CLEAN);



//			printf("DIAGNOSTIC: Channel %d, Peak Flux %f at RA %f DEC %f\n",i,mapmaxpower[i],mappeakRA[i],mappeakDEC[i]);
//			sprintf(cleancomponentfilename,"cleancomponent%03i.dat",i);
//			cleancomponentfile = fopen(cleancomponentfilename,"w");
//			if(!progress_counter)
//			{
//				cleancomponentfile = fopen(cleancomponentfilename,"w");
//				fprintf(cleancomponentfile,"RA DEC FLUX_I\n");
//			}
//			else
//				cleancomponentfile = fopen(cleancomponentfilename,"a+");
//			if(cleancomponentfile == NULL)
//			{
//				printf("Error : Could not open clean component file.\n");
//				exit(EXIT_FAILURE);
//			}
//			fprintf(cleancomponentfile,"%f %f %lf\n",mappeakRA[i],mappeakDEC[i],mapmaxpower[i]);
//			fprintf(annfile,"CROSS W %7.4f %7.4f .1 .1 45\n",mappeakRA[i],mappeakDEC[i]);
			float xmin, xmax, ymin, ymax;
			int n,v,w;
			float *cI = (float*)malloc((n+1)*(n+1) * sizeof(float));
			float *cQ = (float*)malloc((n+1)*(n+1) * sizeof(float));
			float *cU = (float*)malloc((n+1)*(n+1) * sizeof(float));
			float *cV = (float*)malloc((n+1)*(n+1) * sizeof(float));
			for(j = 0;j < NUM_BEAMS;j++)
			{
				if(j == 6)
					continue; //not cleaning beam 6 right now.

				printf("DIAGNOSTIC: Cleaning for beam %d.\n",j);
				//sprintf(beammodelfilename,"../map144/beam%d_model/beam%d_model%03i.dat",j,j,i);
				//sprintf(beammodelfilename,"beam_models/beam%d_fit%04i.dat",j,i);
				sprintf(beammodelfilename,"beam_models/beam1_fit%04i.dat",i);
				//printf("%s\n",beammodelfilename);
//				beammodelfile = fopen(beammodelfilename,"r");
//				if(beammodelfile == NULL)
//				{
//					printf("Error : Could not open beam model file.\n");
//					exit(EXIT_FAILURE);
//				}
//				beammodellinecount = jsd_line_count(beammodelfile);
//				printf("DIAGNOSTIC: Created beam model filename %s\n",beammodelfilename);
			/*for(k = 0;k < numDays;k++)
				{
//					printf("DIAGNOSTIC: Processing day %d.\n",k);
i					if(!progress_counter)
					{
						sprintf(maptimedomainfilename,"%s/beam%d/balance%03i.dat",files[k],j,i);
						maptimedomainfile = fopen(maptimedomainfilename,"r");
						if(maptimedomainfile == NULL)
						{
							printf("Error : Could not open time domain file.\n");
							exit(EXIT_FAILURE);
						}
						sprintf(copyfilename,"%s/beam%d/clean%03i.dat",files[k],j,i);
						copyfile = fopen(copyfilename,"w");
						if(copyfile == NULL)
						{
							printf("Error : Could not open copy file.\n");
							exit(EXIT_FAILURE);
						}
						char ch;
						ch = fgetc(maptimedomainfile);
						while(ch != EOF)
						{
							fputc(ch,copyfile);
							ch = fgetc(maptimedomainfile);
						}
						fclose(copyfile);
						fclose(maptimedomainfile);							
				}*/
				beammodelfile = fopen(beammodelfilename,"r");
				if(beammodelfile == NULL)
				{
					printf("Error : Could not open beam model file %s.\n",beammodelfilename);
				fflush(stdout);
					exit(EXIT_FAILURE);
				}

				fscanf(beammodelfile,"%d %f %f %f %f\n",&n,&xmin,&xmax,&ymin,&ymax);
				printf("%d %f %f %f %f\n",n,xmin,xmax,ymin,ymax);
				fflush(stdout);


				for(v=0; v<=n; v++)	
				{	
				        for(w=0; w<=n; w++)
				        {
				                fscanf(beammodelfile,"%f ",&cI[v*(n+1)+w]);
				//              printf("%2.9f ",c[v*(n+1)+w]);
				        }
				//      printf("\n");
				}
				for(v=0; v<=n; v++)	
				{	
				        for(w=0; w<=n; w++)
				        {
				                fscanf(beammodelfile,"%f ",&cQ[v*(n+1)+w]);
				//              printf("%2.9f ",c[v*(n+1)+w]);
				        }
				//      printf("\n");
				}
				for(v=0; v<=n; v++)	
				{	
				        for(w=0; w<=n; w++)
				        {
				                fscanf(beammodelfile,"%f ",&cU[v*(n+1)+w]);
				//              printf("%2.9f ",c[v*(n+1)+w]);
				        }
				//      printf("\n");
				}
				for(v=0; v<=n; v++)	
				{	
				        for(w=0; w<=n; w++)
				        {
				                fscanf(beammodelfile,"%f ",&cV[v*(n+1)+w]);
				//              printf("%2.9f ",c[v*(n+1)+w]);
				        }
				//      printf("\n");
				}
				fclose(beammodelfile);

				printf("DIAGNOSTIC: Read model for beam %d.\n",j);
				fflush(stdout);
//				beammodellinecount = jsd_line_count(beammodelfile);
//				printf("DIAGNOSTIC: Created beam model filename %s\n",beammodelfilename);
				for(k = 0;k < numDays;k++)
				{

					if(k%7 == j)
					{
					        FluxDayData * daydata = &wappdata->daydata[k];
					        if(daydata->numRecords == 0)
					                continue;
  						for (l=0; l<daydata->numRecords; l++)
					        {
							
  							for (m=0; m<cclinecount; m++)
							{
		                                                if(fabs(ccRA[m] - daydata->records[l].RA) < clean_patch && fabs(ccDEC[m] - daydata->records[l].DEC) < clean_patch \
        		                                        && finite(daydata->records[l].stokes.I))
                		                                {
                        	                                	float dx = daydata->records[l].RA - ccRA[m];
	                        	                                float dy = daydata->records[l].DEC - ccDEC[m];
                        	                                	//float dx = ccRA[m] - daydata->records[l].RA;
	                        	                                //float dy = ccDEC[m] - daydata->records[l].DEC;
									float distance = sqrt(dx*dx+dy*dy)*60;

        	                        	                        dx = CNORMALIZE(dx,xmin,xmax);
                	                        	                dy = CNORMALIZE(dy,ymin,ymax);
									float tmp;
									tmp = chebyshev_eval_surface(dx,dy,cI,n);	
									//printf("%f %f %f\n",dx,dy,tmp);
	                                                	        daydata->records[l].stokes.I = daydata->records[l].stokes.I - (tmp*ccStokesI[m]/105.007);
									tmp = chebyshev_eval_surface(dx,dy,cQ,n);	
        	                                                	daydata->records[l].stokes.Q = daydata->records[l].stokes.Q - (tmp*ccStokesI[m]/105.007);
									tmp = chebyshev_eval_surface(dx,dy,cU,n);	
	                	                                        daydata->records[l].stokes.U = daydata->records[l].stokes.U - (tmp*ccStokesI[m]/105.007);
									tmp = chebyshev_eval_surface(dx,dy,cV,n);	
        	                	                                daydata->records[l].stokes.V = daydata->records[l].stokes.V - (tmp*ccStokesI[m]/105.007);
        	                        	                        dx = CDENORMALIZE(dx,xmin,xmax);
                	                        	                dy = CDENORMALIZE(dy,ymin,ymax);
                        	                                	//daydata->records[l].stokes.I+=exp(-0.254599*distance*distance)*(ccStokesI-mapminpower[i]);
									tmp = exp(-0.254599*distance*distance);
                        	                                	daydata->records[l].stokes.I = daydata->records[l].stokes.I + tmp*ccStokesI[m]/105.007;
								}

							}




						}

					}
					//printf("DIAGNOSTIC: Processing day %d.\n",k);
					/*if(!progress_counter)
					{
						sprintf(maptimedomainfilename,"%s/beam%d/balance%03i.dat",files[k],j,i);
						maptimedomainfile = fopen(maptimedomainfilename,"r");
						if(maptimedomainfile == NULL)
						{
							printf("Error : Could not open time domain file.\n");
							exit(EXIT_FAILURE);
						}
						sprintf(copyfilename,"%s/beam%d/clean%03i.dat",files[k],j,i);
						copyfile = fopen(copyfilename,"w");
						if(copyfile == NULL)
						{
							printf("Error : Could not open copy file.\n");
							exit(EXIT_FAILURE);
						}
						char ch;
						ch = fgetc(maptimedomainfile);
						while(ch != EOF)
						{
							fputc(ch,copyfile);
							ch = fgetc(maptimedomainfile);
						}
						fclose(copyfile);
						fclose(maptimedomainfile);							
					}*/
					//sprintf(maptimedomainfilename,"%s/beam%d/clean%03i.dat",files[k],j,i);
					/*sprintf(maptimedomainfilename,"%s/beam%d/balance%04i.dat",files[k],j,i);
					maptimedomainfile = fopen(maptimedomainfilename,"r");
					if(maptimedomainfile == NULL)
					{
						printf("Error: Could not open clean file.\n");
						exit(EXIT_FAILURE);
					}
					linecount = jsd_line_count(maptimedomainfile);
				        fgets(tempheader, 80, maptimedomainfile);
					FluxRecord *rec;
					rec = (FluxRecord *) malloc (sizeof (FluxRecord) * (linecount-1));
					if(!rec)
					{
						printf("Error: Memory allocation failed.\n");
						exit(EXIT_FAILURE);
					}
					for(l = 0;l < linecount-1; l++)
					{
					//	FluxRecord rec;
						fscanf(maptimedomainfile,"%f %f %f %lf %lf %lf %lf",&rec[l].RA,&rec[l].DEC,&rec[l].AST, \
						&rec[l].stokes.I,&rec[l].stokes.Q,&rec[l].stokes.U, &rec[l].stokes.V);
						if(fabs(mappeakRA[i] - rec[l].RA) < clean_patch && fabs(mappeakDEC[i] - rec[l].DEC) < clean_patch \
						&& finite(rec[l].stokes.I))
						{
						        //fgets(tempheader, 80, beammodelfile);
							//fgets(tempheader, 80, beammodelfile);
							double Flux_I = 0.0,Flux_Q = 0.0,Flux_U = 0.0,Flux_V = 0.0,minRA,minDEC,dx,dy,mindistance = 10.0;
							dx = rec[l].RA - mappeakRA[i];
							dy = rec[l].DEC - mappeakDEC[i];

							dx = NORMALIZE(dx,xmin,xmax);
							dy = NORMALIZE(dy,ymin,ymax);

							rec[l].stokes.I -= chebyshev_eval_surface(dx,dy,cI,n)
							rec[l].stokes.Q -= chebyshev_eval_surface(dx,dy,cQ,n)
							rec[l].stokes.U -= chebyshev_eval_surface(dx,dy,cU,n)
							rec[l].stokes.V -= chebyshev_eval_surface(dx,dy,cV,n)

						*/	/*for(m = 0;m < beammodellinecount - 2;m++)
							{
								float distance1,distance2,ddistance;
								FluxRecord beammodelrec;
								fscanf(beammodelfile,"%f %f %f %lf %lf %lf %lf",&beammodelrec.RA,&beammodelrec.DEC, \
								&beammodelrec.AST,&beammodelrec.stokes.I,&beammodelrec.stokes.Q, &beammodelrec.stokes.U, \
								&beammodelrec.stokes.V);	
								distance1 = dx - beammodelrec.RA;
								distance2 = dy - beammodelrec.DEC;
								ddistance = sqrt(distance1*distance1 + distance2*distance2);
								if(ddistance < mindistance)
								{
									Flux_I = beammodelrec.stokes.I;
									Flux_Q = beammodelrec.stokes.Q;
									Flux_U = beammodelrec.stokes.U;
									Flux_V = beammodelrec.stokes.V;
									mindistance = ddistance;
								}
							}*/
						//}
					//}
//					printf("DIAGNOSTIC: done day %d.\n",k);
				}
//				printf("DIAGNOSTIC: done cleaning for beam %d.\n",j);
			}
			//free(cI);
			//free(cQ);
			//free(cU);
			//free(cV);
			//fclose(cleancomponentfile);
//			printf("DIAGNOSTIC: done cleaning for channel %d.\n",i);
//sguram	}

//		printf("DIAGNOSTIC: RMS sum of map %lf\n",rms);
//		rms = sqrt(rms/rms_counter);
//		printf("DIAGNOSTIC: RMS of map %lf\n",rms);
//		printf("DIAGNOSTIC: RMS counter %ld\n",rms_counter);
//		if(!progress_counter)
//			sprintf(dirtymapname,"beam8_Icube.fits");
//		progress_counter++;
//		stop_rms = 10000.0;
//HACK STARTS
//	FluxWappData *wappdata;
//	char wapp[10];
//	sprintf(wapp,"beam8");
//	printf("Wapp %s.\n",wapp);
	// allocate and initialize the wapp data days
//	numDays1 = numDays*7;
//	wappdata = fluxwappdata_alloc(wapp, files, numDays1);
			//md.n3 = 1;
			//md.n3 = md.highchan - md.lowchan;
//	printf("After start fits.\n");
	//		float *dataI, *dataQ, *dataU, *dataV, *weight;
	//		int numbytes; 
	//		numbytes = md.n1 * md.n2 * sizeof (float);
	//		dataI  = (float *) malloc (numbytes);
	//		dataQ  = (float *) malloc (numbytes);
	//		dataU  = (float *) malloc (numbytes);
	//		dataV  = (float *) malloc (numbytes);
	//		weight  = (float *) malloc (numbytes);
	
	//		if (!dataI || !dataQ || !dataU || !dataV || !weight) {
	//			printf("ERROR: memory allocation of %i bytes failed!\n", numbytes);
	//			return EXIT_FAILURE;
	//		}
//	printf("After psf lookup.\n");
//	fluxwappdata_readchan(wappdata, md.lowchan);
//	printf("After fluxwappdata_readchan.\n");
			grid_data(wappdata, &md, dataI, dataQ, dataU, dataV, weight);
	//	printf("After grid data.\n");
			write_fits_planes(dataI, dataQ, dataU, dataV, weight);
/*	if(showprogress)
	{
		writefits_plane(progressfileI, dataI, &progressfile_hpar);
		writefits_plane(progressfileQ, dataQ, &progressfile_hpar);
		writefits_plane(progressfileU, dataU, &progressfile_hpar);
		writefits_plane(progressfileV, dataV, &progressfile_hpar);
	}*/
//	printf("After write fits.\n");
	//finish_fits_cubes();
//	printf("After finish fits.\n");
			//free(dataI);
			//free(dataQ);
			//free(dataU);
			//free(dataV);
			//free(weight);
			printf("DIAGNOSTIC: Iteration %d complete.\n",i);
		}
		//printf("DIAGNOSTIC: Here.\n");
		//fflush(stdout);
//HACK ENDS
	//}while(progress_counter < 20);
	}while(0);
		//printf("DIAGNOSTIC: Here.\n");
		//fflush(stdout);
		finish_fits_cubes();
		printf("DIAGNOSTIC: Here.\n");
//sss		fflush(stdout);
	
//	fclose(annfile);

//Restoration starts
//	printf("DIAGNOSTIC: Starting restoration.\n");

//temporary stuff
//	FluxWappData *wappdata;
//	char wapp[10];
//	float *dataI, *dataQ, *dataU, *dataV, *weight;
//temporary stuff

/*	for(i = md.lowchan;i < md.highchan;i++)
	{
		sprintf(cleancomponentfilename,"cleancomponent%03i.dat",i);
		cleancomponentfile = fopen(cleancomponentfilename,"r");
		if(cleancomponentfile == NULL)
		{
			printf("Error : Could not open clean component file.\n");
			exit(EXIT_FAILURE);
		}
		linecount = jsd_line_count(cleancomponentfile);
	        fgets(tempheader, 80, cleancomponentfile);
		for(l = 0;l < linecount-1; l++)
		{
			double RA,DEC,power;
			fscanf(cleancomponentfile,"%lf %lf %lf\n",&RA,&DEC,&power);
			printf("Reading clean component file RA %lf, DEC %lf, Power %lf\n",RA,DEC,power);
			for(j = 0;j < NUM_BEAMS;j++)
			{
				for(k = 0;k < numDays;k++)
				{
					int linecount1;
					sprintf(maptimedomainfilename,"%s/beam%d/clean%03i.dat",files[k],j,i);
					maptimedomainfile = fopen(maptimedomainfilename,"r");
					if(maptimedomainfile == NULL)
					{
						printf("Error: Could not open clean file.\n");
						exit(EXIT_FAILURE);
					}
					linecount1 = jsd_line_count(maptimedomainfile);
				        fgets(tempheader, 80, maptimedomainfile);
					FluxRecord *rec;
					rec = (FluxRecord *) malloc (sizeof (FluxRecord) * (linecount1-1));
					for(m = 0;m < linecount1 - 1;m++)
					{
						double dx,dy,distance;
						fscanf(maptimedomainfile,"%f %f %f %lf %lf %lf %lf",&rec[m].RA,&rec[m].DEC,&rec[m].AST, \
						&rec[m].stokes.I,&rec[m].stokes.Q,&rec[m].stokes.U, &rec[m].stokes.V);
						if(fabs(RA - rec[m].RA) < clean_patch && fabs(DEC - rec[m].DEC) < clean_patch && finite(rec[m].stokes.I))
						{
							dx = (RA - rec[m].RA)*60;
							dy = (DEC - rec[m].DEC)*60;
							distance = sqrt(dx*dx+dy*dy);	
							rec[m].stokes.I=rec[m].stokes.I+exp(-0.254599*distance*distance)*loop_gain*(power-mapminpower[i]);
							//HACK
						}
					}
					fclose(maptimedomainfile);
					maptimedomainfile = fopen(maptimedomainfilename,"w");
					fprintf(maptimedomainfile, "RA DEC AST I Q U V\n");

					for(m = 0;m < linecount1 - 1; m++)
					{
						fprintf(maptimedomainfile, "%2.8f %2.8f %8.2f %4.6f %4.6f %4.6f %4.6f\n",rec[m].RA, rec[m].DEC, \
						rec[m].AST, rec[m].stokes.I, rec[m].stokes.Q, rec[m].stokes.U, rec[m].stokes.V);
					}
					fclose(maptimedomainfile);
				}
			}
		}
		fclose(cleancomponentfile);		
	}	
*/
//	sprintf(wapp,"beam8");
//	printf("Wapp %s.\n",wapp);
	// allocate and initialize the wapp data days
//	numDays1 = numDays*7;
//	wappdata = fluxwappdata_alloc(wapp, files, numDays1);
//	md.n3 = 1;
/*	start_fits_cubes(wapp,&md);
//	printf("After start fits.\n");
//	float *dataI, *dataQ, *dataU, *dataV, *weight;
//temporary
	int numbytes; 
	numbytes = md.n1 * md.n2 * sizeof (float);
	dataI  = (float *) malloc (numbytes);
	dataQ  = (float *) malloc (numbytes);
	dataU  = (float *) malloc (numbytes);
	dataV  = (float *) malloc (numbytes);
	weight  = (float *) malloc (numbytes);
	
	if (!dataI || !dataQ || !dataU || !dataV || !weight) {
		printf("ERROR: memory allocation of %i bytes failed!\n", numbytes);
		return EXIT_FAILURE;
	}
//temporary
//	printf("After psf lookup.\n");
//	fluxwappdata_readchan(wappdata, md.lowchan);
//	printf("After fluxwappdata_readchan.\n");
	grid_data(wappdata, &md, dataI, dataQ, dataU, dataV, weight);
//	printf("After grid data.\n");
	write_fits_planes(dataI, dataQ, dataU, dataV, weight);
//	if(showprogress)
//		writefits_plane(progressfileI, dataI, &progressfile_hpar);
//	printf("After write fits.\n");
	finish_fits_cubes();
//	printf("After finish fits.\n");
	free(dataI);
	free(dataQ);
	free(dataU);
	free(dataV);
	printf("DIAGNOSTIC: Restoration complete.\n");*/
/*	if(showprogress)
	{
		//close the fits cube
		writefits_pad_end(progressfileI, &progressfile_hpar);
		fclose(progressfileI);
		fclose(progressfileQ);
		fclose(progressfileU);
		fclose(progressfileV);

		//change the fits header on disk to update the number of planes of the cube
		progressfile_hpar.naxis[2] = progress_counter+1;
		progressfileI = fopen("CLEAN_ProgressI.fits", "r+");
		writefits_header(progressfileI, &progressfile_hpar);
		progressfile_hpar.naxis[2] = progress_counter;
		progressfileQ = fopen("CLEAN_ProgressQ.fits", "r+");
		writefits_header(progressfileQ, &progressfile_hpar);
		progressfileU = fopen("CLEAN_ProgressU.fits", "r+");
		writefits_header(progressfileU, &progressfile_hpar);
		progressfileV = fopen("CLEAN_ProgressV.fits", "r+");
		writefits_header(progressfileV, &progressfile_hpar);
		fclose(progressfileI);
		fclose(progressfileQ);
		fclose(progressfileU);
		fclose(progressfileV);
	}*/
//Restoration ends
	//free(files);
	//fluxwappdata_free(wappdata);
	//free(wappdata);
	printf("DIAGNOSTIC: Multibeam CLEAN complete!\n");
	return EXIT_SUCCESS;
}


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
#include "leakagefit.h"
//#include "clean.h"
int multibeam;
static void print_usage(const char * prog)
{
	printf(
	"Usage: %s <fcen> <fitsmap> <lowchan> <highchan> "
	"<grid patch radius> <clean patch radius> <clean loopgain> <stop rms> <title> <show progress>\n"
	"\n"
	"\t fcen - center freqyency for the band"
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
	int numDays,i,j,k,l,m,beammodellinecount;
	char * dirtymapname;
	char beammodelfilename[64],cleancomponentfilename[64];
	float mappeakRA[MAX_CHANNELS],mappeakDEC[MAX_CHANNELS],mapmaxpower[MAX_CHANNELS],mapminpower[MAX_CHANNELS];
	float clean_patch;
	char *wapp;
	FILE * dirtymap_file, *beammodelfile, *cleancomponentfile;
	header_param_list dirtymap_hpar;

	header_param_list model_hpar;
	int data_len,indx_max;
	float *planeI,*planeQ,*planeU,*planeV,maxI;
	int nx,ny;
	double crpx,crpy,delx,dely,crvx,crvy;
	
	/* Process command line arguments */ 
	/* Convert args into more usable units were required */

	if (argc != 9) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	} else {
		wapp = argv[1];
		md.fcen = (float)atof(argv[2]);
		//printf("%lf\n",md.fcen);
		//exit(1);
		dirtymapname = argv[3];
		md.lowchan = atoi(argv[4]);
		md.highchan = atoi(argv[5]);
		md.patch = atoi(argv[6]);
		clean_patch = (float )atof(argv[7]);
		md.title = argv[8];
		dirtymap_file = fopen(dirtymapname, "r");
		readfits_header(dirtymap_file, &dirtymap_hpar);
		md.beamwidths = md.patch;
		md.cellsize = fabs(dirtymap_hpar.cdelt[0]); 
		md.RAcen = dirtymap_hpar.crval[0];
		md.DECcen = dirtymap_hpar.crval[1];
		md.RArange =  fabs(dirtymap_hpar.naxis[0]*dirtymap_hpar.cdelt[0]);
		md.DECrange = fabs(dirtymap_hpar.naxis[1]*dirtymap_hpar.cdelt[1]);
		md.n1 = dirtymap_hpar.naxis[0];
		md.n2 = dirtymap_hpar.naxis[1];
		//md.n3 = dirtymap_hpar.naxis[2];
		md.n3 = (int)((md.highchan-md.lowchan)/10);
		md.fwhm = 2.0;
		//md.fstart = dirtymap_hpar.crval[2];
		md.fstart = md.fcen - (md.lowchan-(MAX_CHANNELS/2-1)) * (0.042);
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
		printf("Start Frequency: %g\n", md.fstart);

	}

	//numDays = get_date_dirs("./", &files);
	numDays = get_date_dirs("/n/ras/processed/FIELD1/run5/band0", &files);
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

        wappdata = fluxwappdata_alloc(wapp, files, numDays);

        init_psf_map(md.fwhm, md.cellsize, md.beamwidths);
	init_psf_lookup_table(65537, 9.1);
	
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
		//for(i = md.lowchan;i < md.highchan;i++)
		for(i = md.lowchan;i < md.highchan;i=i+10)
		{
// i = 0;
			printf("Channel %d\n",i);
                        fluxwappdata_readchan(wappdata, i, CLEAN);

			float ratio = 0.0;
			float bmwidth = 225.0+i*0.0077;
			//float bmwidth = 225.0+2000*0.0077;
			float bmwidth_true = sqrt(bmwidth*bmwidth-120.0*120.0);
			ratio = bmwidth/bmwidth_true;
			bmwidth_true/=60.0;
			//printf("%f %f %f\n",bmwidth,bmwidth_true,ratio);

			//vleakagefit(wappdata,i);

			float xmin, xmax, ymin, ymax;
			int maxx,maxy;
			int n=12,v,w;
			//float *cI = (float*)malloc((n+1)*(n+1) * sizeof(float));
			//float *cQ = (float*)malloc((n+1)*(n+1) * sizeof(float));
			//float *cU = (float*)malloc((n+1)*(n+1) * sizeof(float));
			//float *cV = (float*)malloc((n+1)*(n+1) * sizeof(float));
			//for(j = 0;j < NUM_BEAMS;j++)
			for(j = 0;j < 7;j++)
			{
				if(j == 6)
					continue; //not cleaning beam 6 right now.

				//printf("DIAGNOSTIC: Cleaning for beam %d.\n",j);
				//sprintf(beammodelfilename,"beam_models/beam%d_fit%04i.dat",j,i);
				//sprintf(beammodelfilename,"beam_models/beam5_fit%04i.dat",i);
				sprintf(beammodelfilename,"beam_models/MODEL_BEAM%d_Icube.fits",j);
				//printf("DIAGNOSTIC: opened %s.\n",beammodelfilename);

				beammodelfile = fopen(beammodelfilename,"r");
				if(beammodelfile == NULL)
				{
					printf("Error : Could not open beam model file %s.\n",beammodelfilename);
					fflush(stdout);
					exit(EXIT_FAILURE);
				}

				readfits_header(beammodelfile, &model_hpar);
			        //allocate memory for plane
			        data_len = model_hpar.naxis[0] * model_hpar.naxis[1];
			        nx = model_hpar.naxis[0];
				ny = model_hpar.naxis[1];
			        crpx = model_hpar.crpix[0];
				crpy = model_hpar.crpix[1];
			        delx = model_hpar.cdelt[0];
				dely = model_hpar.cdelt[1];
			        crvx = model_hpar.crval[0];
				crvy = model_hpar.crval[1];

			        planeI = (float*) calloc(data_len, sizeof (float));
			        planeQ = (float*) calloc(data_len, sizeof (float));
			        planeU = (float*) calloc(data_len, sizeof (float));
			        planeV = (float*) calloc(data_len, sizeof (float));
		                readfits_plane(beammodelfile, planeI, &model_hpar);
				maxI = 0;
				maxx = 0;
				maxy = 0;
                		for (k=0; k<data_len; k++)
		                {
                		        if (!IS_BLANK_PIXEL(planeI[k]) && isfinite(planeI[k]))
		                        {
                                		if(maxI < planeI[k])
                		                {
		                                        maxI = planeI[k];
                                		        indx_max = k;
                		                }
		//                                if(min_response[startchannel+i] > plane_data[j])
	//                                        min_response[startchannel+i] = plane_data[j];
                		        }
		                }

				maxx = (int)(indx_max%nx);
				maxy = (int)(indx_max - maxx)/ny;
				//printf("maxx %d maxy %d indxmax %d maxI %f\n",maxx,maxy,indx_max,maxI);
                		for (k=0; k<data_len; k++)
		                {
                		        if (!IS_BLANK_PIXEL(planeI[k]) && isfinite(planeI[k]))
		                        {
						planeI[k]/=maxI;
					}
				}
				fclose(beammodelfile);

				sprintf(beammodelfilename,"beam_models/MODEL_BEAM%d_Qcube.fits",j);
				//printf("DIAGNOSTIC: opened %s.\n",beammodelfilename);

				beammodelfile = fopen(beammodelfilename,"r");
				if(beammodelfile == NULL)
				{
					printf("Error : Could not open beam model file %s.\n",beammodelfilename);
					fflush(stdout);
					exit(EXIT_FAILURE);
				}
				readfits_header(beammodelfile, &model_hpar);
		                readfits_plane(beammodelfile, planeQ, &model_hpar);
                		for (k=0; k<data_len; k++)
		                {
                		        if (!IS_BLANK_PIXEL(planeQ[k]) && isfinite(planeQ[k]))
		                        {
						planeQ[k]/=maxI;
					}
				}
				fclose(beammodelfile);

				sprintf(beammodelfilename,"beam_models/MODEL_BEAM%d_Ucube.fits",j);
				//printf("DIAGNOSTIC: opened %s.\n",beammodelfilename);

				beammodelfile = fopen(beammodelfilename,"r");
				if(beammodelfile == NULL)
				{
					printf("Error : Could not open beam model file %s.\n",beammodelfilename);
					fflush(stdout);
					exit(EXIT_FAILURE);
				}
				readfits_header(beammodelfile, &model_hpar);
		                readfits_plane(beammodelfile, planeU, &model_hpar);
                		for (k=0; k<data_len; k++)
		                {
                		        if (!IS_BLANK_PIXEL(planeU[k]) && isfinite(planeU[k]))
		                        {
						planeU[k]/=maxI;
					}
				}
				fclose(beammodelfile);

				sprintf(beammodelfilename,"beam_models/MODEL_BEAM%d_Vcube.fits",j);
				//printf("DIAGNOSTIC: opened %s.\n",beammodelfilename);

				beammodelfile = fopen(beammodelfilename,"r");
				if(beammodelfile == NULL)
				{
					printf("Error : Could not open beam model file %s.\n",beammodelfilename);
					fflush(stdout);
					exit(EXIT_FAILURE);
				}
				readfits_header(beammodelfile, &model_hpar);
		                readfits_plane(beammodelfile, planeV, &model_hpar);
                		for (k=0; k<data_len; k++)
		                {
                		        if (!IS_BLANK_PIXEL(planeV[k]) && isfinite(planeV[k]))
		                        {
						planeV[k]/=maxI;
					}
				}
				fclose(beammodelfile);
				/*fscanf(beemmodelfile,"%d %f %f %f %f\n",&n,&xmin,&xmax,&ymin,&ymax);
				printf("%d %f %f %f %f\n",n,xmin,xmax,ymin,ymax);
				fflush(stdout);


				for(v=0; v<=n; v++)	
				{	
				        for(w=0; w<=n; w++)
				        {
				                fscanf(beammodelfile,"%f ",&cI[v*(n+1)+w]);
				                //printf("%2.9f ",cI[v*(n+1)+w]);
				        }
				        //printf("\n");
				}
				for(v=0; v<=n; v++)	
				{	
				        for(w=0; w<=n; w++)
				        {
				                fscanf(beammodelfile,"%f ",&cQ[v*(n+1)+w]);
				                //printf("%2.9f ",cQ[v*(n+1)+w]);
				        }
				        //printf("\n");
				}
				for(v=0; v<=n; v++)	
				{	
				        for(w=0; w<=n; w++)
				        {
				                fscanf(beammodelfile,"%f ",&cU[v*(n+1)+w]);
				                //printf("%2.9f ",cU[v*(n+1)+w]);
				        }
				        //printf("\n");
				}
				for(v=0; v<=n; v++)	
				{	
				        for(w=0; w<=n; w++)
				        {
				                fscanf(beammodelfile,"%f ",&cV[v*(n+1)+w]);
				                //printf("%2.9f ",cV[v*(n+1)+w]);
				        }
				        //printf("\n");
				}
				fflush(stdout);
				fclose(beammodelfile);
				*/
				//printf("DIAGNOSTIC: Read model for beam %d.\n",j);
				fflush(stdout);
				
				//only for leakage maps
				for(k = 0;k < numDays;k++)
				{
					if(k%7 == j)
					{
						FluxDayData * daydata = &wappdata->daydata[k];
						if(daydata->numRecords == 0)
						                continue;
						for (l=0; l<daydata->numRecords; l++)
					        {
	                                		daydata->records[l].stokes.I = 0.0;
	                                		daydata->records[l].stokes.Q = 0.0;
		                                	daydata->records[l].stokes.U = 0.0;
		                                	daydata->records[l].stokes.V = 0.0;
						}
					}
				}


				for(k = 0;k < numDays;k++)
				{


					if(k%7 == j)
					//if(1)
					{
					        FluxDayData * daydata = &wappdata->daydata[k];
					        if(daydata->numRecords == 0)
					                continue;
		
						//FILE * compfile;
						//char compfilename[165];
						//sprintf(compfilename,"%s/beam%d/comp%04i.dat",daydata->mjd,j,i);
  				
						//compfile = fopen(compfilename,"w");
						//if(compfile == NULL)
						//	printf("cant open compfile %s\n",compfilename);
	
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
									int xcd,ycd,arrcd;
									//xcd = (int)(maxx+dx/delx);
									//ycd = (int)(maxy+dy/dely);
									//xcd = (int)(maxx+dx*1.1/delx);
									//ycd = (int)(maxy+dy*1.1/dely);
									xcd = (int)(maxx+dx*ratio/delx);
									ycd = (int)(maxy+dy*ratio/dely);
									arrcd = ycd*nx+xcd;
									//printf("x %d y %d arr %d nx %d ny %d\n",xcd,ycd,arrcd,nx,ny);
									//if(l == 100)
									//	exit(1);
        	                        	                        //dx = CNORMALIZE(dx,xmin,xmax);
                	                        	                //dy = CNORMALIZE(dy,ymin,ymax);
									float tmp;
									tmp = planeI[arrcd];	
									//printf("%f %f %f\n",dx,dy,tmp);
	                                                	        //daydata->records[l].stokes.I = daydata->records[l].stokes.I - (tmp*ccStokesI[m]/105.007);
	                                                	        //daydata->records[l].stokes.I = daydata->records[l].stokes.I - (tmp*ccStokesI[m]/1000.0);
	                                                	        daydata->records[l].stokes.I += tmp*ccStokesI[m]/1000.0;
									//fprintf(compfile,"%2.9f %2.9f %2.9f %2.9f %2.9f\n",daydata->records[l].RA,daydata->records[l].DEC,daydata->records[l].stokes.Q,(tmp*ccStokesI[m]/105.007),(daydata->records[l].stokes.Q - (tmp*ccStokesI[m]/105.007)));	
	
									tmp = planeQ[arrcd];	
        	                                                	//daydata->records[l].stokes.Q = daydata->records[l].stokes.Q - (tmp*ccStokesI[m]/105.007);
        	                                                	//daydata->records[l].stokes.Q = daydata->records[l].stokes.Q - (tmp*ccStokesI[m]*1.0/1000.0);
        	                                                	daydata->records[l].stokes.Q += tmp*ccStokesI[m]*1.0/1000.0;
									tmp = planeU[arrcd];	
	                	                                        //daydata->records[l].stokes.U = daydata->records[l].stokes.U - (tmp*ccStokesI[m]/105.007);
	                	                                        //daydata->records[l].stokes.U = daydata->records[l].stokes.U - (tmp*ccStokesI[m]*1.0/1000.0);
	                	                                        daydata->records[l].stokes.U += tmp*ccStokesI[m]*1.0/1000.0;
									tmp = planeV[arrcd];	
        	                	                                //daydata->records[l].stokes.V = daydata->records[l].stokes.V - (tmp*ccStokesI[m]/105.007);
        	                	                                //daydata->records[l].stokes.V = daydata->records[l].stokes.V - (tmp*ccStokesI[m]*1.0/1000.0);
        	                	                                daydata->records[l].stokes.V += tmp*ccStokesI[m]*1.0/1000.0;
									//fprintf(compfile,"%2.9f %2.9f %2.9f %2.9f %2.9f %2.9f\n",daydata->records[l].RA,daydata->records[l].DEC,daydata->records[l].stokes.V,(tmp*ccStokesI[m]/105.007),(daydata->records[l].stokes.V - (tmp*ccStokesI[m]/105.007)),daydata->records[l].stokes.I);	
        	                        	                        //dx = CDENORMALIZE(dx,xmin,xmax);
                	                        	                //dy = CDENORMALIZE(dy,ymin,ymax);
									float distance = sqrt(dx*dx+dy*dy)*60;
									//float distance2 = dx*dx+dy*dy*60*60;
                        	                                	//daydata->records[l].stokes.I+=exp(-0.254599*distance*distance)*(ccStokesI-mapminpower[i]);
									//tmp = exp(-0.254599*distance*distance);
									//tmp = exp(-0.254599*distance2);
									//tmp = exp(-0.22633*distance2);
									//tmp = exp(-0.22633*distance*distance);
									tmp = exp(-2.772589*(distance/bmwidth_true)*(distance/bmwidth_true));
                        	                                	//daydata->records[l].stokes.I = daydata->records[l].stokes.I + tmp*ccStokesI[m]/105.007;
                        	                                	//daydata->records[l].stokes.I = daydata->records[l].stokes.I + tmp*ccStokesI[m]/1000.0;
                        	                                	daydata->records[l].stokes.I = daydata->records[l].stokes.I - tmp*ccStokesI[m]/1000.0;
								}
								/*else
								{
									fprintf(compfile,"%2.9f %2.9f %2.9f %2.9f %2.9f\n",daydata->records[l].RA,daydata->records[l].DEC,daydata->records[l].stokes.Q,0.0,daydata->records[l].stokes.Q);	
								}*/

							}

						}
						//fclose(compfile);
					}
				}
			}
			//vleakagefit(wappdata,i);
			grid_data(wappdata, &md, dataI, dataQ, dataU, dataV, weight);

//writeout V/I
			/*int ii;
			FILE *VvsI = fopen("V_vs_I.dat","w");
			fprintf(VvsI,"# RA DEC I V\n");
  			for (ii=0; ii<cclinecount; ii++)
			{
				int xx,yy,array;
				xx = (int)(ccRA[ii] - md.ramin)/md.cellsize;				
				yy = (int)(ccDEC[ii] - md.decmin)/md.cellsize;
				array = yy*md.n1 + xx;			
				if(isfinite(dataI[array])) 
				{	
					fprintf(VvsI,"%2.6f %2.6f %2.6f %2.6f\n",ccRA[ii],ccDEC[ii],dataI[array],dataV[array]);
				}
			}
			fclose(VvsI);
			*/
//
			write_fits_planes(dataI, dataQ, dataU, dataV, weight);
			printf("DIAGNOSTIC: Iteration %d complete.\n",i);
		}
	}while(0);
	finish_fits_cubes();
	//printf("DIAGNOSTIC: Here.\n");
	
	printf("DIAGNOSTIC: Multibeam CLEAN complete!\n");
	return EXIT_SUCCESS;
}

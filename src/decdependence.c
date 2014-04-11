#include "decdependence.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "jsd/jsd_fit.h"
#include "stats.h"

/*
 * use a single channel worth of data
 * iterate over each beam
 * iterate over every day
 * put points into DEC bins
 * reject outliers in the bins
 * fit a final curve to the dec bins
 * then subtract the amount of that curve from the data
 *
 */


static void compute_dirty_bins (double **binarrayX, int countX[], double binX[], int num_bins)
{
	int b;
	for (b=0; b<num_bins; b++)
	{
		binX[b] = compute_mean(binarrayX[b], 0,countX[b]);
	}
}

static void compute_clean_bins (double **binarrayX, int *countX, double *binX, int num_bins, float nsigma, int *outliercounts)
{
	int b, i;
	double mean, sigma;
	int outlier;
	double *binarray;
	int count;
	mean = sigma = 0.0;
	memset(outliercounts, 0, num_bins*sizeof(int));
	for (b=0; b<num_bins; b++)
	{
		//repeat until no more outliers
		count = countX[b];
		do {
			outlier = 0;
			binarray = binarrayX[b];
			//if((outliercounts[b] < (int)(0.5*count)) && count != 0)
			if(( count >= 10))
			{
				mean = compute_mean(binarray, 0,count);
				sigma = compute_sigma(binarray, count, mean);

				//do the cleaning by setting outliers to NAN
				for (i=0; i<count; i++) {
					//if (fabs(mean-binarray[i]) > nsigma*sigma) {
					if (binarray[i]-mean > -0.01*sigma) {
						binarray[i] = NAN;
						outlier = 1;
						outliercounts[b]++;
					}
				}
			}
//			else if (count != 0)
//			else 
//				printf("Bad stats bin: %d, count:%d outliers:%d mean:%f\n",b,count,outliercounts[b],mean);
		} while (outlier);
		binX[b] = mean;
	}
}


static void create_dec_bins(FluxWappData * wappdata, double **binarrayI,  double **binarrayQ,  double **binarrayU,  double **binarrayV, int *counts, int beam, float decmin, float decgrain, int num_bins)
{
	int d, r;

	memset(counts, 0, num_bins*sizeof(int));

	for (d=0; d<wappdata->numDays; d++)
	{
		FluxDayData * daydata = &wappdata->daydata[d];

//		if (d%7!=beam)
	//		continue; //TODO: this is a workaround for non-beam aware datastructures

		for (r=0; r<daydata->numRecords; r++)
		{
			double I = daydata->records[r].stokes.I;
			double Q = daydata->records[r].stokes.Q;
			double U = daydata->records[r].stokes.U;
			double V = daydata->records[r].stokes.V;
			double DEC = daydata->records[r].DEC;
			int bin = (int) floor((DEC - decmin) / decgrain);

			if (bin<0 || bin>=num_bins) continue;
			if (isfinite(I)) {
				double *binarray;
				binarray = binarrayI[bin];
				binarray[counts[bin]] = I;
				binarray = binarrayQ[bin];
				binarray[counts[bin]] = Q;
				binarray = binarrayU[bin];
				binarray[counts[bin]] = U;
				binarray = binarrayV[bin];
				binarray[counts[bin]] = V;
				counts[bin]++;
			}
		}
	}
}


#define BIN_ORDER 15
//curve fit the binned data
static void fit_dec_bins(FluxWappData * wappdata, double *binI, double *binQ, double *binU, double *binV, int num_bins, int day, float decmin, float decmax, float decgrain,int chan)
{

	int i, d, r;
	float nsigma = 4.0; //for the final curve fit exclusion
	double chisq;
	double cI[BIN_ORDER+1];
	double cQ[BIN_ORDER+1];
	double cU[BIN_ORDER+1];
	double cV[BIN_ORDER+1];
	double min, max;
	double *binDEC;

	//ssg beam gain corrections
//	float beamgain[7][MAX_CHANNELS/2];
//	float Qleakage[7][MAX_CHANNELS/2];
//	float Uleakage[7][MAX_CHANNELS/2];
//	FILE * beamgainfile;
//	beamgainfile = fopen("beam_all_averageband.dat","r");
//	if (beamgainfile == NULL)
//		printf("ERORR:Unable to open beamgain file.\n");

//	int temp;
//	float tmpflt;
//	int b;
//	printf("**********beam gains*****************\n");
//	for(i=0;i<MAX_CHANNELS/2;i++)
//	{
//		if(i<150 || i>1894)
//		{
//			beamgain[0][i] = 1.0;
//			beamgain[1][i] = 1.0;
//                      beamgain[2][i] = 1.0;
//                        beamgain[3][i] = 1.0;
//                        beamgain[4][i] = 1.0;
 //                       beamgain[5][i] = 1.0;
  //              	beamgain[6][i] = 1.0;
 
//		}
//		else
//		{
//		      fscanf(beamgainfile,"%d %f %f %f %f %f %f",&temp,&beamgain[0][i],&beamgain[1][i],&beamgain[2][i],&beamgain[3][i],&beamgain[4][i],&beamgain[5][i]);
//		      beamgain[6][i] = beamgain[1][i];
//		}
//	        printf("%d %f %f %f %f %f %f %f\n",i,beamgain[0][i],beamgain[1][i],beamgain[2][i],beamgain[3][i],beamgain[4][i],beamgain[5][i],beamgain[6][i]);
//	}
//	fclose(beamgainfile);
	//ssg beam gain coorections

//	FILE * Qleakagefile;
//	Qleakagefile = fopen("QI_spline.dat","r");
//	if (Qleakagefile == NULL)
//		printf("ERORR:Unable to open Qleakage file.\n");

	
//	for(b = 0;b<6;b++)
//	{
//	for(i=0;i<MAX_CHANNELS/2;i++)
//	{
//		if(b==6)
//		{
//			Qleakage[b][i] = 0.0;
//		}
//		if(i<150 || i>1894 || (i > 1369 && i < 1380))
//		{
//			Qleakage[b][i] = 0.0;
 
//		}
//		else
//		{
//		      fscanf(Qleakagefile,"%f %f %f",&tmpflt,&tmpflt,&Qleakage[b][i]);
//		}
//	}
//	}
//	fclose(Qleakagefile);

//	printf("**********Q leakage*****************\n");
//	for(i=0;i<MAX_CHANNELS/2;i++)
//	        printf("%d %f %f %f %f %f %f %f\n",i,Qleakage[0][i],Qleakage[1][i],Qleakage[2][i],Qleakage[3][i],Qleakage[4][i],Qleakage[5][i],Qleakage[6][i]);


//	FILE * Uleakagefile;
//	Uleakagefile = fopen("UI_spline.dat","r");
//	if (Uleakagefile == NULL)
//		printf("ERORR:Unable to open Uleakage file.\n");

	
//	for(b = 0;b<6;b++)
//	{
//	for(i=0;i<MAX_CHANNELS/2;i++)
//	{
//		if(b==6)
//		{
//			Uleakage[b][i] = 0.0;
//		}
//		if(i<150 || i>1894|| (i > 1369 && i < 1380))
//		{
//			Uleakage[b][i] = 0.0;
//		}
//		else
//		{
//		      fscanf(Uleakagefile,"%f %f %f",&tmpflt,&tmpflt,&Uleakage[b][i]);
//		}
//	}
//	}
//	fclose(Uleakagefile);

//	printf("**********U leakage*****************\n");
//	for(i=0;i<MAX_CHANNELS/2;i++)
//	        printf("%d %f %f %f %f %f %f %f\n",i,Uleakage[0][i],Uleakage[1][i],Uleakage[2][i],Uleakage[3][i],Uleakage[4][i],Uleakage[5][i],Uleakage[6][i]);
//	exit(1);
	binDEC = calloc(num_bins, sizeof(double));
	for (i=0; i<num_bins; i++) {
		binDEC[i] = i*decgrain+decmin; //TODO: compute this more intellegently
	}

	jsd_minmax(binDEC, num_bins, &min, &max);
	jsd_normalize(binDEC, num_bins, min, max);
	//printf("decmin:%g min:%g decmax:%g, max: %g\n", decmin, min, decmax, max);
/*	FILE * binfile;
	char binfilename[40];
	sprintf(binfilename,"%s/beam%d/decbins%04i.dat",wappdata->daydata[day].mjd,day%7,chan);
	binfile = fopen(binfilename,"w");
	fprintf(binfile,"#DEC binI binQ binU binV\n");
	for(i=0;i<num_bins;i++)
	{
		fprintf(binfile,"%2.6f %2.6f %2.6f %2.6f %2.6f\n",binDEC[i],binI[i],\
		binQ[i],binU[i],binV[i]);
	}
	fclose(binfile);
*/
	FILE * polyfile;
	char polyfilename[40];
	sprintf(polyfilename,"%s/beam%d/polyfile%04i.dat",wappdata->daydata[day].mjd,day%7,chan);
//	sprintf(polyfilename,"%s/beam%d/polyfile2730.dat",wappdata->daydata[day].mjd,day%7,chan);
	polyfile = fopen(polyfilename,"w");
//	polyfile = fopen(polyfilename,"r");
	jsd_poly_fit(binDEC, binI, num_bins, nsigma, cI, BIN_ORDER, &chisq);
	jsd_poly_fit(binDEC, binQ, num_bins, nsigma, cQ, BIN_ORDER, &chisq);
	jsd_poly_fit(binDEC, binU, num_bins, nsigma, cU, BIN_ORDER, &chisq);
	jsd_poly_fit(binDEC, binV, num_bins, nsigma, cV, 5, &chisq);
	int m;
	for(m =0;m<=BIN_ORDER;m++)
	{
		fprintf(polyfile,"+ %2.6lf*x**%d ",cI[m],m);
		//fscanf(polyfile,"%lf ",&cI[m]);
		//printf("%2.6lf ",cI[m]);
	}
	fprintf(polyfile,"\n");
	for(m =0;m<=BIN_ORDER;m++)
	{
		fprintf(polyfile,"%2.6lf ",cQ[m]);
		//fscanf(polyfile,"%lf ",&cQ[m]);
		//printf("%2.6lf ",cQ[m]);
	}
	fprintf(polyfile,"\n");
	for(m =0;m<=BIN_ORDER;m++)
	{
		fprintf(polyfile,"%2.6lf ",cU[m]);
		//fscanf(polyfile,"%lf ",&cU[m]);
		//printf("%2.6lf ",cU[m]);
	}
	fprintf(polyfile,"\n");
	for(m =0;m<=BIN_ORDER;m++)
	{
		fprintf(polyfile,"%2.6lf ",cV[m]);
		//fscanf(polyfile,"%lf ",&cV[m]);
		//printf("%2.6lf ",cV[m]);
	}
	fprintf(polyfile,"\n");

/*	fprintf(polyfile,"I ");
	jsd_print_poly(polyfile, cI, BIN_ORDER);
	fprintf(polyfile,"Q ");
	jsd_print_poly(polyfile, cQ, BIN_ORDER);
	fprintf(polyfile,"U ");
	jsd_print_poly(polyfile, cU, BIN_ORDER);
	fprintf(polyfile,"V ");
	jsd_print_poly(polyfile, cV, BIN_ORDER);
*///	fclose(polyfile);

/*	FILE * difffile;
	char difffilename[40];
	sprintf(binfilename,"%s/beam%d/polybindiffs%04i.dat",wappdata->daydata[day].mjd,day%7,chan);
	binfile = fopen(binfilename,"w");
	fprintf(binfile,"#DEC diffI diffQ diffU diffV\n");
	for(i=0;i<num_bins;i++)
	{
		double I,Q,U,V;
		I = binI[i] - jsd_poly_eval(binDEC[i], cI, BIN_ORDER);
		Q = binQ[i] - jsd_poly_eval(binDEC[i], cQ, BIN_ORDER);
		U = binU[i] - jsd_poly_eval(binDEC[i], cU, BIN_ORDER);
		V = binV[i] - jsd_poly_eval(binDEC[i], cV, BIN_ORDER);
		fprintf(binfile,"%2.6f %2.6f %2.6f %2.6f %2.6f\n",binDEC[i],I,Q,U,V);
	}
	fclose(binfile);
*/
//	jsd_denormalize(binDEC, num_bins, min, max);
	//modify the data to remove the dec effect according to the polyfit
//	for (d=0; d<wappdata->numDays; d++)
//	{
//		if (d%7!=beam) continue; //TODO: this is a workaround for non-beam aware datastructures
		FluxDayData * daydata;
//		printf("**********chan/2***%d:%f\n",chan/2,beamgain[0][chan/2]);
		daydata = &wappdata->daydata[day];
//	printf("in fitdecbin %s %d\n",daydata->mjd,day%7);
		for (r=0; r<daydata->numRecords; r++)
		{
//	            double DEC = daydata->records[r].DEC;
//	            i = (int) floor((DEC - decmin) / decgrain);
//        	    if(daydata->records[r].stokes.I > binI[i])
//          	    daydata->records[r].stokes.I -= binI[i];
//	            else
//         	    daydata->records[r].stokes.I = 0.0;
//	            daydata->records[r].stokes.Q -= binQ[i];
  //      	    daydata->records[r].stokes.U -= binU[i];
//        	    daydata->records[r].stokes.V -= binV[i];
			double DEC = NORMALIZE(daydata->records[r].DEC, min, max);
			daydata->records[r].stokes.I -= jsd_poly_eval(DEC, cI, BIN_ORDER);
			daydata->records[r].stokes.Q -= jsd_poly_eval(DEC, cQ, BIN_ORDER);
			daydata->records[r].stokes.U -= jsd_poly_eval(DEC, cU, BIN_ORDER);
			daydata->records[r].stokes.V -= jsd_poly_eval(DEC, cV, 5);

			if(day%7==1)
			{
			//daydata->records[r].stokes.Q-=daydata->records[r].stokes.I*Qleakage[1][chan/2];
			//daydata->records[r].stokes.U-=daydata->records[r].stokes.I*Uleakage[1][chan/2];
			//daydata->records[r].stokes.I/=beamgain[1][chan/2];
			daydata->records[r].stokes.I /=0.746917 ;
			}
			if(day%7==2)
			{
			//daydata->records[r].stokes.Q-=daydata->records[r].stokes.I*Qleakage[2][chan/2];
			//daydata->records[r].stokes.U-=daydata->records[r].stokes.I*Uleakage[2][chan/2];
			//daydata->records[r].stokes.I/=beamgain[2][chan/2];
			daydata->records[r].stokes.I /=0.667072 ;
			}
			if(day%7==3)
			{
			//daydata->records[r].stokes.Q-=daydata->records[r].stokes.I*Qleakage[3][chan/2];
			//daydata->records[r].stokes.U-=daydata->records[r].stokes.I*Uleakage[3][chan/2];
			//daydata->records[r].stokes.I/=beamgain[3][chan/2];
			daydata->records[r].stokes.I /=0.658215 ;
			}
			if(day%7==4)
			{
			//daydata->records[r].stokes.Q-=daydata->records[r].stokes.I*Qleakage[4][chan/2];
			//daydata->records[r].stokes.U-=daydata->records[r].stokes.I*Uleakage[4][chan/2];
			//daydata->records[r].stokes.I/=beamgain[4][chan/2];
			daydata->records[r].stokes.I /=0.641912 ;
			}
			if(day%7==5)
			{
			//daydata->records[r].stokes.Q-=daydata->records[r].stokes.I*Qleakage[5][chan/2];
			//daydata->records[r].stokes.U-=daydata->records[r].stokes.I*Uleakage[5][chan/2];
			//daydata->records[r].stokes.I/=beamgain[5][chan/2];
			daydata->records[r].stokes.I /=0.688649 ;
			}
			if(day%7==6)
			{
			//daydata->records[r].stokes.Q-=daydata->records[r].stokes.I*Qleakage[6][chan/2];
			//daydata->records[r].stokes.U-=daydata->records[r].stokes.I*Uleakage[6][chan/2];
			//daydata->records[r].stokes.I/=beamgain[6][chan/2];
			daydata->records[r].stokes.I /=0.650000 ;
			}
		}
//	}

	free(binDEC);
}

//hack
#define MAX_BIN_SIZE 8000

/*
 * decgrain - granularity of the declination bins, ie the bin width, in degrees
 */
static void day_dec_dependence(FluxWappData * wappdata, int day, float decmin, float decmax, float decgrain, int chan)
//static void beam_dec_dependence(FluxWappData * wappdata, int beam, float decmin, float decmax, float decgrain, int chan)
{
	int i;
	int num_bins;
	double **binarrayI;
	double **binarrayQ;
	double **binarrayU;
	double **binarrayV;
	double *cleanbinI;
	double *cleanbinQ;
	double *cleanbinU;
	double *cleanbinV;
	int * counts, *outliercounts;
	FILE * decfile;
	char filename[32+1];
	float nsigma = 1.5; //for the per bin data point exclusion
	num_bins = (int) ceil((decmax-decmin)/decgrain);
	decgrain = (decmax-decmin)/num_bins; //actual grain (bin) size
	binarrayI = (double**) malloc(num_bins * sizeof(double*));
	binarrayQ = (double**) malloc(num_bins * sizeof(double*));
	binarrayU = (double**) malloc(num_bins * sizeof(double*));
	binarrayV = (double**) malloc(num_bins * sizeof(double*));
	for (i=0; i<num_bins; i++)
	{
		binarrayI[i] = (double*) calloc(MAX_BIN_SIZE,sizeof(double));
		binarrayQ[i] = (double*) calloc(MAX_BIN_SIZE,sizeof(double));
		binarrayU[i] = (double*) calloc(MAX_BIN_SIZE,sizeof(double));
		binarrayV[i] = (double*) calloc(MAX_BIN_SIZE,sizeof(double));
	}
	counts = calloc(num_bins,sizeof(int));
	outliercounts = calloc(num_bins,sizeof(double));
	cleanbinI = calloc(num_bins,sizeof(double));
	cleanbinQ = calloc(num_bins,sizeof(double));
	cleanbinU = calloc(num_bins,sizeof(double));
	cleanbinV = calloc(num_bins,sizeof(double));

//	printf("Creating bins\n");
	//bin up the values in declination
//	create_dec_bins(wappdata, binarrayI, binarrayQ, binarrayU, binarrayV, counts, beam, decmin, decgrain, num_bins);

//	printf("Cleaning bins\n");
//	compute_clean_bins(binarrayI, counts, cleanbinI, num_bins, nsigma, outliercounts);
//	printf("Cleaning bins\n");
//	compute_clean_bins(binarrayQ, counts, cleanbinQ, num_bins, nsigma, outliercounts);
//	printf("Cleaning bins\n");
//	compute_clean_bins(binarrayU, counts, cleanbinU, num_bins, nsigma, outliercounts);
//	printf("Cleaning bins\n");
//	compute_clean_bins(binarrayV, counts, cleanbinV, num_bins, nsigma, outliercounts);

	//print bin data
/*	snprintf(filename, 32, "decbins_beam%i_chan%i.dat", beam, chan);
	decfile = fopen(filename, "w");
	fprintf(decfile, "#bin DEC counts cleanI cleanQ cleanU cleanV\n");
	for (i=0; i<num_bins; i++) {
		fprintf(decfile, "%i %f %i %g %g %g %g\n", i, i*decgrain+decmin, counts[i], cleanbinI[i], cleanbinQ[i], cleanbinU[i], cleanbinV[i]);
	}
	fclose(decfile);
*/
	//modify the data
//	printf("Fitting bins\n");
//	fit_dec_bins(wappdata, cleanbinI, cleanbinQ, cleanbinU, cleanbinV, num_bins, beam, decmin, decmax, decgrain);

	//day by day

	int d, r;

	memset(counts, 0, num_bins*sizeof(int));

	FluxDayData * daydata = &wappdata->daydata[day];

	if(daydata->numRecords == 0)
		return;
//	printf("in daydecdep %s %d\n",daydata->mjd,day%7);
	for (r=0; r<daydata->numRecords; r++)
	{
		double I = daydata->records[r].stokes.I;
		double Q = daydata->records[r].stokes.Q;
		double U = daydata->records[r].stokes.U;
		double V = daydata->records[r].stokes.V;
		double DEC = daydata->records[r].DEC;
		int bin = (int) floor((DEC - decmin) / decgrain);

		if (bin<0 || bin>=num_bins) continue;
		if (isfinite(I))
		{
			double *binarray;
			binarray = binarrayI[bin];
			binarray[counts[bin]] = I;
			binarray = binarrayQ[bin];
			binarray[counts[bin]] = Q;
			binarray = binarrayU[bin];
			binarray[counts[bin]] = U;
			binarray = binarrayV[bin];
			binarray[counts[bin]] = V;
			counts[bin]++;
		}
	}	

	int b;
	double mean, sigma;
	int outlier;
	double *binarray;
	int count;
	memset(outliercounts, 0, num_bins*sizeof(int));
	/*for (b=0; b<num_bins; b++)
	{
				//repeat until no more outliers
		do {
			outlier = 0;
			binarray = binarrayI[b];
			count = counts[b];
			mean = compute_mean(binarray, 0,count);
			sigma = compute_sigma(binarray, count, mean);

			//do the cleaning by setting outliers to NAN
			for (i=0; i<count; i++)
			{
				if (fabs(mean-binarray[i]) > nsigma*sigma)
				{
					binarray[i] = NAN;
					outlier = 1;
					outliercounts[b]++;
				}
			}

		} while (outlier);
		cleanbinI[b] = mean;
	}*/
///		printf("Cleaning bins day:%s beam%d I\n",daydata->mjd,day%7);
		compute_clean_bins(binarrayI, counts, cleanbinI, num_bins, nsigma, outliercounts);
//		printf("Cleaning bins Q\n");
		compute_clean_bins(binarrayQ, counts, cleanbinQ, num_bins, nsigma, outliercounts);
//		printf("Cleaning bins U\n");
		compute_clean_bins(binarrayU, counts, cleanbinU, num_bins, nsigma, outliercounts);
//		printf("Cleaning bins V\n");
		compute_clean_bins(binarrayV, counts, cleanbinV, num_bins, nsigma, outliercounts);


//			{
//				if (d%7!=beam) continue; //TODO: this is a workaround for non-beam aware datastructures
//				FluxDayData * daydata;

//				daydata = &wappdata->daydata[d];

	FILE * binfile;
	char binfilename[40];
	sprintf(binfilename,"%s/beam%d/decbins.dat",daydata->mjd,day%7);
	binfile = fopen(binfilename,"w");
	fprintf(binfile,"#DEC binI binQ binU binV\n");
	for(i=0;i<num_bins;i++)
	{
		fprintf(binfile,"%2.6f %2.6f %2.6f %2.6f %2.6f\n",decmin+i*decgrain+decgrain/2,cleanbinI[i],\
		cleanbinQ[i],cleanbinU[i],cleanbinV[i]);
	}
	fclose(binfile);

	fit_dec_bins(wappdata, cleanbinI, cleanbinQ, cleanbinU, cleanbinV, num_bins, day, decmin, decmax, decgrain,chan);

/*	for (r=0; r<daydata->numRecords; r++)
	{
		double DEC = daydata->records[r].DEC;
         	i = (int) floor((DEC - decmin) / decgrain);
 //     if(daydata->records[r].stokes.I > binI[i])
	//		daydata->records[r].stokes.I -= binI[i];
		//else
//	       	    daydata->records[r].stokes.I -= cleanbinI[i];
//		    daydata->records[r].stokes.Q -= cleanbinQ[i];
//		    daydata->records[r].stokes.U -= cleanbinU[i];
//		    daydata->records[r].stokes.V -= cleanbinV[i];
			double DEC = NORMALIZE(daydata->records[r].DEC, min, max);
			daydata->records[r].stokes.I -= jsd_poly_eval(DEC, cI, BIN_ORDER);
/			daydata->records[r].stokes.Q -= jsd_poly_eval(DEC, cQ, BIN_ORDER);
			daydata->records[r].stokes.U -= jsd_poly_eval(DEC, cU, BIN_ORDER);
			daydata->records[r].stokes.V -= jsd_poly_eval(DEC, cV, BIN_ORDER);
	}
*/
//day by day

	//do it again just so we can print out the bin data and see the residuals
/*	create_dec_bins(wappdata, binarrayI, binarrayQ, binarrayU, binarrayV, counts, beam, decmin, decgrain, num_bins);
	compute_clean_bins(binarrayQ, counts, cleanbinQ, num_bins, nsigma, outliercounts);
	compute_clean_bins(binarrayU, counts, cleanbinU, num_bins, nsigma, outliercounts);
	compute_clean_bins(binarrayV, counts, cleanbinV, num_bins, nsigma, outliercounts);
	compute_clean_bins(binarrayI, counts, cleanbinI, num_bins, nsigma, outliercounts);

	snprintf(filename, 32, "decbinsnew_beam%i_chan%i.dat", beam, chan);
	decfile = fopen(filename, "w");
	fprintf(decfile, "#bin DEC counts cleanI cleanQ cleanU cleanV\n");
	for (i=0; i<num_bins; i++) {
		fprintf(decfile, "%i %f %i %g %g %g %g\n", i, i*decgrain+decmin, counts[i], cleanbinI[i], cleanbinQ[i], cleanbinU[i], cleanbinV[i]);
	}
	fclose(decfile);
*/

	//cleanup
	for (i=0; i<num_bins; i++)
	{
		free(binarrayI[i]);
		free(binarrayQ[i]);
		free(binarrayU[i]);
		free(binarrayV[i]);
	}
	free(binarrayI);
	free(binarrayQ);
	free(binarrayU);
	free(binarrayV);
	free(cleanbinI);
	free(cleanbinQ);
	free(cleanbinU);
	free(cleanbinV);
	free(counts);
	free(outliercounts);
}


void calculate_dec_dependence(FluxWappData * wappdata, float decmin, float decmax, float decgrain, int chan)
{
//	int beam;
	int d;
//	for (beam=0; beam<7; beam++) {
//		beam_dec_dependence(wappdata, beam, decmin, decmax, decgrain, chan);
//	}
	for (d=0; d<wappdata->numDays; d++) {
		day_dec_dependence(wappdata, d, decmin, decmax, decgrain, chan);
	}
}

void average_dec_dependence(FluxWappData * wappdata, int lowchan, int highchan)
{
	int i,d,r,m;
	for (d=0; d<wappdata->numDays; d++) 
	{
		FluxDayData * daydata;
		double avgI[BIN_ORDER+1];
		double avgQ[BIN_ORDER+1];
		double avgU[BIN_ORDER+1];
		double avgV[BIN_ORDER+1];
		for(m =0;m<=BIN_ORDER;m++)
		{
			avgI[m]=0.0;
			avgQ[m]=0.0;
			avgU[m]=0.0;
			avgV[m]=0.0;
		}
		daydata = &wappdata->daydata[d];
		FILE *polyfile;
		char polyfilename[41];
		for(r = lowchan;r<=highchan;r++)
		{
//			printf("Reading polynomials chan %d\n",r);
			double cI[BIN_ORDER+1];
			double cQ[BIN_ORDER+1];
			double cU[BIN_ORDER+1];
			double cV[BIN_ORDER+1];
			sprintf(polyfilename,"%s/beam%d/polyfile%04i.dat",daydata->mjd,d%7,r);
			polyfile = fopen(polyfilename,"r");
			for(i = 0;i<=BIN_ORDER;i++)
			{
				fscanf(polyfile,"%lf ",&cI[i]);
				avgI[i] += cI[i];
//				printf("%f ",cI[i]);
			}	
//			printf("\n");
			for(i = 0;i<=BIN_ORDER;i++)
			{
				fscanf(polyfile,"%lf ",&cQ[i]);
				avgQ[i] += cQ[i];
//				printf("%f ",cQ[i]);
			}	
//			printf("\n");
			for(i = 0;i<=BIN_ORDER;i++)
			{
				fscanf(polyfile,"%lf ",&cU[i]);
				avgU[i] += cU[i];
//				printf("%f ",cU[i]);
			}	
//			printf("\n");
			for(i = 0;i<=BIN_ORDER;i++)
			{
				fscanf(polyfile,"%lf ",&cV[i]);
				avgV[i] += cV[i];
//				printf("%f ",cV[i]);
			}	
//			printf("\n");
			fclose(polyfile);
		}
		sprintf(polyfilename,"%s/beam%d/polyfile%04i.dat",daydata->mjd,d%7,1);
		polyfile = fopen(polyfilename,"w");
//		printf("Writing avg polynomials\n",r);
		for(m =0;m<=BIN_ORDER;m++)
		{
			fprintf(polyfile,"%2.6lf ",avgI[m]/(highchan-lowchan+1));
		}
		fprintf(polyfile,"\n");
		for(m =0;m<=BIN_ORDER;m++)
		{
			fprintf(polyfile,"%2.6lf ",avgQ[m]/(highchan-lowchan+1));
		}
		fprintf(polyfile,"\n");
		for(m =0;m<=BIN_ORDER;m++)
		{
			fprintf(polyfile,"%2.6lf ",avgU[m]/(highchan-lowchan+1));
		}
		fprintf(polyfile,"\n");
		for(m =0;m<=BIN_ORDER;m++)
		{
			fprintf(polyfile,"%2.6lf ",avgV[m]/(highchan-lowchan+1));
		}
		fprintf(polyfile,"\n");
		fclose(polyfile);		
	}
}	

void remove_dec_dependence(FluxWappData * wappdata, float decmin, float decmax, int chan)
{
	int i,d,r;
	for (d=0; d<wappdata->numDays; d++) 
	{
		FluxDayData * daydata;
		daydata = &wappdata->daydata[d];
		double cI[BIN_ORDER+1];
		double cQ[BIN_ORDER+1];
		double cU[BIN_ORDER+1];
		double cV[BIN_ORDER+1];
		FILE *polyfile;
		char polyfilename[41];
		sprintf(polyfilename,"%s/beam%d/polyfile%04i.dat",daydata->mjd,d%7,chan);
		polyfile = fopen(polyfilename,"r");
		for(i = 0;i<=BIN_ORDER;i++)
			fscanf(polyfile,"%lf ",&cI[i]);	
		for(i = 0;i<=BIN_ORDER;i++)
			fscanf(polyfile,"%lf ",&cQ[i]);	
		for(i = 0;i<=BIN_ORDER;i++)
			fscanf(polyfile,"%lf ",&cU[i]);	
		for(i = 0;i<=BIN_ORDER;i++)
			fscanf(polyfile,"%lf ",&cV[i]);	
		fclose(polyfile);
		for (r=0; r<daydata->numRecords; r++)
		{
			//double DEC = daydata->records[r].DEC;
			double DEC = NORMALIZE(daydata->records[r].DEC, decmin, decmax);
			daydata->records[r].stokes.I -= jsd_poly_eval(DEC, cI, BIN_ORDER);
			daydata->records[r].stokes.Q -= jsd_poly_eval(DEC, cQ, BIN_ORDER);
			daydata->records[r].stokes.U -= jsd_poly_eval(DEC, cU, BIN_ORDER);
			daydata->records[r].stokes.V -= jsd_poly_eval(DEC, cV, BIN_ORDER);
		}
	}
}

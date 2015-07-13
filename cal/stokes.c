#include <math.h>
#include "stokes.h"
#include "rfi.h"
#include "cal.h"
#include "fluxdata.h"
#include "denoising.h"
#include "chebyshev.h"
#include "errno.h"//SSG
//----------------------------------------------------------------------------------------------------------------------------------------
extern int errno; //SSG
extern int multibeam; //SSG
//----------------------------------------------------------------------------------------------------------------------------------------
typedef struct 
	{
    float x[MAX_CHANNELS];
    float y[MAX_CHANNELS];
    float phi[MAX_CHANNELS];
	}GainSet;

typedef struct 
	{
	float xx;
	float yy;
	float xy;
	float yx;
	}PolParams;
//----------------------------------------------------------------------------------------------------------------------------------------
 
// Computes the gains values for x, y and phi using the cal values and supplied Tcal corrections.
// The Tcal values are expected to range from 0..highchan-lowchan instead of lowchan..highchan.

static void compute_gains(const SpecRecord * pRec, GainSet * pGain,
	const int lowchan, const int highchan, const float Tcalx[],
	const float Tcaly[])
{
	int i;
	double calU, calV;
	for (i=lowchan; i<highchan; i++)
		{
		pGain->x[i] = sqrt(Tcalx[i]/pRec->cal.xx[i]);
		pGain->y[i] = sqrt(Tcaly[i]/pRec->cal.yy[i]);
		}
}
//----------------------------------------------------------------------------------------------------------------------------------------
static void compute_phases(const SpecRecord * pRec, GainSet * pGain,
	const int lowchan, const int highchan, const float Tcalx[],
	const float Tcaly[])
{
	int i;
	double calU, calV;

	for (i=lowchan; i<highchan; i++)
	{
		calU = pRec->cal.xy[i] + pRec->cal.yx[i];
		calV = pRec->cal.xy[i] - pRec->cal.yx[i];
		pGain->phi[i] = atan2(calV, calU);
		if(i >= lowchan+1 && pGain->phi[i-1] > M_PI/2.0 && pGain->phi[i] < -0.0)
		{
			//printf("inside if i is %d\n",i); 
			pGain->phi[i] += 2*M_PI;
		}
		else if(i >= lowchan+1 && pGain->phi[i-1] < -M_PI/2.0 && pGain->phi[i] > 0.0)
		{
			//printf("inside else if i is %d\n",i); 
			pGain->phi[i] -= 2*M_PI;
		}

	}
}

// Uses obsStokes to calibrate the calStokes

static void calibrate_stokes(StokesSet * calStokes, const GainSet * gain,
	const StokesSet * obsStokes, const int lowchan, const int highchan)
{
	int i;

	for (i=lowchan; i<highchan; i++)
		{
		calStokes->I[i] =  ((gain->x[i]*gain->x[i] + gain->y[i]*gain->y[i])/4.0) * obsStokes->I[i] + ((gain->x[i]*gain->x[i] - gain->y[i]*gain->y[i])/4.0) * obsStokes->Q[i];

		calStokes->Q[i] =  ((gain->x[i]*gain->x[i] - gain->y[i]*gain->y[i])/4.0) * obsStokes->I[i] + ((gain->x[i]*gain->x[i] + gain->y[i]*gain->y[i])/4.0) * obsStokes->Q[i];


		//A factor of 2 might be missing in these calculations below giving bigger U and V that are fixed later. Need
		// to confirm. The original assumption that the spectrometer was spitting values a factor of 2 too large
		// may be incorrect. nonetheless we are dividing by 2 later and hence the output is all good. but it may
		//be better to fix it here if the problem is actually here.
		calStokes->U[i] = (gain->x[i]*gain->y[i]) * (cos(gain->phi[i]) * obsStokes->U[i] + sin(gain->phi[i])*obsStokes->V[i]);

		calStokes->V[i] = (gain->x[i]*gain->y[i]) * (-sin(gain->phi[i]) * obsStokes->U[i] + cos(gain->phi[i])*obsStokes->V[i]);
		}
}
//----------------------------------------------------------------------------------------------------------------------------------------

static void print_stokes(FILE * file, const StokesSet * obs,
	const StokesSet * cal, const GainSet * gain, const int lowchan,
	const int highchan)
{
	int i;
	for (i=lowchan; i<highchan; i++)
		{
		fprintf(file, "%4d %7.4f %7.4f %7.4f %7.4f %9.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f\n",
				i,      
				obs->I[i], obs->Q[i], obs->U[i], obs->V[i], cal->I[i], cal->Q[i], cal->U[i], cal->V[i],
				gain->x[i], gain->y[i], gain->phi[i]);
		}
}
//----------------------------------------------------------------------------------------------------------------------------------------
static void print_gain(FILE * filex, FILE * filey, FILE * filep,
	const SpecRecord * pRec, const GainSet * gain, const int lowchan,
	const int highchan)
{
	int i;
	fprintf(filex, "%f ", pRec->AST);
	fprintf(filey, "%f ", pRec->AST);
	fprintf(filep, "%f ", pRec->AST);
	for (i=lowchan; i<highchan; i++)
		{
		fprintf(filex, "%f ", gain->x[i]);
		fprintf(filey, "%f ", gain->y[i]);
		fprintf(filep, "%f ", gain->x[i]);
		}
	fprintf(filex, "\n");
	fprintf(filey, "\n");
	fprintf(filep, "\n");
	
}
//----------------------------------------------------------------------------------------------------------------------------------------
// NOTE: The ObsCal uses the raw calon-caloff, not the smoothed cal field of the record
// The ObsSky uses caloff + calon - smoothed cal.

static void compute_observed_stokes(const SpecRecord * pRec, StokesSet * ObsCal,
	StokesSet * ObsSky, const int lowchan, const int highchan)
{
	int i;
	PolParams cal; //calculated values of the cal
	PolParams data; //calculated values of the cal

	for (i=lowchan; i<highchan; i++)
		{
/*		cal.xx = pRec->calon.xx[i] - pRec->caloff.xx[i];
		cal.yy = pRec->calon.yy[i] - pRec->caloff.yy[i];
		cal.xy = pRec->calon.xy[i] - pRec->caloff.xy[i];
		cal.yx = pRec->calon.yx[i] - pRec->caloff.yx[i];

		ObsCal->I[i] = cal.xx + cal.yy;
		ObsCal->Q[i] = cal.xx - cal.yy;      
		ObsCal->U[i] = cal.xy + cal.yx;
		ObsCal->V[i] = cal.xy - cal.yx;*/

		//Use smoothed cal
                ObsCal->I[i] = pRec->cal.xx[i] + pRec->cal.yy[i];
                ObsCal->Q[i] = pRec->cal.xx[i] - pRec->cal.yy[i];
                ObsCal->U[i] = pRec->cal.xy[i] + pRec->cal.yx[i];
                ObsCal->V[i] = pRec->cal.xy[i] - pRec->cal.yx[i];

	 	data.xx = pRec->caloff.xx[i] + (pRec->calon.xx[i] - pRec->cal.xx[i]);
	 	data.xy = pRec->caloff.xy[i] + (pRec->calon.xy[i] - pRec->cal.xy[i]);
	 	data.yx = pRec->caloff.yx[i] + (pRec->calon.yx[i] - pRec->cal.yx[i]);
	 	data.yy = pRec->caloff.yy[i] + (pRec->calon.yy[i] - pRec->cal.yy[i]);

		ObsSky->I[i] = data.xx + data.yy;
		ObsSky->Q[i] = data.xx - data.yy;
		ObsSky->U[i] = data.xy + data.yx;
		ObsSky->V[i] = data.xy - data.yx;
		/*ObsSky->I[i] = data.xx;
		ObsSky->Q[i] = data.yy;
		ObsSky->U[i] = data.xy;
		ObsSky->V[i] = data.yx;
*/
		}

}
//----------------------------------------------------------------------------------------------------------------------------------------

// Loads a particular datapoint pRec with data from the 
// provided TrueSky stokes set.  RFI flags are respected if ignoreRFI is 0.
// NAN is loaded into the pRec if a value cannot be computed.
// Average stokes values are packed into channel 0.

static void compute_final_stokes(SpecRecord * pRec, const StokesSet * TrueSky,
	const int RFIF, const int lowchan, const int highchan)
{
	int i;

//	FILE *f = fopen("linfitphase.dat","w");
	//FILE *f = fopen("smoothphase.dat","w");

	for(i=lowchan; i<=highchan; i++)
	{
		if((!RFIF || pRec->flagRFI[i] == RFI_NONE) && isfinite(TrueSky->I[i]))
		{
			pRec->stokes.I[i] = TrueSky->I[i];
			pRec->stokes.Q[i] = TrueSky->Q[i];
			pRec->stokes.U[i] = TrueSky->U[i];
			pRec->stokes.V[i] = TrueSky->V[i];
		}
		else
		{
			pRec->stokes.I[i] = NAN;
			pRec->stokes.Q[i] = NAN;
			pRec->stokes.U[i] = NAN;
			pRec->stokes.V[i] = NAN;
		}
//		fprintf(f,"%04i %2.6f %2.6f %2.6f %2.6f\n",i,pRec->stokes.I[i],pRec->stokes.Q[i],pRec->stokes.U[i],pRec->stokes.V[i]);
	}
//	fclose(f);
	//exit(0);
	//printf("Chan 1000 Stokes I %f\n",pRec->stokes.I[1000]);
}
//----------------------------------------------------------------------------------------------------------------------------------------

// Does a straight average of the channel data in the given channel range.
// Results are packed into the 'channel 0' records for convienence.
// This allows a single channel to represent the average data later in the pipeline.

void average_stokes(SpecRecord dataset[], const int size, const int lowchan,
	const int highchan, const float hidrogenfreq, const float hidrogenband,
	const float freq[])
{
	float minf=hidrogenfreq-hidrogenband, maxf=hidrogenfreq+hidrogenband;
	int i, chan;
	float I, Q, U, V;
	int count;
	SpecRecord * pRec;

	//iterate over each time step
	for(i=0; i<size; i++) 
		{
		pRec = &(dataset[i]);
		if(pRec->flagBAD) 
			{
			//put in NAN if the record is bad
			pRec->stokes.I[0] = NAN;
			pRec->stokes.Q[0] = NAN;
			pRec->stokes.U[0] = NAN;
			pRec->stokes.V[0] = NAN;
			continue;
			}

		count = 0;
		I = Q = U = V = 0.0;
		for(chan=lowchan; chan<=highchan; chan++)
			{
			if(isfinite(pRec->stokes.I[chan]) && (freq[chan] < minf || freq[chan] > maxf))
				{
				I += pRec->stokes.I[chan];
				Q += pRec->stokes.Q[chan];
				U += pRec->stokes.U[chan];
				V += pRec->stokes.V[chan];
				count++;
				}
			}

		//pack the averages into channel 0
		pRec->stokes.I[0] = I/count;
		pRec->stokes.Q[0] = Q/count;
		pRec->stokes.U[0] = U/count;
		pRec->stokes.V[0] = V/count;
	}
}
//----------------------------------------------------------------------------------------------------------------------------------------

//external function to calculate the stokes parameters

void calculate_stokes(SpecRecord dataset[], const int size, const int lowchan,
	const int highchan, const int RFIF, const int calskyfiles,
	const float Tcalx[], const float Tcaly[], const int uvDenoising,
	const float uvDenoisingTau, const float uvDenoisingLambda, const int start,
	const int end)
{
	int i, n;

	SpecRecord * pRec;

	StokesSet CalCal;
	StokesSet ObsCal;
	StokesSet TrueSky;
	StokesSet ObsSky;
	FILE * skyfile;
	FILE * calfile;
	GainSet gain;

	if(calskyfiles)
	{
		printf("[WARNING] Writing cal/sky files not supported by openmp branch.");
	}

	//Uses average phases over all time datapoints for better S/N
	float av_phi[MAX_CHANNELS];
	for(n=lowchan; n<highchan; n++)
	{
		av_phi[n] = 0.0;
	}

	int count = 0;
	for (i=0; i<size; i++)
	{
		pRec = &(dataset[i]);
		if (pRec->flagBAD)
			continue;

		compute_phases(pRec, &gain, lowchan, highchan, Tcalx, Tcaly);

		for(n=lowchan; n<highchan; n++)
		{
			av_phi[n] += gain.phi[n];
		}
		count++;
	}

	for(n=lowchan; n<highchan; n++)
	{
		av_phi[n] /= count;
		gain.phi[n] = av_phi[n];
	}

	FILE *err = fopen("phase.dat","w");

	if(uvDenoising)
	{
		for(n=lowchan; n<highchan; n++)
		{
			fprintf(err,"%04i %2.6f\n",n,gain.phi[n]);
		}
	}

	fclose(err);

	#pragma omp parallel for \
		firstprivate(gain) \
		private(i, ObsCal, ObsSky, TrueSky, CalCal) \
		shared(dataset, Tcalx, Tcaly) \
		schedule(static)

	for(i=start; i<end; i++)
	{
		pRec = &(dataset[i]);
		if (pRec->flagBAD)
			continue;
		compute_observed_stokes(pRec, &ObsCal, &ObsSky, lowchan, highchan);
		compute_gains(pRec, &gain, lowchan, highchan, Tcalx, Tcaly);
		calibrate_stokes(&TrueSky, &gain, &ObsSky, lowchan, highchan);
		compute_final_stokes(pRec, &TrueSky, RFIF, lowchan, highchan);
	}
}
//----------------------------------------------------------------------------------------------------------------------------------------

void correct_beamgains(SpecRecord dataset[], const int size, const int lowchan,
	const int highchan, const int beam)
{
	int n;
	double I, Q, U, V;
	char filename[32+1];
	int chan;
	int found;
	char header[81];
	int readchan; 
	FILE * beamgainfile;
	float gain[7];
	
	beamgainfile = fopen("../../beamgains.dat","r");
	if(beamgainfile == NULL)
		{
		printf("ERROR: unable to open gain file %s\n", filename);
		return;
		}
		else fgets(header,80,beamgainfile);
	found = 0;
	for(chan=lowchan; chan<highchan; chan++) 
		{
		while(!found)//TODO : algo for searching not very elegant yet.
			{
			fscanf(beamgainfile,"%d",&readchan);
			if(readchan == chan)
				{
				fscanf(beamgainfile,"%f %f %f %f %f %f %f",&gain[0],&gain[1],&gain[2],&gain[3],&gain[4],&gain[5],&gain[6]);
				found = 1;
				}
			else
				{
				fgets(header,80,beamgainfile);
				}
			}
		found = 0;
		for(n=0; n<size; n++)
			{
			SpecRecord * pRec = &(dataset[n]);
			if (pRec->flagBAD || pRec->flagRFI[chan] != RFI_NONE) 
				{
				I = Q = U = V = NAN;
				}
				else 
					{
					I = pRec->stokes.I[chan]/gain[beam];
					Q = pRec->stokes.Q[chan];
					U = pRec->stokes.U[chan];
					V = pRec->stokes.V[chan];
					}
			} 
		}
	fclose(beamgainfile);
}
//----------------------------------------------------------------------------------------------------------------------------------------
/*
void write_channel_data(SpecRecord dataset[], int size, int lowchan, int highchan)
{
	int n, i;
	double I, Q, U, V; 
	FILE *fluxfile;
	char filename[32+1];
	int chan;
	
	for(chan=lowchan; chan<highchan; chan++) 
		{
		char *buffer = (char*)calloc(size*14*32, sizeof(char));
		i = sprintf(buffer, "# RA DEC AST I Q U V\n"); 
		for(n=0; n<size; n++)
			{
			if(dataset[n].flagBAD || dataset[n].flagRFI[chan] != RFI_NONE) 
				{
				I = Q = U = V = NAN;
				} 
				else 
					{
					I = dataset[n].stokes.I[chan];
					Q = dataset[n].stokes.Q[chan];
					U = dataset[n].stokes.U[chan];
					V = dataset[n].stokes.V[chan];
					}
			i +=sprintf(buffer + i,"%2.8f %2.8f %8.2f %4.6f %4.6f %4.6f %4.6f\n", dataset[n].RA, dataset[n].DEC, dataset[n].AST, I, Q, U, V);
			}
		snprintf(filename, 32, "fluxtime%04i.dat", chan);
		fluxfile = fopen(filename, "wb");
		if(fluxfile == NULL)
			{
			printf("ERROR: unable to open file %s\n", filename);
			printf("DIAGNOSTIC: errno %d\n",errno);//SSG
			return;
			}			
		fprintf(fluxfile, "%s\n", buffer);
		fclose(fluxfile);
		free(buffer);
		printf("%f\r", (chan - lowchan + 1)*100.0/(highchan - lowchan));
		}
	printf("\n");
}
*/
//----------------------------------------------------------------------------------------------------------------------------------------
void write_channel_data(const SpecRecord dataset[], const int size,
	const int lowchan, const int highchan)
{
	int n;
	double I, Q, U, V; 
	FILE *fluxfile;
	char filename[32+1];
	int chan;
	
	for(chan=lowchan; chan<highchan; chan++) 
		{
		snprintf(filename, 32, "fluxtime%04i.dat", chan);
		fluxfile = fopen(filename, "w");
		if(fluxfile == NULL)
			{
			printf("ERROR: unable to open file %s\n", filename);
			printf("DIAGNOSTIC: errno %d\n",errno);//SSG
			return;
			}
		fprintf(fluxfile, "# RA DEC AST I Q U V\n");
		for(n=0; n<size; n++)
			{
			if(dataset[n].flagBAD || dataset[n].flagRFI[chan] != RFI_NONE) 
				{
				I = Q = U = V = NAN;
				} 
				else 
					{
					I = dataset[n].stokes.I[chan];
					Q = dataset[n].stokes.Q[chan];
					U = dataset[n].stokes.U[chan]/2.0; // UV correction
					V = dataset[n].stokes.V[chan]/2.0; // UV correction
					}
			fprintf(fluxfile,"%2.8f %2.8f %8.2f %4.6f %4.6f %4.6f %4.6f\n", dataset[n].RA, dataset[n].DEC, dataset[n].AST, I, Q, U, V);
			}
		fclose(fluxfile);
		printf("%f\r", (chan - lowchan + 1)*100.0/(highchan - lowchan));
		}
	printf("\n");
}
//----------------------------------------------------------------------------------------------------------------------------------------
void write_binary_channel_data(const SpecRecord dataset[], const int size,
	const int lowchan, const int highchan)
{
	int n;
	float I, Q, U, V, fRA, fDEC, fAST; 
	FILE *fluxfile;
	char filename[32+1];
	int chan;
	
	for(chan=lowchan; chan<highchan; chan++) 
		{
		// average channel is still written out the old way
		// but lets name it something better
		if( highchan == 1 && lowchan == 0 ) {
			snprintf(filename, 32, "average.dat");
		} else {
			// if we get here, we aren't calling the new single file code method
			snprintf(filename, 32, "fluxtime%04i.dat", chan);
		}
		fluxfile = fopen(filename, "wb");
		if(fluxfile == NULL)
			{
			printf("ERROR: unable to open file %s\n", filename);
			printf("DIAGNOSTIC: errno %d\n",errno);//SSG
			return;
			}
		fwrite(&size, sizeof(int), 1, fluxfile);
		for(n=0; n<size; n++)
			{
			fRA = dataset[n].RA; 
			fDEC = dataset[n].DEC;
			fAST = dataset[n].AST;
			if(dataset[n].flagBAD || dataset[n].flagRFI[chan] != RFI_NONE) 
				{
				I = Q = U = V = NAN;
				} 
				else 
					{
					I = dataset[n].stokes.I[chan];
					Q = dataset[n].stokes.Q[chan];
					U = dataset[n].stokes.U[chan]/2.0; // UV correction
					V = dataset[n].stokes.V[chan]/2.0; // UV correction
					}
			fwrite(&fRA, sizeof(float), 1, fluxfile);
			fwrite(&fDEC, sizeof(float), 1, fluxfile);
			fwrite(&fAST, sizeof(float), 1, fluxfile);
			fwrite(&I, sizeof(float), 1, fluxfile);
			fwrite(&Q, sizeof(float), 1, fluxfile);
			fwrite(&U, sizeof(float), 1, fluxfile);
			fwrite(&V, sizeof(float), 1, fluxfile);
			}
		fclose(fluxfile);
		//printf("%f\r", (chan - lowchan + 1)*100.0/(highchan - lowchan));
		}
	printf("\n");
}
//----------------------------------------------------------------------------------------------------------------------------------------
void write_binary_channel_data_single_file(const SpecRecord dataset[],
	const int numRecords, const int lowchan, const int highchan)
{
	int n;
	float I, Q, U, V, fRA, fDEC, fAST;
	FILE *fluxfile;
	FILE *fluxconfig;
	char filename[32+1];
	char configfilename[32+1];
	int chan;
	int startRecord;

	snprintf(filename, 32, "fluxtime.dat");
	fluxfile = fopen(filename, "wb");
	if(fluxfile == NULL)
		{
		printf("ERROR: unable to open file %s\n", filename);
		printf("DIAGNOSTIC: errno %d\n",errno);//SSG
		return;
		}

	char buf[4096];
	setbuffer(fluxfile,buf,4096);

	snprintf(configfilename, 32, "fluxtime.dat_cfg");
	fluxconfig = fopen(configfilename, "wb");
	if(fluxconfig == NULL)
	{
		printf("ERROR: unable to open file %s\n", filename);
		printf("DIAGNOSTIC: errno %d\n",errno);//SSG
		return;
	}

	startRecord = 0;

	for(chan=lowchan; chan<highchan; chan++)
	{
		// write config file entry
		fprintf(fluxconfig, "%d %d %d\n", chan, startRecord, numRecords);

		for(n=0; n<numRecords; n++)
		{
			fRA = dataset[n].RA;
			fDEC = dataset[n].DEC;
			fAST = dataset[n].AST;

			if(dataset[n].flagBAD || dataset[n].flagRFI[chan] != RFI_NONE)
			{
				I = NAN;
				Q = NAN;
				U = NAN;
				V = NAN;
			}
			else
			{
				I = dataset[n].stokes.I[chan];
				Q = dataset[n].stokes.Q[chan];
				U = dataset[n].stokes.U[chan]/2.0; // UV correction
				V = dataset[n].stokes.V[chan]/2.0; // UV correction
			}
			fwrite(&fRA, sizeof(float), 1, fluxfile);
			fwrite(&fDEC, sizeof(float), 1, fluxfile);
			fwrite(&fAST, sizeof(float), 1, fluxfile);
			fwrite(&I, sizeof(float), 1, fluxfile);
			fwrite(&Q, sizeof(float), 1, fluxfile);
			fwrite(&U, sizeof(float), 1, fluxfile);
			fwrite(&V, sizeof(float), 1, fluxfile);
			startRecord++;
		}
	//printf("%f\r", (chan - lowchan + 1)*100.0/(highchan - lowchan));
	}
	fclose(fluxfile);
	fclose(fluxconfig);
	printf("\n");
}

void write_rfi_data(const SpecRecord dataset[], const int numRecords,
	const int lowchan, const int highchan )
{
	FILE *freqRFI = fopen( "rfifrq.dat", "w");
	FILE *timeRFI = fopen( "rfitime.ann", "w");

	if (freqRFI == NULL || timeRFI == NULL ) {
		printf("ERROR: unable to open rfi files\n");
		return;
	}

	int chan, n, i, freqCount=0, timeCount=0;
	float fRA, fDEC;

	fprintf( timeRFI, "COLOUR GREEN\n");

	for(n=0; n<numRecords; n++)
	{
		if( dataset[n].flagBAD )
		{
			fRA = dataset[n].RA;
			fDEC = dataset[n].DEC;

			if( dataset[n].flagRFI[lowchan] == RFI_OUTOFBAND )
			{
				timeCount++;
				fprintf( timeRFI, "CROSS %f %f 0.2 0.2\n", fRA, fDEC );
			}

		}

		for(chan=lowchan; chan<highchan; chan++)
		{
//			if( dataset[n].flagRFI[chan] ==  RFI_CALOFF_XX ||
//					dataset[n].flagRFI[chan] ==  RFI_CALOFF_XY ||
//					dataset[n].flagRFI[chan] ==  RFI_CALOFF_YX  ||
//					dataset[n].flagRFI[chan] ==  RFI_CALOFF_YY )
			if( dataset[n].flagRFI[chan] != RFI_NONE)
			{
				fprintf( freqRFI, "%d %d %f %f\n", n, chan, dataset[n].RA, dataset[n].DEC );
				freqCount++;
			}
		}
	}

	//printf("After, found freq  = %d and time = %d", freqCount, timeCount);

	fclose( freqRFI );
	fclose( timeRFI );
}

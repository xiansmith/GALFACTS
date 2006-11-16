#include "stokes.h"
#include "rfi.h"
#include "cal.h"
#include <math.h>
#include "fluxdata.h"
#include "errno.h"//SSG

extern int errno; //SSG
extern int multibeam; //SSG


typedef struct {
    float x[MAX_CHANNELS];
    float y[MAX_CHANNELS];
    float phi[MAX_CHANNELS];
} GainSet;



typedef struct {
	float xx;
	float yy;
	float xy;
	float yx;
} PolParams;



/* 
 * Computes the gains values for x, y and phi using the cal values and supplied Tcal corrections.
 * The Tcal values are expected to range from 0..highchan-lowchan instead of lowchan..highchan.
 */
static void compute_gains(const SpecRecord * pRec, GainSet * pGain, int lowchan, int highchan, float Tcalx[], float Tcaly[])
{
	int i;
	double calU, calV;

	for (i=lowchan; i<=highchan; i++)
	{
		pGain->x[i] = sqrt(Tcalx[i]/pRec->cal.xx[i]);
		pGain->y[i] = sqrt(Tcaly[i]/pRec->cal.yy[i]);
		calU = pRec->cal.xy[i] + pRec->cal.yx[i];
		calV = pRec->cal.xy[i] - pRec->cal.yx[i];
		pGain->phi[i] = atan2(calV, calU);

		//if (!finite(pGain->x[i]))
		//	printf("%8.2f %i %6.4f %6.4f\n", pRec->AST, i, pRec->cal.xx[i], pRec->cal.yy[i]);
	}
}



//uses obsStokes to calibrate the calStokes
static void calibrate_stokes(StokesSet * calStokes, const GainSet * gain, const StokesSet * obsStokes, int lowchan, int highchan)
{
	int i;

	// calibrate the cal signal
	for (i=lowchan; i<=highchan; i++)
	{
		calStokes->I[i] =  ((gain->x[i]*gain->x[i] + gain->y[i]*gain->y[i])/4.0) * obsStokes->I[i] 
			+ ((gain->x[i]*gain->x[i] - gain->y[i]*gain->y[i])/4.0) * obsStokes->Q[i];

		calStokes->Q[i] =  ((gain->x[i]*gain->x[i] - gain->y[i]*gain->y[i])/4.0) * obsStokes->I[i] 
			+ ((gain->x[i]*gain->x[i] + gain->y[i]*gain->y[i])/4.0) * obsStokes->Q[i];

		calStokes->U[i] = (gain->x[i]*gain->y[i]) * (cos(gain->phi[i]) * obsStokes->U[i]
			+ sin(gain->phi[i])*obsStokes->V[i]);

		calStokes->V[i] = (gain->x[i]*gain->y[i]) * (-sin(gain->phi[i]) * obsStokes->U[i]
			+cos(gain->phi[i])*obsStokes->V[i]);
	}
}


static void print_stokes(FILE * file, const StokesSet * obs, const StokesSet * cal, const GainSet * gain, int lowchan, int highchan)
{
	int i;
	for (i=lowchan; i<=highchan; i++)
	{
		fprintf(file, "%4d %7.4f %7.4f %7.4f %7.4f %9.4f %7.4f %7.4f %7.4f %8.3f %6.3f %6.3f\n",
				i,      
				obs->I[i], obs->Q[i], obs->U[i], obs->V[i], cal->I[i], cal->Q[i], cal->U[i], cal->V[i],
				gain->x[i], gain->y[i], gain->phi[i]);
	}
	fprintf(file, "\n");
}


//NOTE: The ObsCal uses the raw calon-caloff, not the smoothed cal field of the record
//The ObsSky uses caloff + calon - smoothed cal.
static void compute_observed_stokes(const SpecRecord * pRec, StokesSet * ObsCal, StokesSet * ObsSky)
{
	int i;
	PolParams cal; //calculated values of the cal
	PolParams data; //calculated values of the cal

	for (i=0; i<MAX_CHANNELS; i++)
	{
		cal.xx = pRec->calon.xx[i] - pRec->caloff.xx[i];
		cal.yy = pRec->calon.yy[i] - pRec->caloff.yy[i];
		cal.xy = pRec->calon.xy[i] - pRec->caloff.xy[i];
		cal.yx = pRec->calon.yx[i] - pRec->caloff.yx[i];

		ObsCal->I[i] = cal.xx + cal.yy;
		ObsCal->Q[i] = cal.xx - cal.yy;      
		ObsCal->U[i] = cal.xy + cal.yx;
		ObsCal->V[i] = cal.xy - cal.yx;

	 	data.xx = pRec->caloff.xx[i] + (pRec->calon.xx[i] - pRec->cal.xx[i]);
	 	data.xy = pRec->caloff.xy[i] + (pRec->calon.xy[i] - pRec->cal.xy[i]);
	 	data.yx = pRec->caloff.yx[i] + (pRec->calon.yx[i] - pRec->cal.yx[i]);
	 	data.yy = pRec->caloff.yy[i] + (pRec->calon.yy[i] - pRec->cal.yy[i]);

		ObsSky->I[i] = data.xx + data.yy;
		ObsSky->Q[i] = data.xx - data.yy;
		ObsSky->U[i] = data.xy + data.yx;
		ObsSky->V[i] = data.xy - data.yx;
	}
}

/* Loads a particular datapoint pRec with data from the 
 * provided TrueSky stokes set.  RFI flags are respected if ignoreRFI is 0.
 * NAN is loaded into the pRec if a value cannot be computed.
 * Average stokes values are packed into channel 0.
 */
static void compute_final_stokes(SpecRecord * pRec, StokesSet * TrueSky, int ignoreRFI, int lowchan, int highchan)
{
	int i;
    int count;
	double I, Q, U, V;

	count = 0;
	I = Q = U = V = 0.0;
	for (i=lowchan; i<=highchan; i++)
	{
		if ((ignoreRFI || pRec->flagRFI[i] == RFI_NONE) && isfinite(TrueSky->I[i]))
		{
			pRec->stokes.I[i] = TrueSky->I[i];
			pRec->stokes.Q[i] = TrueSky->Q[i];
			pRec->stokes.U[i] = TrueSky->U[i];
			pRec->stokes.V[i] = TrueSky->V[i];
			I += TrueSky->I[i];
			Q += TrueSky->Q[i];
			U += TrueSky->U[i];
			V += TrueSky->V[i];
			count++;
		}
		else
		{
			pRec->stokes.I[i] = NAN;
			pRec->stokes.Q[i] = NAN;
			pRec->stokes.U[i] = NAN;
			pRec->stokes.V[i] = NAN;
		}
	}
	//pack the averages into channel 0
	pRec->stokes.I[0] = I/count;
	pRec->stokes.Q[0] = Q/count;
	pRec->stokes.U[0] = U/count;
	pRec->stokes.V[0] = V/count;
}


//external function to calculate the stokes parameters
void calculate_stokes(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI, float Tcalx[], float Tcaly[])
{
	int i;

	SpecRecord * pRec;

	StokesSet CalCal;
	StokesSet ObsCal;
	StokesSet TrueSky;
	StokesSet ObsSky;

	GainSet gain;

	FILE * skyfile;
	FILE * calfile;
	FILE * gainfile;
/*SSG
	calfile = fopen("calfile.dat", "w");
	fprintf(calfile, "# Chan ObsCalI ObsCalQ ObsCalU ObsCalV CalCalI CalCalQ CalCalU CalCalV "
			"GainX GainY GainPhi\n");

	skyfile = fopen("skyfile.dat", "w");
	fprintf(skyfile, "# Chan ObsSkyI ObsSkyQ ObsSkyU ObsSkyV CalSkyI CalSkyQ CalSkyU CalSkyV "
			"GainX GainY GainPhi\n");

	gainfile = fopen("gainfile.dat", "w");
	fprintf(gainfile, "# Chan ObsSkyI ObsSkyQ ObsSkyU ObsSkyV CalSkyI CalSkyQ CalSkyU CalSkyV "
			"GainX GainY GainPhi\n");
SSG*/
	//iterate over each time step
	for (i=0; i<size; i++) 
	{
		pRec = &(dataset[i]);
		if (pRec->flagBAD) continue;

		//printf("compute observed stokes ...\n");
		compute_observed_stokes(pRec, &ObsCal, &ObsSky);

		//printf("compute gains ...\n");
		compute_gains(pRec, &gain, lowchan, highchan, Tcalx, Tcaly);

		// calibrate the cal signal
		//printf("calibrate the cal ...\n");
		calibrate_stokes(&CalCal, &gain, &ObsCal, lowchan, highchan);

		// calibrate the sky signal
		//printf("calibrate the sky ...\n");
		calibrate_stokes(&TrueSky, &gain, &ObsSky, lowchan, highchan);

		//printf("compute final stokes ...\n");
		compute_final_stokes(pRec, &TrueSky, ignoreRFI, lowchan, highchan);

//		print_stokes(calfile, &ObsCal, &CalCal, &gain, lowchan, highchan);
//		print_stokes(skyfile, &ObsSky, &TrueSky, &gain, lowchan, highchan);
	}

//	fclose(calfile);
//	fclose(skyfile);
//	fclose(gainfile);
}


/*
 * Writes all the stokes data for a single channel to file.
 */
//void write_channel_data(SpecRecord dataset[], int size, int lowchan, int highchan) //SSG
void write_channel_data(SpecRecord dataset[], int size, int lowchan, int highchan, int beam)
{
	int n;
	double I, Q, U, V; //working copies
	double RA, DEC, AST;
	FILE * fluxfile;
	char filename[32+1];
	int chan;
	int found;//SSG
	char header[81];//SSG
	int readchan; //SSG
	FILE * beamgainfile; //SSG
	float gain[7];//SSG
	printf("DIAGNOSTIC:Before opening gainfile.\n");//SSG
	if(multibeam) //SSG
	{
		beamgainfile = fopen("../../beamgains.dat","r");//SSG
		if (beamgainfile == NULL) {//SSG
			printf("ERROR: unable to open gain file %s\n", filename);//SSG
			return;//SSG
		}//SSG

	else//SSG
		printf("DIAGNOSTIC:Opened gainfile.\n");//SSG
		fgets(header,80,beamgainfile);//SSG
		found = 0;//SSG
	}
	for (chan=lowchan; chan<highchan; chan++) 
	{
		snprintf(filename, 32, "fluxtime%03i.dat", chan);
		fluxfile = fopen(filename, "w");
		if (fluxfile == NULL) {
			printf("ERROR: unable to open file %s\n", filename);
			printf("DIAGNOSTIC: errno %d\n",errno);//SSG
			return;
		}
		fprintf(fluxfile, "# RA DEC AST I Q U V\n");
		//SSG
		if(multibeam)
		{
			while(!found)
			{
				fscanf(beamgainfile,"%d",&readchan);
//				printf("DIAGNOSTIC:%d %d\n",chan,readchan);
				if (readchan == chan)
				{
					fscanf(beamgainfile,"%f %f %f %f %f %f %f",&gain[0],&gain[1],&gain[2],&gain[3],&gain[4],&gain[5],&gain[6]);
		//			printf("Gain values: %f %f %f %f %f %f %f\n",gain[0],gain[1],gain[2],gain[3],gain[4],gain[5],gain[6]);
					found = 1;
				}
				else
				{
					fgets(header,80,beamgainfile);
				}
			}
			found = 0;
		}
		//SSG	
		for (n=0; n<size; n++)
		{
			SpecRecord * pRec = &(dataset[n]);
			RA = pRec->RA;
			DEC = pRec->DEC;
			AST = pRec->AST;

//			if (pRec->flagBAD || pRec->flagRFI[chan] != RFI_NONE) {
			if (pRec->flagBAD || pRec->flagRFI[chan] != RFI_NONE /*|| !isfinite(pRec->stokes.I[chan])*/) {//ssg
				I = Q = U = V = NAN;

			} else {
				//SSG
				if(multibeam)
				{
					I = pRec->stokes.I[chan]/gain[beam];
					Q = pRec->stokes.Q[chan];
					U = pRec->stokes.U[chan];
					V = pRec->stokes.V[chan];
				}
				//SSG
				else
				{
					I = pRec->stokes.I[chan];
					Q = pRec->stokes.Q[chan];
					U = pRec->stokes.U[chan];
					V = pRec->stokes.V[chan];
				}
			}

			fprintf(fluxfile,"%2.8f %2.8f %8.2f %4.6f %4.6f %4.6f %4.6f\n", 
					RA, DEC, AST, I, Q, U, V);

		} //end for timestep

		fclose(fluxfile);

	} //end for chan
	if(multibeam)//SSG
		fclose(beamgainfile);//SSG
}


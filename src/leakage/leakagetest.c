#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MAX_CHANNELS 256


  typedef int integer;
  typedef unsigned int logical;
  typedef float real;
  #define FORTRAN_LOGICAL_TRUE	1
  #define FORTRAN_LOGICAL_FALSE	0

//Desh's FORTRAN routine.
//extern void fetch_from_best_fit_model__ (integer *pixel_num, integer *pol_ch_num, real *frequency, char *file_of_model_sets, logical *use_only_latest_model, real *model_val, real *scaled_model_val, real *fractional_error, logical *model_was_used);
extern void fetch_from_best_fit_model__ (integer *pixel_num, integer *pol_ch_num, real *frequency, logical *use_only_latest_model, real *model_val, real *scaled_model_val, real *fractional_error, logical *model_was_used);
                                                                                                static void get_leakage_factors(int beam, float fcen, int lowchan, int highchan, float* L)
{
    int i, chan;
    float bw = 100.0/MAX_CHANNELS;
    int numPlanes;
    float F[MAX_CHANNELS];

    numPlanes = highchan-lowchan;

    //get the leakage factors from Desh's fortran routine
	for (i=0, chan=lowchan; chan<highchan; i++, chan++) {
		F[i] = fcen + bw*(chan-128);
    }
    for (i=0; i<numPlanes; i++)
    {
        real frequency;
        integer pol_chan = 1; //the same for both polarisations
        logical latest_model = FORTRAN_LOGICAL_TRUE;
        real model_val; //the result
        real ignored; //the scaled value does not apply
        real fractional_error;
        logical model_was_used; 

        frequency = F[i];
		fetch_from_best_fit_model__(&beam, &pol_chan, &frequency, &latest_model, 
            &model_val, &ignored, &fractional_error, &model_was_used);
        L[i] = model_val;

		if (model_was_used == FORTRAN_LOGICAL_FALSE) {
			printf("WARN: model was not used!\n");
		}

	}

}

                                                                                                
int main()
{
	float F[MAX_CHANNELS];

	float L0[MAX_CHANNELS];
	float L1[MAX_CHANNELS];
	float L2[MAX_CHANNELS];
	float L3[MAX_CHANNELS];
	float L4[MAX_CHANNELS];
	float L5[MAX_CHANNELS];
	float L6[MAX_CHANNELS];

	int i,chan;
	FILE *file;
	float fcen = 1440.0;
	int lowchan = 0;
	int highchan = MAX_CHANNELS;
	float bw = 100.0/MAX_CHANNELS;

	float U[MAX_CHANNELS];
	float V[MAX_CHANNELS];

	for (i=0, chan=lowchan; chan<highchan; i++, chan++) {
		F[i] = fcen + bw*(chan-128);
    }

	file = fopen("leakage.dat", "w");
	fprintf(file, "#Leakage factor from Desh's code for beam0\n");
	fprintf(file, "#chan frequency L zeta L*cos(zeta) L*sin(zeta)\n");
	get_leakage_factors(0, fcen, lowchan, highchan, L0);
	for (i=lowchan; i<highchan; i++) {
		float zeta = atan2(-V[i], U[i]);
		float L =  L0[i];
		fprintf(file, "%i %f %f %f %f %f\n", 
			i, F[i], zeta, L, L*cos(zeta), L*sin(zeta));
	}
	
}

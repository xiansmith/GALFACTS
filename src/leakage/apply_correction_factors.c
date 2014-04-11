#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/chardefs.h"
#include "programs/mathdefs.h"
#include "programs/misc.c"
#include "programs/misc_math.c"
#include "programs/coord_utils.c"
#include "programs/fitsio.c"


#define MAX_CHANNELS 256

//stuff to cope with fortran
  typedef int integer;
  typedef unsigned int logical;
  typedef float real;
  #define FORTRAN_LOGICAL_TRUE	1
  #define FORTRAN_LOGICAL_FALSE	0

//Desh's FORTRAN routine.
//extern void fetch_from_best_fit_model__ (integer *pixel_num, integer *pol_ch_num, real *frequency, char *file_of_model_sets, logical *use_only_latest_model, real *model_val, real *scaled_model_val, real *fractional_error, logical *model_was_used);
extern void fetch_from_best_fit_model__ (integer *pixel_num, integer *pol_ch_num, real *frequency, logical *use_only_latest_model, real *model_val, real *scaled_model_val, real *fractional_error, logical *model_was_used);
                                                                                                     


static void get_leakage_factors(int beam, float fcen, float* L)
{
    int i;
    float bw = 100.0/256.0;
    float F[MAX_CHANNELS];


    //get the leakage factors from Desh's fortran routine
	for (i=0; i<MAX_CHANNELS; i++) {
		F[i] = fcen + bw*(i-128);
    }
    for (i=0; i<MAX_CHANNELS; i++)
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


static void get_zeta_values(int beamid, float Z[])
{
	char filename[32+1];
	float zeta;
	int chan;
	FILE *zetafile;

	sprintf(filename, "beam%i_zeta.dat", beamid);
	zetafile = fopen(filename, "r");
	if (!zetafile) {
		printf("ERROR: unable to open '%s'\n", filename);
		return;
	}

	fscanf(zetafile, "%i %f\n", &chan, &zeta);
	while (!feof(zetafile)) {
		Z[chan] = zeta;
		fscanf(zetafile, "%i %f\n", &chan, &zeta);
	}
}

static void get_gamma_values(int beamid, float G[])
{
	char filename[32+1];
	float gamma;
	int chan;
	FILE *gammafile;

	sprintf(filename, "beam%i_gamma.dat", beamid);
	gammafile = fopen(filename, "r");
	if (!gammafile) {
		printf("ERROR: unable to open '%s'\n", filename);
		return;
	}

	fscanf(gammafile, "%i %f\n", &chan, &gamma);
	while (!feof(gammafile)) {
		G[chan] = gamma;
		fscanf(gammafile, "%i %f\n", &chan, &gamma);
	}

}


static void correct_cubes(int beamid, int lowchan, int highchan, float L[], float Z[], float G[])
{
    int i, j;
    char filename[32+1];
    float *I, *Q, *U, *V; //planes of data from the cubes
    int data_len; //number of points in a plane
    FILE *Icube_file, *Qcube_file, *Ucube_file, *Vcube_file;
    FILE *Icube_file_new, *Qcube_file_new, *Ucube_file_new, *Vcube_file_new;
    header_param_list Icube_hpar, Qcube_hpar, Ucube_hpar, Vcube_hpar;

    //open the existing fits cubes I, Q, U and V
    sprintf(filename, "beam%i_Icube.fits", beamid);
	Icube_file = fopen(filename, "r");
    sprintf(filename, "beam%i_Qcube.fits", beamid);
	Qcube_file = fopen(filename, "r");
    sprintf(filename, "beam%i_Ucube.fits", beamid);
	Ucube_file = fopen(filename, "r");
    sprintf(filename, "beam%i_Vcube.fits", beamid);
	Vcube_file = fopen(filename, "r");

    //sanity check
    if (!Icube_file || !Qcube_file || !Ucube_file || !Vcube_file) {
        printf("ERROR: unable to open input fits file\n");
        return;
    }

	//read the headers
	readfits_header(Icube_file, &Icube_hpar);
	readfits_header(Qcube_file, &Qcube_hpar);
	readfits_header(Ucube_file, &Ucube_hpar);
	readfits_header(Vcube_file, &Vcube_hpar);

	//allocate memory for a plane of each cube
	data_len = Icube_hpar.naxis[0] * Icube_hpar.naxis[1];
	I = (float*) calloc(data_len, sizeof (float));
	Q = (float*) calloc(data_len, sizeof (float));
	U = (float*) calloc(data_len, sizeof (float));
	V = (float*) calloc(data_len, sizeof (float));


	//open the new cubes
    sprintf(filename, "beam%i_Icube_corrected.fits", beamid);
	Icube_file_new = fopen(filename, "w");
    sprintf(filename, "beam%i_Qcube_corrected.fits", beamid);
	Qcube_file_new = fopen(filename, "w");
    sprintf(filename, "beam%i_Ucube_corrected.fits", beamid);
	Ucube_file_new = fopen(filename, "w");
    sprintf(filename, "beam%i_Vcube_corrected.fits", beamid);
	Vcube_file_new = fopen(filename, "w");

    //sanity check the new files
    if (!Icube_file_new || !Qcube_file_new || !Ucube_file_new || !Vcube_file_new) {
        printf("ERROR: unable to open output fits file\n");
        return;
    }

    //write the headers for the new corrected cubes
	writefits_header(Icube_file_new, &Icube_hpar);
	writefits_header(Qcube_file_new, &Qcube_hpar);
	writefits_header(Ucube_file_new, &Ucube_hpar);
	writefits_header(Vcube_file_new, &Vcube_hpar);

	//iterate over the cubes and write out the corrected data
	for (i=0; i<Icube_hpar.naxis[2]; i++) 
	{
        float zeta;
		float gamma, gamma_sq_inv;
		int chan = i+lowchan;
		printf("channel %i\n", chan);

		zeta = Z[chan];
		gamma = G[chan];
		gamma_sq_inv = 1.0/(gamma*gamma);

		//read a plane
		readfits_plane(Icube_file, I, &Icube_hpar);
		readfits_plane(Qcube_file, Q, &Qcube_hpar);
		readfits_plane(Ucube_file, U, &Ucube_hpar);
		readfits_plane(Vcube_file, V, &Vcube_hpar);

        //iterate over the points in the plane
        for (j=0; j<data_len; j++) 
        {
			//cache the old values so we work with a consistent set
			float Iold = I[j];
			float Qold = Q[j];
			float Uold = U[j];
			float Vold = V[j];

			//do the leakage correction (involving L and zeta) 
			//and the gamma correction
			I[j] = (Iold/2)*(1+gamma_sq_inv) + (Qold/2)*(1-gamma_sq_inv);
			Q[j] = (Iold/2)*(1-gamma_sq_inv) + (Qold/2)*(1+gamma_sq_inv);
            U[j] = (Uold - Iold*L[chan]*cos(zeta))/gamma;
            V[j] = (Vold + Iold*L[chan]*sin(zeta))/gamma;
        }

        //write the planes
        writefits_plane(Icube_file_new, I, &Icube_hpar);
        writefits_plane(Qcube_file_new, Q, &Qcube_hpar);
        writefits_plane(Ucube_file_new, U, &Ucube_hpar);
        writefits_plane(Vcube_file_new, V, &Vcube_hpar);
	}

    //pad and close the new cubes
    writefits_pad_end(Icube_file_new, &Icube_hpar);
    writefits_pad_end(Qcube_file_new, &Qcube_hpar);
    writefits_pad_end(Ucube_file_new, &Ucube_hpar);
    writefits_pad_end(Vcube_file_new, &Vcube_hpar);

    //cleanup
    close(Icube_file);
    close(Qcube_file);
    close(Ucube_file);
    close(Vcube_file);
    close(Icube_file_new);
    close(Qcube_file_new);
    close(Ucube_file_new);
    close(Vcube_file_new);
	free(I);
	free(Q);
	free(U);
	free(V);

}                                                           

int main(int argc, char *argv[])
{
    int beamid;
    int lowchan;
    int highchan;
    float L[MAX_CHANNELS];
    float Z[MAX_CHANNELS];
    float G[MAX_CHANNELS];


    if (argc != 4) {
        printf("Usage: %s <beamid> <lowchan> <highchan>\n", argv[0]);
        printf("ex: %s 0 30 230\n", argv[0]);
        return EXIT_FAILURE;
    }
    else
    {
        beamid = atoi(argv[1]);
        lowchan = atoi(argv[2]);
        highchan = atoi(argv[3]);
    }

    printf("get leakage factors ...\n");
    get_leakage_factors(beamid, 1440, L);

	printf("get zeta values from file ...\n");
	get_zeta_values(beamid, Z);

	printf("get gamma values from file ...\n");
	get_gamma_values(beamid, G);

    printf("correct cubes for leakage ...\n");
    correct_cubes(beamid, lowchan, highchan, L, Z, G);

	printf("done!\n");

    return EXIT_SUCCESS;

}



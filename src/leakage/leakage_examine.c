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

static void read_zeta(float Z)
{
	zeta = atan2(-V[j], U[j]);
}

static void get_leakage_factors(int beam, float fcen, int lowchan, int highchan, float* L)
{
    int i, chan;
    float bw = 100.0/256.0;
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


static void correct_cubes_for_leakage(int beamid, float L[])
{
    int i, j;
    char filename[32+1];
    float *I, *U, *V; //planes of data from the cubes
    int data_len; //number of points in a plane
    FILE *Icube_file, *Ucube_file, *Vcube_file;
    FILE *Ucube_file_new, *Vcube_file_new;
    header_param_list Icube_hpar, Ucube_hpar, Vcube_hpar;

    //open the fits cubes I, U and V
    sprintf(filename, "beam%i_Icube.fits", beamid);
	Icube_file = fopen(filename, "r");
	readfits_header(Icube_file, &Icube_hpar);

    sprintf(filename, "beam%i_Ucube.fits", beamid);
	Ucube_file = fopen(filename, "r");
	readfits_header(Ucube_file, &Ucube_hpar);

    sprintf(filename, "beam%i_Vcube.fits", beamid);
	Vcube_file = fopen(filename, "r");
	readfits_header(Vcube_file, &Vcube_hpar);

    //sanity check
    if (!Icube_file || !Ucube_file || !Vcube_file) {
        printf("ERROR: unable to open input fits file\n");
        return;
    }

	//allocate memory for a plane of each cube
	data_len = Icube_hpar.naxis[0] * Icube_hpar.naxis[1];
	I = (float*) calloc(data_len, sizeof (float));
	U = (float*) calloc(data_len, sizeof (float));
	V = (float*) calloc(data_len, sizeof (float));


    //write the headers for the new leakage corrected cubes
    sprintf(filename, "beam%i_Ucube_Lzeta.fits", beamid);
	Ucube_file_new = fopen(filename, "w");
    //strcat(Ucube_hpar.object, " L");
	writefits_header(Ucube_file_new, &Ucube_hpar);

    sprintf(filename, "beam%i_Vcube_Lzeta.fits", beamid);
	Vcube_file_new = fopen(filename, "w");
    //strcat(Ucube_hpar.object, " L");
	writefits_header(Vcube_file_new, &Vcube_hpar);

    //sanity check
    if (!Ucube_file_new || !Vcube_file_new) {
        printf("ERROR: unable to open output fits file\n");
        return;
    }

	//iterate over the cube and write out the corrected data
	for (i=0; i<Icube_hpar.naxis[2]; i++) 
	{
        float zeta;

		//read a plane
        printf("read plane %i\n", i);
		readfits_plane(Icube_file, I, &Icube_hpar);
		readfits_plane(Ucube_file, U, &Ucube_hpar);
		readfits_plane(Vcube_file, V, &Vcube_hpar);

        //iterate over the points in the plane
        for (j=0; j<data_len; j++) 
        {
            //select the leakage correction factor
            //printf("leakage correction angle %f\n", zeta);
            zeta = atan2(-V[j],U[j]);
            U[j] = L[i]*cos(zeta);
            V[j] = L[i]*sin(zeta);
        }

        //write the planes
        writefits_plane(Ucube_file_new, U, &Ucube_hpar);
        writefits_plane(Vcube_file_new, V, &Vcube_hpar);
	}

    //pad and close the new cubes
    writefits_pad_end(Ucube_file_new, &Ucube_hpar);
    writefits_pad_end(Vcube_file_new, &Vcube_hpar);

    //cleanup
    close(Icube_file);
    close(Ucube_file);
    close(Vcube_file);
    close(Ucube_file_new);
    close(Vcube_file_new);
	free(I);
	free(U);
	free(V);

}                                                           

int main(int argc, char *argv[])
{
    int beamid;
    int lowchan;
    int highchan;
    float L[MAX_CHANNELS];


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
    get_leakage_factors(beamid, 1440, lowchan, highchan, L);

    printf("correct cubes for leakage ...\n");
    correct_cubes_for_leakage(beamid, L);

    return EXIT_SUCCESS;

}



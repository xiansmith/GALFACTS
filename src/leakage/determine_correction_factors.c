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
#include "float.h"

#define MAX_CHANNELS 256

static int find_max_amplitude(FILE * cube_file)
{
	float *plane_data, *avg_data, *count_data;
	header_param_list cube_hpar;
	int datalen;
	float maxval;
	int maxpos;
	int data_len;
	int i, j;

	readfits_header(cube_file, &cube_hpar);

	data_len = cube_hpar.naxis[0] * cube_hpar.naxis[1];

	plane_data = (float*) calloc(data_len, sizeof (float));
	avg_data = (float*) calloc(data_len, sizeof (float));
	count_data = (float*) calloc(data_len, sizeof (float));

	//iterate over the cube
	//calculate average
	for (i=0; i<cube_hpar.naxis[2]; i++) {
		//read a plane
		readfits_plane(cube_file, plane_data, &cube_hpar);
		for (j=0; j<data_len; j++) {
			//add to average and increment count if its not blank
			if (plane_data[j] != BLANK_PIXEL) {
				avg_data[j] += plane_data[j];
				count_data[j]++;
			}
		}
	}
	free(plane_data);
	for (j=0; j<data_len; j++) {
		if (count_data[j] == 0) {
			avg_data[j] = FLT_MIN;
		} else {
			avg_data[j] /= count_data[j];
		}
	}
	free(count_data);

	//find max amplitude positions
	maxval = FLT_MIN;
	for (j=0; j<data_len; j++) {
			if (avg_data[j] > maxval) {
				maxval = avg_data[j];
				maxpos = j;
			}
	}
	free(avg_data);
	
	//return results
	printf("max amplitude: %f\n", maxval);
	return maxpos;
}

static void determine_gamma_correction(int beamid, int pos, float G[])
{
    int i, j;
    char filename[32+1];
    float *I, *Q; //planes of data from the cubes
    int data_len; //number of points in a plane
    FILE *Icube_file, *Qcube_file;
    header_param_list Icube_hpar, Qcube_hpar;

	//open the U and V cubes
    sprintf(filename, "beam%i_Icube.fits", beamid);
	Icube_file = fopen(filename, "r");
	readfits_header(Icube_file, &Icube_hpar);

    sprintf(filename, "beam%i_Qcube.fits", beamid);
	Qcube_file = fopen(filename, "r");
	readfits_header(Qcube_file, &Qcube_hpar);

    //sanity check
    if (!Icube_file || !Qcube_file) {
        printf("ERROR: unable to open input fits file\n");
        return;
    }

	//allocate memory for a plane of each cube
	data_len = Icube_hpar.naxis[0] * Icube_hpar.naxis[1];
	I = (float*) calloc(data_len, sizeof (float));
	Q = (float*) calloc(data_len, sizeof (float));

	//iterate over the cube and write out zeta
	for (i=0; i<Icube_hpar.naxis[2]; i++) 
	{
		float Ix, Iy;

		//read a plane
		readfits_plane(Icube_file, I, &Icube_hpar);
		readfits_plane(Qcube_file, Q, &Qcube_hpar);
		
		//calculate and store gamma
		Ix = (I[i]+Q[i])/2.0;
		Iy = (I[i]-Q[i])/2.0;
		G[i] = sqrt(Iy/Ix);
	}

    //cleanup
    close(Icube_file);
    close(Qcube_file);
	free(I);
	free(Q);
}

static void determine_leakage_zeta(int beamid, int pos, float Z[])
{
    int i, j;
    char filename[32+1];
    float *U, *V; //planes of data from the cubes
    int data_len; //number of points in a plane
    FILE *Ucube_file, *Vcube_file;
    header_param_list Ucube_hpar, Vcube_hpar;

	//open the U and V cubes
    sprintf(filename, "beam%i_Ucube.fits", beamid);
	Ucube_file = fopen(filename, "r");
	readfits_header(Ucube_file, &Ucube_hpar);

    sprintf(filename, "beam%i_Vcube.fits", beamid);
	Vcube_file = fopen(filename, "r");
	readfits_header(Vcube_file, &Vcube_hpar);

    //sanity check
    if (!Ucube_file || !Vcube_file) {
        printf("ERROR: unable to open input fits file\n");
        return;
    }

	//allocate memory for a plane of each cube
	data_len = Ucube_hpar.naxis[0] * Ucube_hpar.naxis[1];
	U = (float*) calloc(data_len, sizeof (float));
	V = (float*) calloc(data_len, sizeof (float));

	//iterate over the cube and write out zeta
	for (i=0; i<Ucube_hpar.naxis[2]; i++) 
	{
        float zeta;

		//read a plane
		readfits_plane(Ucube_file, U, &Ucube_hpar);
		readfits_plane(Vcube_file, V, &Vcube_hpar);

		//printf("leakage correction angle %f\n", zeta);
		Z[i] = atan2(-V[pos], U[pos]);
	}

    //cleanup
    close(Ucube_file);
    close(Vcube_file);
	free(U);
	free(V);
}                                                           

int main(int argc, char *argv[])
{
	int i;
    int beamid;
    int lowchan, highchan;
    float Z[MAX_CHANNELS];
    float G[MAX_CHANNELS];
	FILE * zetafile;
	FILE * gammafile;
	char filename[32+1];
	FILE *Icube_file;
	int pos;

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

	//open the output file
	sprintf(filename, "beam%i_zeta.dat", beamid);
	zetafile = fopen(filename, "w");
	sprintf(filename, "beam%i_gamma.dat", beamid);
	gammafile = fopen(filename, "w");
	if (!zetafile || !gammafile) {
		printf("ERROR: unable to open output file\n");
		return EXIT_FAILURE;
	}

	printf("find the position of max I amplitude ...\n");
    sprintf(filename, "beam%i_Icube.fits", beamid);
	Icube_file = fopen(filename, "r");
	pos = find_max_amplitude(Icube_file);
    close(Icube_file);


	printf("determine leakage zeta ...\n");
    determine_leakage_zeta(beamid, pos, Z);

	printf("determine gamma correction ...\n");
    determine_gamma_correction(beamid, pos, G);

	//write the results
	printf("write the data files ...\n");
	for (i=0; i<lowchan; i++) {
		fprintf(zetafile, "%i %f\n", i, 0.0);
		fprintf(gammafile, "%i %f\n", i, 0.0);
	}
	for (i=lowchan; i<highchan; i++) {
		fprintf(zetafile, "%i %f\n", i, Z[i-lowchan]);
		fprintf(gammafile, "%i %f\n", i, G[i-lowchan]);
	}
	for (i=highchan; i<MAX_CHANNELS; i++) {
		fprintf(zetafile, "%i %f\n", i, 0.0);
		fprintf(gammafile, "%i %f\n", i, 0.0);
	}

	close(zetafile);
	close(gammafile);
	printf("done!\n");

    return EXIT_SUCCESS;
}



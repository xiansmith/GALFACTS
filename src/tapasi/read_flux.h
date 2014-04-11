/*
These data structures deal with the Stokes flux for a particular
channel.  

FluxRecord - a single line in a fluxtime.dat file
FluxData - a collection of flux records for a particular day

*/

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double I;
    double Q;
    double U;
    double V;
} StokesParams;

typedef struct {
	float RA;
	float DEC;
	float AST;
	StokesParams stokes;
} FluxRecord;

typedef struct {
	char mjd[8+1];
	int numRecords;
	float RAmin;
	float RAmax;
	FluxRecord * records;
} FluxData;

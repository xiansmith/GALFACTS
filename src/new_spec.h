/*
File: spec.h
The data file format produced from an Arecibo observation consist of two key files:
.cfg and .spec.  The .cfg contains the meta-data about the observation and is plain text.
The .spec is a binary file that contains the spectral and pointing information.

Reference is made to the FORTRAN definition of the code that procuces these data files:
http://www2.naic.edu/~tghosh/a1947/continuum_pointing.inc

*/

#ifndef _NEW_SPEC_H
#define _NEW_SPEC_H

#include "common.h"
#define MAX_CHANNELS1 4096

/* 
Object to contain the data from the .cfg file.
Describes the metadata from the observation.
*/
typedef struct {
	float integrationTime; //miliseconds
	int stokesSetSize; //#bytes per stokes set
	float centerMHz; //MHz
	float bandwitdhkHz; //kHz
	int startChanNum; //start channel number
	int samplesPerStokesSet; //no. of samples in a Stokes set
	int stokesProducts; //no. of Stokes products
	int numStokesSets; //and no of such sets corresponding a given integration (i.e. cal_on,cal_off)
	int crap; //undocumented
	char observingTag[32]; //Source-name or observing-mode tag, WCAL tag to indicate winking cal data 
	int day; //date dd mm yyyy
	char observatoryCode[16]; //Observatory code
} ConfigData;


/* 
CentralBeamBlock
Contains the pointing information for the primary beam, (beam0), along with 
some common info that is general to the position of the receiver on the sky,
as opposed to a particular beam.
*/
typedef struct {
    double raj_true_in_degrees;
    double decj_true_in_degrees;
    double arecibo_local_mean_sidereal_time_in_sec;
    double mjd_last_five_digits;
    double az_cur_in_degrees;
    double za_cur_in_degrees;
    double raj_requested_in_degrees;
    double decj_requested_in_degrees;
    double alfa_rotation_angle_in_degrees;
} SpecPointingBlock;

/*
OuterBeamBlock
Contains the pointing information for an outer beam.   Some information from
the CentralBeamBlock will be required to construct the precise position
of any particular beam.
*/

/*
The entire SpecPointingBlock consists of the central beam block and
7 outer beams.  outerBeams[0] has data from the central beam to allow
iterating directly over the beams.
The size of this structure must be 512 bytes.
*/

#define SIZEOF_SPECPOINTINGBLOCK 72

//TODO: to port to different architectures, ensure that
//sizeof(SpecPointingBlock) == SIZEOF_SPECPOINTINGBLOCK

/*
PolSet
A set of the auto and cross correlated polarisation values.
*/
typedef struct {
	unsigned long int A;
	unsigned long int B;
	long int U;
	long int V;
	long int fft_weight;
} PolAvg;

/*
PolSet
A set of the auto and cross correlated polarisation values.
*/
typedef struct {
	unsigned int A[MAX_CHANNELS1];
	unsigned int B[MAX_CHANNELS1];
	int U[MAX_CHANNELS1];
	int V[MAX_CHANNELS1];
	int fft_weight;
} PolSet;

/*
StokesSet
A set of the four stokes parameters.
*/
typedef struct {
	float I[MAX_CHANNELS1];
	float Q[MAX_CHANNELS1];
	float U[MAX_CHANNELS1];
	float V[MAX_CHANNELS1];
} StokesSet;

/*
SpecRecord
This is the core data structure for operating on data in the time domain: ie
calibration.  RA,DEC,AST is the pointing information.  calon and caloff is the
raw spectral data, cal is the computed value of the cal signal.  flagRFI
contain the RFI bit flags on a per channel basis.  flagBAD is non-zero when
the entire spectra is marked bad.  stokes contains the spectra of the computed,
calibrated stokes parameters.
*/
typedef struct {
	double RA;
	double DEC;
	double AST;
	PolSet calon;
	PolSet caloff;
	PolSet cal;
	enum RFIFlags flagRFI[MAX_CHANNELS1];
	char  flagBAD;
	StokesSet stokes;
} SpecRecord;


int read_datafile(FILE * pFile, SpecRecord ** pDataset);
int read_cfgfile(FILE * pFile, ConfigData * pCfg);

void print_data(SpecRecord dataset[], int numRecords, int lowchan, int highchan, float freq[]);


#endif //_NEW_SPEC_H


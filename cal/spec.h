/*
File: spec.h
The data file format produced from an Arecibo observation consist of two key files:
.cfg and .spec.  The .cfg contains the meta-data about the observation and is plain text.
The .spec is a binary file that contains the spectral and pointing information.

Reference is made to the FORTRAN definition of the code that procuces these data files:
http://www2.naic.edu/~tghosh/a1947/continuum_pointing.inc

*/

#ifndef _SPEC_H
#define _SPEC_H

#include "common.h"


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
	int month;
	int year;
	int hour; // AST start time hh mm ss.ssssss
	int minute;
	float second;
	char observatoryCode[16]; //Observatory code
	int bandFlip; //flag to indicate band-flip (-1 does indicate band inversion, +1: no flip)
} ConfigData;


/* 
CentralBeamBlock
Contains the pointing information for the primary beam, (beam0), along with 
some common info that is general to the position of the receiver on the sky,
as opposed to a particular beam.
*/
typedef struct {
    double raj_true_in_hours;
    double decj_true_in_degrees;
    char epoch[8];
    double atlantic_solar_time_now_in_sec;
    double arecibo_local_mean_sidereal_time_in_sec;
    double mjd_last_five_digits;
    double az_cur_in_degrees;
    double za_cur_in_degrees;
    double az_encoder_model_offset_in_degrees;
    double za_encoder_model_offset_in_degrees;
    double az_requested_in_degrees;
    double za_requested_in_degrees;
    double heliocentric_velocity_in_km_per_second;
    double geocentric_velocity_in_km_per_second;
    double parallactic_angle_now_in_dregrees;
    double alfa_rotation_angle_in_degrees;
} CentralBeamBlock;

/*
OuterBeamBlock
Contains the pointing information for an outer beam.   Some information from
the CentralBeamBlock will be required to construct the precise position
of any particular beam.
*/
typedef struct {
    double raj_true_in_hours;
    double decj_true_in_degrees;
    double az_cur_in_degrees;
    double za_cur_in_degrees;
    double az_alfa_offset_now_in_arcsec;
    double za_alfa_offset_now_in_arcsec;
    double az_offset_for_unrotated_alfa_in_arcsec;
    double za_offset_for_unrotated_alfa_in_arcsec;
} OuterBeamBlock;

/*
The entire SpecPointingBlock consists of the central beam block and
7 outer beams.  outerBeams[0] has data from the central beam to allow
iterating directly over the beams.
The size of this structure must be 512 bytes.
*/
typedef struct {
    CentralBeamBlock centralBeam; //central beam
    OuterBeamBlock outerBeams[6]; //beam at postion i
} SpecPointingBlock;

#define SIZEOF_SPECPOINTINGBLOCK 512

//TODO: to port to different architectures, ensure that
//sizeof(SpecPointingBlock) == SIZEOF_SPECPOINTINGBLOCK

/*
PolSet
A set of the auto and cross correlated polarisation values.
*/
typedef struct {
	double xx;
	double yy;
	double xy;
	double yx;
} PolAvg;

/*
PolSet
A set of the auto and cross correlated polarisation values.
*/
typedef struct {
	float xx[MAX_CHANNELS];
	float yy[MAX_CHANNELS];
	float xy[MAX_CHANNELS];
	float yx[MAX_CHANNELS];
} PolSet;

/*
StokesSet
A set of the four stokes parameters.
*/
typedef struct {
	float I[MAX_CHANNELS];
	float Q[MAX_CHANNELS];
	float U[MAX_CHANNELS];
	float V[MAX_CHANNELS];
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
	//double freq;  add this once fluxdata loads and computes frequency
	PolSet calon;
	PolSet caloff;
	PolSet cal;
	enum RFIFlags flagRFI[MAX_CHANNELS];
	int  flagBAD;
	StokesSet stokes;
} SpecRecord;


int read_datafile(FILE * pFile, SpecRecord ** pDataset, int beam,float RAmin,float RAmax,float DECmin,float DECmax);
int read_cfgfile(FILE * pFile, ConfigData * pCfg);

void print_data(SpecRecord dataset[], int numRecords, int lowchan, int highchan, float freq[]);


#endif //_SPEC_H


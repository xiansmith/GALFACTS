#include "common.h"
#include "spec.h"

#ifndef _STOKES_H
#define _STOKES_H

void calculate_stokes(SpecRecord dataset[], int size, int lowchan, int highchan, int RFIF, int calskyfiles, float Tcalx[], float Tcaly[], int uvDenoising, float uvDenoisingTau, float uvDenoisingLambda, int start, int end);
void write_channel_data(SpecRecord dataset[], int size, int lowchan, int highchan);
void write_binary_channel_data(SpecRecord dataset[], int size, int lowchan, int highchan);
void average_stokes(SpecRecord dataset[], int size, int lowchan, int highchan, float hidrogenfreq, float hidrogenband, float freq[]);
void correct_beamgains(SpecRecord dataset[], int size, int lowchan, int highchan, int beam);
void write_binary_channel_data_single_file(SpecRecord dataset[], int numRecords, int lowchan, int highchan);
void write_binary_average_data(SpecRecord dataset[], int size, int lowchan, int highchan);
void write_rfi_data( SpecRecord dataset[], int numRecords, int lowchan, int highchan);
#endif //_STOKES_H


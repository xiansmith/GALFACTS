#include "common.h"
#include "spec.h"

#ifndef _STOKES_H
#define _STOKES_H

void calculate_stokes(SpecRecord dataset[], const int size, const int lowchan,
	const int highchan, const int RFIF, const int calskyfiles,
	const float Tcalx[], const float Tcaly[], const int uvDenoising,
	const float uvDenoisingTau, const float uvDenoisingLambda, const int start,
	const int end);

void write_channel_data(const SpecRecord dataset[], const int size, const int lowchan,
	const int highchan);

void write_binary_channel_data(const SpecRecord dataset[], const int size,
	const int lowchan, const int highchan);

void average_stokes(SpecRecord dataset[], const int size, const int lowchan,
	const int highchan, const float hidrogenfreq, const float hidrogenband,
	const float freq[]);

void correct_beamgains(SpecRecord dataset[], const int size, const int lowchan,
	const int highchan, const int beam);

void write_binary_channel_data_single_file(const SpecRecord dataset[],
	const int numRecords, const int lowchan, const int highchan);

void write_binary_average_data(const SpecRecord dataset[], const int size,
	const int lowchan, const int highchan);

void write_rfi_data(const SpecRecord dataset[], const int numRecords,
	const int lowchan, const int highchan);

#endif //_STOKES_H

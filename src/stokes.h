#include "common.h"
#include "spec.h"

#ifndef _STOKES_H
#define _STOKES_H

void calculate_stokes(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI, float Tcalx[], float Tcaly[]);
void write_channel_data(SpecRecord dataset[], int size, int lowchan, int highchan);
void average_stokes(SpecRecord dataset[], int size, int lowchan, int highchan);
void correct_beamgains(SpecRecord dataset[], int size, int lowchan, int highchan, int beam);

#endif //_STOKES_H


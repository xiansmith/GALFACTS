#ifndef _BEAMGAINS_H
#define _BEAMGAINS_H
#include "common.h"
#define NUM_BEAMS 7

void get_peak_power_response(char * filename,float peak_response[MAX_CHANNELS],int startchannel,int endchannel);

#endif


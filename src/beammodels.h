#ifndef _BEAMMODELS_H
#define _BEAMMODELS_H
//#include <stdio.h>
#include "common.h"
//#define NUM_BEAMS 7

void get_avg(char * filename,float avg[MAX_CHANNELS],int startchannel,int endchannel);
//void make_beam_model(FILE * beammodelfile,int channel,int beamno,float RA,float DEC,float radius,float max_response,float min_response,float avgQ,\
float avgU,float avgV);
void make_beam_model(int channel,int beamno,float RA,float DEC,float radius);
void get_peak_power_coord(char * filename,float RA[MAX_CHANNELS],float DEC[MAX_CHANNELS],float max_response[MAX_CHANNELS],float min_response[MAX_CHANNELS],int startchannel,int endchannel);
//void get_maxmin_power_response(char * filename,float max_response[MAX_CHANNELS],float min_response[MAX_CHANNELS],int startchannel,int endchannel);

#endif


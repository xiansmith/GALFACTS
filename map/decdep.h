#include "common.h"
#include "fluxdata.h"

#ifndef _DECDEPENDENCE_H
#define _DECDEPENDENCE_H

void calculate_dec_dependence(FluxWappData * wappdata, int order, int chan, float *cIc, float *cQc, float *cUc, float *cVc, int avg);
void beam_gain_calibration(FluxWappData * wappdata);
void beam_gain_calibration_table(FluxWappData * wappdata, int cal_low, int cal_high, float cal_table[][7], int chan);

#endif


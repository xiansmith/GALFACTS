#include "common.h"
#include "spec.h"

#ifndef _CALIBRATE_H
#define _CALIBRATE_H

void compute_raw_cal(SpecRecord dataset[], int size, int lowchan, int highchan);
void linear_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int RFIF);
void simple_smooth_cal(SpecRecord dataset[], int start, int end, int records, int lowchan, int highchan, int twindow, int fwindow);
void rolling_smooth_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int window, int cwindow);
void write_band(float* Fxx,float* Fyy, float *Fxy, float *Fyx,int r,int lowchan, int highchan);
#endif //_CALIBRATE_H


#include "common.h"
#include "spec.h"

#ifndef _CALIBRATE_H
#define _CALIBRATE_H

void compute_raw_cal(SpecRecord dataset[], int size, int lowchan, int highchan);
void linear_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int RFIF);
void smooth_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int window);

#endif //_CALIBRATE_H


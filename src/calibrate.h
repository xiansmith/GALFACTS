#include "common.h"
#include "spec.h"

#ifndef _CALIBRATE_H
#define _CALIBRATE_H

void compute_raw_cal(SpecRecord dataset[], int size);
void linear_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI);
void write_cal_fits(SpecRecord dataset[], int size, float fcen, float df);
void smooth_cal_bandaverage(SpecRecord dataset[], int size, int lowchan, int highchan, int smooth_width, float nsigma);

#endif //_CALIBRATE_H


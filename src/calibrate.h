#include "common.h"
#include "spec.h"

#ifndef _CALIBRATE_H
#define _CALIBRATE_H

void compute_raw_cal(SpecRecord dataset[],int lowchan,int highchan, int size);
//ssg void compute_raw_cal(SpecRecord dataset[], int size);
void print_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI);
void linear_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI);
void write_cal_fits(SpecRecord dataset[], int size, float fcen, float df);

#endif //_CALIBRATE_H


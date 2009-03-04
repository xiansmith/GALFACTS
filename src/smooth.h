#ifndef _SMOOTH_H
#define _SMOOTH_H

#include "common.h"
#include "spec.h"

void perform_freq_smoothing(SpecRecord dataset[], int numRecords, int lowchan, int highchan);

#endif //_SMOOTH_H


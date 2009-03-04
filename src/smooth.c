#include "smooth.h"

static inline float hanning_3(float * data, int i)
{
	return data[i-1]*0.25 + data[i]*0.50 + data[i+1]*0.25;
}

void perform_freq_smoothing(SpecRecord dataset[], int numRecords, int lowchan, int highchan)
{
	int n, chan;
	PolSet smooth_on;
	PolSet smooth_off;

	for (n=0; n<numRecords; n++)
	{
		SpecRecord * pRec = &(dataset[n]);	

//		for (chan=1; chan<MAX_CHANNELS-1; chan++) {
		for (chan=lowchan+1; chan<highchan-1; chan++) {
			smooth_on.xx[chan] = hanning_3(pRec->calon.xx, chan);
			smooth_on.xy[chan] = hanning_3(pRec->calon.xy, chan);
			smooth_on.yx[chan] = hanning_3(pRec->calon.yx, chan);
			smooth_on.yy[chan] = hanning_3(pRec->calon.yy, chan);
			smooth_off.xx[chan] = hanning_3(pRec->caloff.xx, chan);
			smooth_off.xy[chan] = hanning_3(pRec->caloff.xy, chan);
			smooth_off.yx[chan] = hanning_3(pRec->caloff.yx, chan);
			smooth_off.yy[chan] = hanning_3(pRec->caloff.yy, chan);
		}
//		for (chan=1; chan<MAX_CHANNELS-1; chan++) {
		for (chan=lowchan+1; chan<highchan-1; chan++) {
			pRec->calon.xx[chan] = smooth_on.xx[chan];
			pRec->calon.xy[chan] = smooth_on.xy[chan];
			pRec->calon.yx[chan] = smooth_on.yx[chan];
			pRec->calon.yy[chan] = smooth_on.yy[chan];
			pRec->caloff.xx[chan] = smooth_off.xx[chan];
			pRec->caloff.xy[chan] = smooth_off.xy[chan];
			pRec->caloff.yx[chan] = smooth_off.yx[chan];
			pRec->caloff.yy[chan] = smooth_off.yy[chan];
		}
	}	
}


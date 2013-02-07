#ifndef _RFI_H
#define _RFI_H

#include "common.h"
#include "spec.h"

void mark_bad_channels(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[], int *badchannels);

void rfi_detection_frequency_domain(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[]);
void rfi_detection_time_domain1(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[]);
void rfi_detection_time_domain2(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[]);
void rfi_detection_time_domain3(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[]);
void outofbandrfi_ann(SpecRecord dataset[], int size, int lowchan, int highchan);
void rfi_ann(SpecRecord dataset[], int size, int lowchan, int highchan, float freq[]); 
void strongsource_rfi_exclusions( SpecRecord dataset[], int numRecords, int lowchan, int highchan );
#endif


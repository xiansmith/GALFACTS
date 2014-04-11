#ifndef _COMMON_H
#define _COMMON_H

#include <stdio.h>

//The maximum number of channels in the dataset
//#define MAX_CHANNELS 256
//#define MAX_CHANNELS 2048
#define MAX_CHANNELS 4096

//The maximum number of datapoints in a single scan
#define MAX_SCAN_LEN 4000

//The maximum number of days
#define MAX_NUM_DAYS 250

//The maximum number of scans in one direction in a single day
#define MAX_NUM_SCANS 40

#define NUM_BEAMS 7 //SSG
//enum STOKES {STOKES_I, STOKES_Q, STOKES_U, STOKES_V};
#define MULTIBEAM 8 //ssg

enum RFIFlags {	RFI_NONE =		0x0000,
				RFI_CALOFF_XX =	0x0001,
				RFI_CALOFF_YY =	0x0002,
				RFI_CALON_XX =	0x0004,
				RFI_CALON_YY =	0x0008,
				RFI_CALOFF_XY =	0x0010,
				RFI_CALOFF_YX =	0x0020,
				RFI_CALON_XY =	0x0040,
				RFI_CALON_YX =	0x0080,
				RFI_OVERLIMIT = 0x0100,
				RFI_SIGMA_EXCEEDED = 0x0200,
				RFI_SPAN = 0x1000,
				RFI_OUTOFBAND = 0x2000,
};


int get_date_dirs(const char * dir, char ** datedirs[]);
int get_flux_files(const char * dir, char ** fluxfiles[]);


#endif //_COMMON_H

#ifndef GALFACTSLIBH
#define GALFACTSLIBH 

#include <fitsio.h>

// MAX length filename (include directory)
#define	GFL_MAXFNAME 200
#define	GFL_FTYPE_STR_DECTM   "DECTIME"
#define	GFL_FTYPE_STR_DECFRQ  "DECFREQ"
#define GFL_NUM_SHORTS_1STAT  10

/*
 * fits pdev status structure 
 * Note: the fits files contain 10 shorts. the fill bytes are not written.
*/
typedef struct {
    unsigned short seqNum; // counts modulo 64K
    unsigned short fftAccum; // we actually accumulated (includes blanking)
    unsigned short cal     ; // non zero cal on sometime during accuumulation
    unsigned short ovrflAdc; //uses 4 bits. num errors about 2^(n-1)  upto 2^12.
    unsigned short ovrflPfb; // poly phase filter, bufferfly. these wrap around
    unsigned short ovrflVShift; //upshift voltage before compute power
    unsigned short ovrflAccS2S3; //accum S2,S3
    unsigned short ovrflAccS0S1;// accum S0,S1
    unsigned short ovrflAshS2S3;// upshift after accum S2,S3
    unsigned short ovrflAshS0S1;// upshift after accum S0,S1
} GFLIB_PDEVSTAT;
// ----------------------------------------------------------------------
//  dec Time info
//
//   column numbers in fits file
//
typedef struct {
	int	dataCon;		// cal on
	int	dataCoff;		// cal off
	int	statCon;		// stat cal on
	int	statCoff;		// stat cal off
} GFLIB_DTM_COLNUMS;

typedef struct {
	GFLIB_DTM_COLNUMS colNums; // column numbers in fits file
    int     calCyclesSum; // cal cycles added time decimation
    int     sumSeq[2];    // nSum,Nskip each calOn,calOff,[19,1,19,1]
	int		*pdataOn;		  // malloc pon data buffer here
	int		*pdataOff;		 // malloc poff data buffer here
	GFLIB_PDEVSTAT *pstatOn;  // malloc pon  stat buffer here
	GFLIB_PDEVSTAT *pstatOff; // malloc poff stat buffer here
}  GFLIB_DTM_INFO;

// ----------------------------------------------------------------------
//    dec freq info
//
typedef struct {
	int	data;		// cal on
	int	stat;		// cal off
} GFLIB_DFRQ_COLNUMS;

typedef struct {
	GFLIB_DFRQ_COLNUMS colNums; // column numbers in fits file
	int	    frqChnSum; // freq Chan  added freq decimation
	int		*pdata;		  // malloc  data buffer here
	GFLIB_PDEVSTAT *pstat;  // malloc stat buffer here
}  GFLIB_DFRQ_INFO;

// ----------------------------------------------------------------------
// General prog info
//
typedef struct {
	int	crval1;		//  center freq band (hz)
	int	cdelt1;		//  channel width (hz). < 0 --> band flipped
	int	crval2;		//  ra j2000 start of row (deg)
	int	crval3;     //  dec j2000 start of row (deg)
	int mjd_obs;    // mjd start of row
} GFLIB_COLNUMS;

typedef struct {
	fitsfile *fptrIn;      // used by fits i/o 
	int	    status;		   // used by fits i/o

    int     naxis2;		   // number of rows in the file
	char	gfVersion[20]; // 1.0,1.1 ...
	int		decTmFile;	   // true decimateTime, false decFreq
	int	    shorts1Stat;   // 10

//  the following are values after decimation.

    int     nchn;      // each spectra    .. 4096
    int     npol;      // each summed rec .. 4
    int     nsmpRow;   // number summed time samples 1 row ..3
	int	    nelemRowData; // nchan*npol*nsmp
	int	    nelemRowStat; // nshort1stat*nsmp
//  use the following to malloc input buffers
    int     bytesData; // for calOn or calOff 1 row..
    int     bytesStat; // for calOn or calOff 1 row..
	int		bytesDataType;// 4 byte ints
	int		bytesStatType;// 1 byte ushort

	int	    dataTypeCode; // used when reading data
	int	    statTypeCode; // used when reading status
	GFLIB_COLNUMS colNums;// general col numbers..
// 
//  store decimate time files and decimate in freq files separately
//
//
	GFLIB_DTM_INFO  dtmI;    // info for decimate in time
	GFLIB_DFRQ_INFO dfrqI;  // info for decimate in freq
//
// Hold some of the header info read each row
//
 	double	crval1;			// cfr hz
 	double	cdelt1;			// channel width hz. (< 0 --> band flipped)
 	double	crval2;			// ra J2000 start of row (deg)
 	double	crval3;			// decJ2000 start of row (deg)
 	double	mjd_obs;		// start of first spectra this row
	} GFLIB_INFO;
/*
 * prototypes
*/
int gfLclose(GFLIB_INFO *pgfLI);
void gfLfits_err(char *msg,GFLIB_INFO *pgfLI);
int  gfLopen(char *filename,GFLIB_INFO *pgfLI);
int  gfLtdmRead(GFLIB_INFO *pgfLI,int row);
#endif


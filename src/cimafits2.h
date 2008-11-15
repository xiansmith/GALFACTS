/*
File: spec.h
The data file format produced from an Arecibo observation consist of two key files:
.cfg and .spec.  The .cfg contains the meta-data about the observation and is plain text.
The .spec is a binary file that contains the spectral and pointing information.

Reference is made to the FORTRAN definition of the code that procuces these data files:
http://www2.naic.edu/~tghosh/a1947/continuum_pointing.inc

*/

#include "common.h"

#ifndef FITSHEADERH
#define FITSHEADERH
#include    "fits_des.h"
/* AV20041118 $Id$
 *
 * 2005/02/03   lerner
 * Modified a number of comments
 *
 * 2005/02/03   lerner
 * Removed restfrqv, cdelt1v, crpix1v, and specsysv, renamed crval1v to
 * req_vel, renamed cunit1v to req_vel_unit, renamed ctype1v to req_vel_type,
 * added restfrqg, added specsysg, added vel_bary and vel_geo and modified
 * some comments
 *
 * 2005/02/01   lerner
 * Changed lst from hours to seconds, renamed elevation to elevatio, renamed
 * crval2 to req_ra and crval3 to req_dec, renamed equinox to req_equinox and
 * radesys to req_radesys, renamed crval2c to req_raj and crval3c to req_decj,
 * renamed crval2a to crval2 and crval3a to crval3, added restfrq and specsys,
 * changed croff2 from hours to degrees, renamed croff2 to off_ra and croff3
 * to off_dec, renamed croff2b to off_az and croff3b to off_za and modified
 * some comments
 *
 * 2005/01/24   lerner
 * Renamed cur_tol to cur_err and req_tol to allowed_err, renamed specsys to
 * specsysv, renamed crval2b to azimuth and crval3b to elevation, split up tcal
 * in tcal_frq and tcal_val, renamed obs_name to scantype, renamed pattern_scan
 * to pattern_id, renamed scan_number to scan_id, renamed pattern_number to
 * subscan, renamed total_pattern to total_subscans, added az_err and el_err,
 * added tdim1, added radesys, changed enc_time and off_time from AST to UTC
 * and modified some comments
 *
 * 2005/01/17   lerner
 * Changed unit 's' -> 'arcsec' for cur_tol and req_tol, renamed gain to
 * ampgain, changed crval2a and crval2c from h to deg, removed units from
 * crval1v and cdelt1v, added keyword bad_data, changed enc_time from int to
 * double, changed format for lags_in, standardised unit descriptions for
 * mjdxxobs, lst, croff2, rate_c1, rate_c2 and enc_time, rearranged some
 * keywords and modified some comments
 *
 *
 * $Log$
 * Revision 1.2  2008/11/15 09:11:53  sguram
 * *** empty log message ***
 *
 * Revision 1.1  2008/11/14 06:48:33  sguram
 * *** empty log message ***
 *
 * Revision 1.4  2004/12/19 03:01:58  arun
 * ctype1v, specsys TDISP spec per ML
 *
 * Revision 1.3  2004/12/10 04:40:31  arun
 * cunit1v, crpix1v field defs
 *
 * Revision 1.2  2004/12/02 03:50:40  arun
 * Fixed units in header per MN
 * NOTE NEED TO RUN make_fits_init
 *
 * Revision 1.1  2004/11/18 16:36:18  arun
 * Initial revision
 *
 */

/*  Note that the 'comment' for each entry will be parsed by 'make_fits_init'
    when building 'fits_init.c' - it is thus critical to keep the comments
    in the correct format! Also make sure that everything is properly
    8-byte aligned!  */

#define MAIN_HEADER 2880
#define BINTABLE_HEADER 20*2880
#define DUMPS_PER_ROW 500
#define RAW_CHANNELS 2048
/*  Comment-format:       unit format comment */
typedef struct  {
  FITS_ARRAY_DES datapointer; /* x x Pointer to the data in the heap */
  FITS_ARRAY_DES statpointer; /* x x Ptr to status for each spectraj */
  char tdim1[32];      /* x A32 Dimensions of data pointed to in the heap */
  char tdim2[32];      /* x A32 Dimensions of status pointed to in the heap */
  char object[16];     /* x A16 Name of source observed */
//COMMENT axis 1 is the frequency axis
  double crval1;       /* Hz D13.5 Center frequency */
  double cdelt1;       /* Hz D13.5 Frequency interval */
  double crpix1;       /* x D13.5 Pixel of center frequency */
//COMMENT CRVAL2,3 are the measured J2000 positions on the sky
//COMMENT For arrays they are the start of the first pixel
  double crval2;       /* deg D13.5 Actual RA J2000 pointing this beam on sky */
  double crval3;       /* deg D13.5 Actual DEC J2000 pointing this beam on sky */
  double crval4;       /* x D13.5 Polarization (neg -> Pol, pos -> Stokes) */
//COMMENT for multi spectra/row crval5 is the start of the 1st spectra
  double crval5;       /* s D13.5 Seconds since midnight from obsdate (UTC) */
  double cdelt5;	   /* s D13.5 time difference (wall) between spectra in a row */
//COMMENT az,el,glon,glat are alternate coord positions for crval2,3
//COMMENT azimuth,elevation do not contain model corrections.
  double azimuth;      /* deg D13.5 Actual AZ pointing this beam on sky   */
  double elevatio;     /* deg D13.5 Actual EL pointing this beam on sky   */
  double glon  ;      /* deg D13.5 Actual galactic l pointing this beam on sky */
  double glat;        /* deg D13.5 Actual galactic b pointing this beam on sky */
  char datexxobs[16];  /* x A16 Start of this observation (UTC) - YYYY-MM-DD */
  double mjdxxobs;     /* d D13.5 MJD number at exposure start of 1st spc of row */
  double lst;          /* s D13.5 Local mean sidereal time 1st spc of row*/
  double exposure;     /* s D13.5 Exposure (integration) time of current record (dump) */
  double duration;    /* s D13.5 Duration of the dump (Wall time)*/
  double tsys;         /* K D13.5 Last computed Tsys - set to 1.0 if unavailable */
  double bandwid;      /* Hz D13.5 Overall bandwidth of spectrum */
  double restfrq;      /* Hz D13.5 Rest frequency used for doppler correction */
  char specsys[8];      /* x A8 Velocity frame for CRVAL1 (always TOPOCENT for AO data) */
//COMMENT req_ are parameters REQUESTEd by the observer.
  char req_sys[8];      /* x A8 User-requested velocity frame */
  double req_vel;      /* x D13.5 Requested velocity or z in frame REQ_SYS */
  char req_vel_unit[8]; /* x A8 Specifies units of REQ_VEL: either m/s or Z */
  char req_vel_type[8]; /* x A8 Velocity type for REQ_VEL */
//COMMENT req_raj/decj is the requested resultant position on the sky 
//COMMENT it includes of _pos,_off, and _rate 
  double req_raj;      /* deg D13.5 Requested RA J2000 position */
  double req_decj;     /* deg D13.5 Requested Dec J2000 position */
//COMMENT _pos,_off, and _rate are used for mapping. 
//COMMENT _pos is the center,_off is the offset, and _rate is the rate
//COMMENT req_raj/decj are the combination of the three.
  double req_posc1;    /* deg D13.5 Requested long (usually RA) in REQ_COORDSYS */
  double req_posc2;    /* deg D13.5 Requested lat (usually DEC) in REQ_COORDSYS */
  char   req_cspos[8];   /* x A8 Coordinate system used for REQ_POSC1 and REQ_POSC2 */
  double req_offc1;     /* deg D13.5 requested offset from posc1*/
  double req_offc2;     /* deg D13.5 requested offset from posc2 */
  char   req_csoff[8];   /* x A8 Coordinate system used for req_offc1,2 */
  double req_pnttime;  /* s D13.5 UTC secMidnite. tmStamp uninterpolated pntData*/
  double req_ratec1;    /* deg/s D13.5 requested rate for c1 */
  double req_ratec2;    /* deg/s D13.5 requested rate for c2 */
  char   req_csrate[8]; /* x A8 Coordinate system used for req_ratec1,2 */
  double req_equinox;  /* x D13.5 Equinox of  REQ_POSC1 and REQ_POSC2 */
  double rate_dur;     /* s D13.5 How long has req_ratec1,2  been applied        */
//COMMENT enc_ are the raw encoder values read at enc_time
  double enc_time;     /* s D13.5 Time when encoders were read out (UTC) */
  double enc_azimuth;  /* deg D13.5 Azimuth encoder read-out at ENC_TIME */
  double enc_elevatio; /* deg D13.5 Elevation encoder read-out at ENC_TIME */
  double enc_altel;    /* deg D13.5 Elevation encoder of other carriage house */
//COMMENT cur_err is current great circle tracking error at encoder
  double cur_err;      /* deg D13.5 Actual great circle tracking error */
  double allowed_err;  /* deg D13.5 Maximum allowed tracking error   */
  double az_err;       /* deg D13.5 Azimuth tracking error (actualPos-requested) */
  double el_err;       /* deg D13.5 Elevation tracking error (actualPos-requested) */
//COMMENT model_offxx are the model offsets for the measured position.
//COMMENT user_offxx are additional offset requested by the user
//COMMENT they should be added to model_offxx 
//COMMENT Sign is: encPosUsed=computedPos + xxOffset
  double model_offaz;  /* deg D13.5 Pointing model offset AZ at curpos  */
  double model_offza;  /* deg D13.5 Pointing model offset ZA at curpos  */
  double user_offaz;   /* deg D13.5 User selected pointing offset AZ (great circle) */
  double user_offza;   /* deg D13.5 User selected pointing offset ZA    */
  double beam_offaz;   /* deg D13.5 ALFA unrotated offset AZ            */
  double beam_offza;   /* deg D13.5 ALFA unrotated offset ZA            */
// ?? what is this.. the great circle offsets pnt cor setoff az,za ??
//    is it included in the requested positions ??
  double rfeed_offaz;  /* deg D13.5 Rotated offset this beam AZ         */
// ?? from center of alfa or center of prfeed.. comment should specify
  double rfeed_offza;  /* deg D13.5 Rotated offset this beam ZA         */
  double prfeed_offaz; /* deg D13.5 Offset to center prfeed beam AZ     */
  double prfeed_offza; /* deg D13.5 Offset to center prfeed beam ZA     */
  double beam_offraj;  /* deg D13.5 Total RA offset to this beam        */
  double beam_offdecj; /* deg D13.5 Total DEC offset to this beam       */
//COMMENT the map offsets are the measured offsets from req_raj,req_decj positions
  double map_offra;   /* deg D13.5 Actual RA J2000 offset req_raj/decj */
  double map_offdec;  /* deg D13.5 Actual DEC J2000 offset to req_raj/decj*/
  double map_offaz;   /* deg D13.5 Actual AZ offset to req_raj/decj (great circle?)*/
  double map_offza;   /* deg D13.5 Actual ZA offset to req_raj/decj */
  double alfa_ang;    /* deg D13.5 ALFA rotation angle */
  double para_ang;    /* deg D13.5 Parallactic angle */
  double vel_bary;    /* m/s D13.5 Projected barycentric velocity (incl. VEL_GEO) */
  double vel_geo;     /* m/s D13.5 Projected geocentric velocity */
  double vel_obs;	/* m/s D13.5 Observers projected velocity in req_sys*/
  char frontend[8];    /* x A8 Receiver name */
  char backend[8];	/* x A8 Backend name */	
  char backendmode[24]; /* x A24 Backend mode description */
  char caltype[8];     /* x A8 diode calibration mode */
  char obsmode[8];     /* x A8 Name of observation pattern (e.g. ONOFF) */
  char scantype[8];    /* x A8 Type of scan (as part of pattern - e.g. ON OFF) */
  int pattern_id;      /* x I9 Unique number for observation pattern YDDDnnnnn */
  int scan_id;         /* x I9 Unique number for scan YDDDnnnnn */
  int recnum;         /* x I8 Sequential number of current record (dump) */
  int total_recs;     /* x I8 Requested total number of records (dumps) in scan */
  int chan_in;         /* x I8 Number of freq channels input this spc*/
  unsigned int bad_data; /* x I8 0->good data <>0->problem (see COMMENT) */
  double plat_powerA;   /* x D13.5 platform power meter reading polA */
  double plat_powerB;   /* x D13.5 platform power meter reading polb */
  double cntrl_powerA;  /* x D13.5 control room power meter polA */
  double cntrl_powerB;  /* x D13.5 control room power meter polB */
  double syn1;         /* Hz D13.5 Platform synthesizer */
  double syn2;        /* Hz D13.5 2nd lo synth value this band */
  double tot_power;   /* x D13.5 Scaled power in zero-lag. -1 for fft bkends */
  int    tcal_numCoef;   /* x I8 Number of coef. in polyFit for tcal */
  int    fill1;       /* x I8 for 8byte alignment */
  double tcal_coefA[4]; /* GHz D13.5 Polyfit tcal polA. order(deg) 0,1,2,3 */
  double tcal_coefB[4]; /* GHz D13.5 polyfit tcal polB */
  unsigned char backendmask[8]; /* x B1 Which backend boards enabled */
// ?? how is this coded
//    pdev has have 7 beams and with two fpgas/beam. you can use 
//    the highband, lowband or both for each beam.
// --
    unsigned char num_beams;   /* x B1 Number of beams in this observation (1, 7 or 8) */
    unsigned char num_ifs;     /* x B1 Number of IFs for this beam (1 - 8 or 1 - 2) */
//  ?? is an IF a subband (section of frequency space?)
//     - We don't know how many IFs there are for this
//       fits file. Will we always have 1 IF per fits file
//       with pdev (even in single pixel?)
    unsigned char num_pols; /* x B1 Number of pols for this IF and this beam (1, 2 or 4) */
//  ?? 1,2,4 how does this differ from crval4
    unsigned char beam;        /* x B1 Number of this beam (1, 0 - 6 or 0 - 7) */
    unsigned char ifn;	      /* x B1 Number of this IF (1 - 8 or 1 - 2) */
//  ?? if is a reserved identifier and may not be used
    unsigned char pol;	     /* x B1 Number of this pol (1, 2 or 4) */
//  ??  should be 1,2,3,4
  unsigned char prfeed;   /* x B1 ALFA beam used as pointing center */
  unsigned char input_id; /* x B1 spectrometer board number (engineering parameter) */
// --
// ?? how will this be numbered for pdev
//    will high,low bands be different numbers 
  unsigned char uppersb;  /* x B1 True if spectrum flipped */
// ?? does the sign of cdelt1  tell us the same info?
  unsigned char attn_cor; /* x B1 Correlator attenuator: 0-15 dB */
// ?? pdev also has an attenuator setting
//    should the be changed to attn_bkend or something
  unsigned char master;   /* x B1 0=Gregorian dome 1=Carriage house */
  unsigned char onsource; /* x B1 True if on-source at ENC_TIME */
  unsigned char blanking; /* x B1 Blanking turned on */
  unsigned char lbwhyb;   /* x B1 LBandWide Hybrid is in (for circular pol) */
  unsigned char shcl;     /* x B1 True if receiver shutter closed */
  unsigned char sbshcl;   /* x B1 True if S-band receiver shutter closed */
// --
  unsigned char rfnum;    /* x B1 Platform position of the receiver selector */
  unsigned char calrcvmux; /* x B1 Platform cal selector */
// ?? do we have the caltype highcal correlatored cal etc somewhere
  unsigned char zmnormal; /* x B1 Platform transfer switch to reverse channels, true normal */
  unsigned char rfattn[2];  /* x B1 Platform attenuator positions */
  unsigned char if1sel; /* x B1 Platform IF selector, 1/300 2/750, 3/1500, 4/10GHz1500, 5-thru */
  unsigned char ifattn[2];  /* x B1 Platform IF attenuator positions */
// --
  unsigned char fiber;      /* x B1 True if platform fiber chosen (usually true) */
// ??  remove.. bill removed ability to switch to coax.
  unsigned char ac2sw;      /* x B1 Platform AC power to various instruments etc. */
  unsigned char phbsig;     /* x B1 Platform converter combiner signal ph adjust */
  unsigned char hybrid;     /* x B1 Platform converter combiner hybrid */
  unsigned char phblo;      /* x B1 Platform converter combiner LO phase adjust */
  unsigned char xfnormal;   /* x B1 Control room transfer switch true = default */
  unsigned char ampgain[2]; /* x B1 Gain of control room amplifiers */
// --
  unsigned char noise;      /* x B1 Control room noise on */
  unsigned char inpfrq;     /* x B1 Control room input distributor position */
  unsigned char mixer;   /* x B1 Control room mixer source switches */
  unsigned char vlbainp;    /* x B1 Control room VLBA input switch position */
  unsigned char syndest; /* x B1 Control room synthesizer destination for this board*/
// ?? ampinp no longer needed. determined by inpfrq
  unsigned char calsrc;     /* x B1 Control room cal source bit */
  unsigned char cal;        /* x B1 Is cal bit turned on */
  unsigned char vis30mhz;   /* x B1 Control room greg 1 ch 0 */
// --
//  unsigned char pwrmet;     /* x B1 no longer used */
  unsigned char blank430;   /* x B1 Control room 430 blanking on */
  unsigned char fill[7];    /* x B1 Round up header to even 8 bytes */
} cimafits_row;
#endif

/* Status heap structure */
typedef struct {
	short int seqNum;
	short int fftAccum;
	short int calOn;
	short int adcOverFlow;
	short int pfbOverFlow;
	short int satCntVShift;
	short int satCntAccS2S3;
	short int satCntAccS0S1;
	short int satCntAshftS2S3;
	short int satCntAshftS0S1;
}pdev_stat;

/*Data help struture */
typedef struct {
	unsigned short int pol_A[RAW_CHANNELS];
	unsigned short int pol_B[RAW_CHANNELS];
	short int stokes_U[RAW_CHANNELS];
	short int stokes_V[RAW_CHANNELS];
}pdev_datum;

/* 
Object to contain the data from the .cfg file.
Describes the metadata from the observation.
*/
typedef struct {
	float integrationTime; //miliseconds
	int stokesSetSize; //#bytes per stokes set
	float centerMHz; //MHz
	float bandwitdhkHz; //kHz
	int startChanNum; //start channel number
	int samplesPerStokesSet; //no. of samples in a Stokes set
	int stokesProducts; //no. of Stokes products
	int numStokesSets; //and no of such sets corresponding a given integration (i.e. cal_on,cal_off)
	int crap; //undocumented
	char observingTag[32]; //Source-name or observing-mode tag, WCAL tag to indicate winking cal data 
	int day; //date dd mm yyyy
	int month;
	int year;
	int hour; // AST start time hh mm ss.ssssss
	int minute;
	float second;
	char observatoryCode[16]; //Observatory code
	int bandFlip; //flag to indicate band-flip (-1 does indicate band inversion, +1: no flip)
} ConfigData;


/* 
CentralBeamBlock
Contains the pointing information for the primary beam, (beam0), along with 
some common info that is general to the position of the receiver on the sky,
as opposed to a particular beam.
*/
typedef struct {
    double raj_true_in_degrees;
    double decj_true_in_degrees;
    double arecibo_local_mean_sidereal_time_in_sec;
    double mjd_last_five_digits;
    double az_cur_in_degrees;
    double za_cur_in_degrees;
    double raj_requested_in_degrees;
    double decj_requested_in_degrees;
    double alfa_rotation_angle_in_degrees;
} nSpecPointingBlock;

typedef struct {
	unsigned  int *A;
	unsigned  int *B;
	int *U;
	int *V;
}Spec_PolSet;

static inline void cnvrt_end_sint(short int *x);
static inline void cnvrt_end_int(int *x);
static inline void cnvrt_end_db(double *x);

/*
SpecRecord
This is the core data structure for operating on data in the time domain: ie
calibration.  RA,DEC,AST is the pointing information.  calon and caloff is the
raw spectral data, cal is the computed value of the cal signal.  flagRFI
contain the RFI bit flags on a per channel basis.  flagBAD is non-zero when
the entire spectra is marked bad.  stokes contains the spectra of the computed,
calibrated stokes parameters.
*/


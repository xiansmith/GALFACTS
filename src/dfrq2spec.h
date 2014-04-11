#include "galfactsLib.h"
#define MAIN_HEADER 2880
#define BINTABLE_HEADER 20*2880
#define RAW_CHANNELS 128
#define DUMPS_PER_ROW 600
#define LINELEN 80
typedef struct
{
	int A[RAW_CHANNELS];
	int B[RAW_CHANNELS];
	unsigned int U[RAW_CHANNELS];
	unsigned int V[RAW_CHANNELS];
}GFLIB_DATA;

typedef struct
{
	GFLIB_DATA data[DUMPS_PER_ROW];
	GFLIB_PDEVSTAT stat[DUMPS_PER_ROW];
	char fill11[64];
	char object[16];
	double cf; // center frequency
	double fdelt; //delta frequency
	double crpix;
	double RA;
	double DEC;
	double crval4; //?? not used
	double UTC;
	double tdelt; //time delta ?? integration time per point ??
	double azimuth;      /* deg D13.5 Actual AZ pointing this beam on sky   */
	double elevatio;     /* deg D13.5 Actual EL pointing this beam on sky   */
	double glon  ;      /* deg D13.5 Actual galactic l pointing this beam on sky */
	double glat;        /* deg D13.5 Actual galactic b pointing this beam on sky */
	char datexxobs[16];  /* x A16 Start of this observation (UTC) - YYYY-MM-DD */
	double mjdxxobs;     /* d D13.5 MJD number at exposure start of 1st spc of row */
	double lst;
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
	  unsigned char fill[7];
	//char fill[696]; //just read in the rest of the row, values currently not used for anything
}GFLIB_ROW;

static inline void cnvrt_end_sint(short int *x);
static inline void cnvrt_end_int(int *x);
static inline void cnvrt_end_db(double *x);


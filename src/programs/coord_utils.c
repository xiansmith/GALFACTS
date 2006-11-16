/*****************************************************************************/
/**									    **/
/**	COORD_UTILS.C							    **/
/**									    **/
/**	A collection of astronomical coordinate manipulation functions.	    **/
/**									    **/
/*****************************************************************************/

/** functions contained in this collection:

 double ang_sep (ra1, dec1, ra2, dec2)
 double ang_sep_err (ra1, ra1_sig, dec1, dec1_sig, ra2, ra2_sig, dec2,dec2_sig)
 double arc_sep (ra1, dec1, ra2, dec2, ra0, dec0)
 double arc_sep_err (ra1, ra1_err, dec1, dec1_err, ra2, ra2_err, dec2, ... )
   void get_xi_eta (ra0, dec0, ra, dec, xi_ptr, eta_ptr)
   void get_ra_dec (ra0, dec0, xi, eta, ra_ptr, dec_ptr)
   void precess (ra0, dec0, ra1_ptr, dec1_ptr, from_equinox, to_equinox)
   void precess_2000_1950 (ra0, dec0, ra1_ptr, dec1_ptr)
   void precess_1950_2000 (ra0, dec0, ra1_ptr, dec1_ptr)
   void precess_1950_1900 (ra0, dec0, ra1_ptr, dec1_ptr)
   void precess_1900_1950 (ra0, dec0, ra1_ptr, dec1_ptr)
   void eq1950_to_gal (ra1950, dec1950, glon_ptr, glat_ptr)
   void eq2000_to_gal (ra2000, dec2000, glon_ptr, glat_ptr)
   void gal_to_eq1950 (glon, glat, ra1950_ptr, dec1950_ptr)
   void gal_to_eq2000 (glon, glat, ra2000_ptr, dec2000_ptr)
   void ecliptic_to_equatorial (lambda, beta, ra_ptr, dec_ptr, equinox)
   void equatorial_to_ecliptic (ra, dec, lambda_ptr, beta_ptr, equinox)
 double ecliptic_obliquity (equinox)
   void read_coords (coordstr, ra_ptr, dec_ptr)
   void h2hms (h_real, h_int_ptr, m_int_ptr, s_real_ptr)
   void hms2h (h_int, m_int, s_real, h_real_ptr)
 double smart_atan (y, x)
   void limit_celestial_coords (clon_ptr, clat_ptr)
   void limit_celestial_coords_2 (clon_ptr, clat_ptr)
   void limit_ra (ra_ptr)

**/

/** The following definition files must be included prior to COORD_UTILS.C: 

	mathdefs.h
	misc_math.c
**/

/* Data format type flags, used by READ_COORDS and similar routines */

#define SEX_HRS		 1	/*  HH MM SS.SS  DD MM SS.S */
#define SEX_DEG		 2	/* DDD MM SS.SS  DD MM SS.S */
#define DEC_HRS		 3	/*  HH.HHHHHH    DD.DDDDD   */
#define DEC_DEG		 4	/* DDD.DDDDD     DD.DDDDD   */

/*upper limit angle for using approximation in angular separation computation*/
#define MAX_APPROX_ANGLE  (0.01 * D2RAD)	/* 0.01 deg */

/*****************************************************************************/
/** compute angular separation between two points, given RAs & DECs 	    **/
/** all angles are in radians 						    **/
/*****************************************************************************/
double ang_sep (ra1, dec1, ra2, dec2)
    double ra1, dec1, ra2, dec2;
{
  double sep, sep_approx;

  sep_approx = sqrt (cos(dec1)*cos(dec2)*(ra2-ra1)*(ra2-ra1) 
						    + (dec2-dec1)*(dec2-dec1));
  if (fabs (cos (sep_approx)) > cos (MAX_APPROX_ANGLE)) {
    /** use approximation -- angle very close to 0 or 180 **/
    sep = sep_approx;
  } else {
    /** use full forumlation **/
    sep = acos (sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra2-ra1));
  }
  return sep;
}

/*****************************************************************************/
/** compute angular separation ERROR between two points, given RA, DEC, err **/
/** all angles are in radians 						    **/
/*****************************************************************************/
double ang_sep_err (ra1, ra1_sig, dec1, dec1_sig, ra2, ra2_sig, dec2, dec2_sig)
  double ra1, ra1_sig, dec1, dec1_sig, ra2, ra2_sig, dec2, dec2_sig;
{
  double fsqr();
  double sep_approx, denominator, numerator, sep_sig;

  sep_approx = sqrt (cos(dec1)*cos(dec2)*(ra2-ra1)*(ra2-ra1) 
						+ (dec2-dec1)*(dec2-dec1));
  if (fabs (cos (sep_approx)) > cos (MAX_APPROX_ANGLE)) {
    /** use approximation -- angle very close to 0 or 180 **/
    denominator = fsqr (cos(0.5*(dec1+dec2))*(ra2-ra1)) + fsqr (dec2-dec1);
    numerator = (fsqr (ra1_sig) + fsqr (ra2_sig)) 
			* fsqr(fsqr(cos(0.5*dec1+dec2))*(ra2-ra1))
		+ fsqr (dec1_sig * 
			(0.25*sin(dec1+dec2)*fsqr(ra2_sig-ra1_sig)+dec1-dec2))
		+ fsqr (dec2_sig * 
			(0.25*sin(dec1+dec2)*fsqr(ra2_sig-ra1_sig)-dec1+dec2));
  } else {
    /** use full forumlation **/
    denominator = 1.0 - fsqr (sin(dec1)*sin(dec2) 
				+ cos(dec1)*cos(dec2)*cos(ra1-ra2));
    numerator = (fsqr (ra1_sig) + fsqr (ra2_sig))
			* fsqr (cos(dec1) * cos(dec2) * sin(ra1-ra2))
		+ fsqr (dec1_sig *
			(cos(dec1)*sin(dec2)-sin(dec1)*cos(dec2)*cos(ra1-ra2)))
		+ fsqr (dec2_sig *
		       (sin(dec1)*cos(dec2)-cos(dec1)*sin(dec2)*cos(ra1-ra2)));
  }
  if (denominator != 0.0) {
    sep_sig = sqrt (numerator / denominator);
  } else {
    printf ("#%%% Warning: denominator == 0 in ang_sep_err! %%%\n");
    sep_sig = 0.0;
  }
  return sep_sig;
}

/*****************************************************************************/
/** Compute "arc separation" of two points in reference to a third, ie.,    **/
/** the angle between lines drawn from the reference to each of the two     **/
/** points, measured at the reference point.  A positive sign indicates     **/
/** that (ra1,dec1) is counterclockwise of (ra2,dec2). 			    **/
/** all angles are in radians 						    **/
/*****************************************************************************/
double arc_sep (ra1, dec1, ra2, dec2, ra0, dec0)
  double ra1, dec1, ra2, dec2, ra0, dec0;
{
  double ang_sep(), smart_atan();
  double sep, rad1, rad2, asep, asep_cos, flatarc1, flatarc2, diff;

  sep  = ang_sep (ra1, dec1, ra2, dec2);
  rad1 = ang_sep (ra1, dec1, ra0, dec0);
  rad2 = ang_sep (ra2, dec2, ra0, dec0);

  if (rad2 == 0.0) {
    asep = rad1;
  } else {
    if (rad1 == 0.0) {
      asep = rad2;
    } else {
      asep_cos = (cos(sep) - cos(rad1)*cos(rad2)) / (sin(rad1)*sin(rad2));
      if (fabs (asep_cos) <= 1.0) {
	asep = acos (asep_cos);
	/** determine whether angle is positive or negative **/
	flatarc1 = smart_atan (dec1-dec0, ra0-ra1);
	flatarc2 = smart_atan (dec2-dec0, ra0-ra2);
	diff = flatarc1 - flatarc2;
	if (diff > PI) diff -= 2.0*PI;
	if (diff < 0.0) {
	  asep *= -1.0;
	}
      } else {
	if (asep_cos < -1.0) asep = acos (-1.0);
	if (asep_cos >  1.0) asep = acos ( 1.0);
      }
    }
  }
  return asep;
}

/*****************************************************************************/
/** Compute "arc separation error". 					    **/
/**									    **/
/** NOTE: The errors produced by this routine are suspect.                  **/
/** They seem too small by at least an order of magnitude.                  **/
/** The assumptions used in their derivation may not have been appropriate. **/
/**									    **/
/** all angles are in radians						    **/
/*****************************************************************************/
double arc_sep_err (ra1, ra1_err, dec1, dec1_err, ra2, ra2_err, dec2, dec2_err,
								ra0, dec0)
  double ra1, ra1_err, dec1, dec1_err, ra2, ra2_err, dec2, dec2_err, ra0, dec0;
{
  void reverse_video(), normal_video(), wait_for_CR();
  double ang_sep(), ang_sep_err(), fsqr();
  double sep, rad1, rad2, sep_sig, rad1_sig, rad2_sig,
	asep_sig, flatarc1, flatarc2, diff, denominator, numerator,
	p, pfunc, dasep_drad1, dasep_drad2, dasep_dsep;

  sep  = ang_sep (ra1, dec1, ra2, dec2);
  rad1 = ang_sep (ra1, dec1, ra0, dec0);
  rad2 = ang_sep (ra2, dec2, ra0, dec0);
  sep_sig  = ang_sep_err (ra1,ra1_err,dec1,dec1_err,ra2,ra2_err,dec2,dec2_err);
  rad1_sig = ang_sep_err (ra1,ra1_err,dec1,dec1_err,ra0,  0.0,  dec0,  0.0);
  rad2_sig = ang_sep_err (ra2,ra2_err,dec2,dec2_err,ra0,  0.0,  dec0,  0.0);
  if (rad1 == 0.0 || rad2 == 0.0 || sep == 0.0) {
    printf ("#%%% Warning: zero separation detected in arc_sep_err! %%%\n");
    asep_sig = 0.0;
  } else {

    /****** REVISED FORMULATION BASED UPON PLANAR APPROXIMATION ************/

    p = (fsqr(rad1) + fsqr(rad2) - fsqr(sep)) / (2.0 * rad1 * rad2);
    numerator =   fsqr (rad1_sig * (1.0/rad2 - p/rad1))
		+ fsqr (rad2_sig * (1.0/rad1 - p/rad2))
		+ fsqr (sep_sig  * sep / (rad1 * rad2));
    denominator = 1.0 - fsqr(p);

    /***********************************************************************/

    if (denominator != 0.0) {
      asep_sig = sqrt (numerator / denominator);
    } else {
      printf ("#%%% Warning: denominator == 0 in arc_sep_err! %%%\n");
      asep_sig = 0.0;
    }
  }
  return asep_sig;
}

/*****************************************************************************/
/** compute (Xi, Eta) gnomonic projection coordinates from (RA, DEC)        **/
/** Equations in this and the following routine are taken from              **/
/**     Spherical Astronomy, Robin M. Green, 1985, Cambridge U. Press       **/
/** and Textbook on Spherical Astronomy / 6e, W. M. Smart, 1977, Cambridge  **/
/** all angles are in radians                                               **/
/*****************************************************************************/
void get_xi_eta (ra0, dec0, ra, dec, xi_ptr, eta_ptr)
    double ra0, dec0, ra, dec, *xi_ptr, *eta_ptr;
{
    double denominator;

    denominator = (sin(dec0)*sin(dec) + cos(dec0)*cos(dec)*cos(ra-ra0));
    *xi_ptr = cos(dec)*sin(ra-ra0) / denominator;
    *eta_ptr = (cos(dec0)*sin(dec) - sin(dec0)*cos(dec)*cos(ra-ra0)) 
		/ denominator;
}

/*****************************************************************************/
/** compute (RA, DEC) from (Xi, Eta) gnomonic projection coordinates 	    **/
/** all angles are in radians 						    **/
/*****************************************************************************/
void get_ra_dec (ra0, dec0, xi, eta, ra_ptr, dec_ptr)
    double ra0, dec0, xi, eta, *ra_ptr, *dec_ptr;
{
    int sign_xi = 1, sign_tan_rel_ra = 1;
    double denominator, tan_rel_ra, rel_ra;

    /* not sure how to handle sign ambiguities here, but it looks like there
	won't be a problem for the Pleiades mosaic coordinate range */

    denominator = cos(dec0) - eta*sin(dec0);
    tan_rel_ra = xi / denominator;
    if (xi < 0) sign_xi = -1;
    if (tan_rel_ra < 0) sign_tan_rel_ra = -1;
    if (sign_xi != sign_tan_rel_ra) {
/**/
	printf ("%%% Warning: sign_xi != sign_tan_rel_ra %%%\n");
/**/
    }
    rel_ra = atan (tan_rel_ra);
    *ra_ptr = rel_ra + ra0;
    *dec_ptr = atan ((sin(dec0) + eta*cos(dec0)) / (cos(dec0) - eta*sin(dec0))
		     * cos (rel_ra));
}

/*****************************************************************************/
/* Precess RA & DEC from one of (1900, 1950, 2000) to something else.	     */
/* Specify both equinoxes after coordinate pairs in argument list.	     */
/* 									     */
/* A kludge is used to avoid NaN RA's near poles, but even so, this	     */
/* formulation is not trustworthy there.				     */
/* 									     */
/* 	Formulae & Constants in this and other precession routines, plus     */
/*	conversions between equatorial, galactic, and ecliptic coordinates,  */
/*	have been taken directly from Duffett-Smith (1981), 		     */
/*	_Practical Astronomy with Your Calculator_, section 33, pages 58-59  */
/*	They are not as accurate as the rigorous methods espoused in         */
/*	the _Astronomical Almanac_, but should do for most purposes.         */
/* 									     */
/* 									     */
/* 	THIS ROUTINE HAS BEEN OUTMODED.  IT PRODUCES ASYMMETRIC RESULTS	     */
/* 	AT THE 10" LEVEL IN THE PLEIADES FOR 1950 --> 2000 --> 1950.	     */
/* 	IT HAS BEEN REPLACED WITH A BETTERN ALGORITHM DESCRIBED IN THE	     */
/* 	_ASTRONOMICAL ALMANAC_.  [1999 OCTOBER 14]			     */
/* 									     */
/* All coordinates & angles input & output in RADIANS.			     */
/*****************************************************************************/
void precess_old (ra0, dec0, ra1_ptr, dec1_ptr, from_equinox, to_equinox)
  int from_equinox, to_equinox;
  double ra0, dec0, *ra1_ptr, *dec1_ptr;
{
  void exit(), limit_celestial_coords();
  double m, n, nprime, big_N, dec_fudge;

  switch (from_equinox) {
    case (1900) : m = 3.07234; n = 1.33646; nprime = 20.0468; break;
    case (1950) : m = 3.07327; n = 1.33617; nprime = 20.0426; break;
    case (2000) : m = 3.07420; n = 1.33589; nprime = 20.0383; break;
    default: 
      fprintf (stderr, "conversion from equinox %d not supported\n",
								from_equinox);
      exit(1);
  }
  big_N = 1.0 * (to_equinox - from_equinox);

  /* limit approach of input DEC to poles to prevent tan(90) errors */
  dec_fudge = 0.1*PI/180.0;
  if ( 0.5*PI - dec0 < dec_fudge) dec0 = 0.5*PI - dec_fudge;
  if ( 0.5*PI + dec0 < dec_fudge) dec0 = dec_fudge - 0.5*PI;

  *ra1_ptr  = ra0 + (m + n * sin(ra0) * tan(dec0)) * big_N * PI/(12.0*3600.0);
  *dec1_ptr = dec0 + (nprime * cos(ra0)) * big_N * PI/(180.0*3600.0);
  limit_celestial_coords (ra1_ptr, dec1_ptr);
}

/*****************************************************************************/
/* Precess RA & DEC either to or from an equinox of J2000.           	     */
/* A shell routine calls this twice to transform between arbitrary dates.    */
/* 									     */
/* 	Formulae & Constants are from the 1998 _Astronomical Almanac_,       */
/*	page B19.  These are the "approximate" (non-matrix) formulae.        */
/*	They are properly symmetric, unlike those of Duffett-Smith (1981),   */
/*	which are discrepant in 1950 -> 2000 -> 1950 transformations at      */
/*	the 10" level.  The _Almanac_ gives no indication of the level of    */
/*	accuracy of these formulae, but presumably it is higher.             */
/* 									     */
/* All coordinates & angles input & output in RADIANS.			     */
/*****************************************************************************/
#define TO_2000    1
#define FROM_2000 -1
void precess_2000 (ra0, dec0, ra1_ptr, dec1_ptr, equinox, FROM_TO_FLAG)
  int equinox, FROM_TO_FLAG;
  double ra0, dec0, *ra1_ptr, *dec1_ptr;
{
  void exit(), limit_celestial_coords();
  double dec_fudge, T, M, N, ra_m, dec_m;

  /* limit approach of input DEC to poles to prevent tan(90) errors */
  dec_fudge = 0.1 * AS2RAD;
  if ( 0.5*PI - dec0 < dec_fudge) dec0 = 0.5*PI - dec_fudge;
  if ( 0.5*PI + dec0 < dec_fudge) dec0 = dec_fudge - 0.5*PI;

  T = 0.01 * (equinox - 2000);
  M = (1.2812323 * T + 0.0003879 * T*T + 0.0000101 * T*T*T) * D2RAD;
  N = (0.5567530 * T - 0.0001185 * T*T - 0.0000116 * T*T*T) * D2RAD;
  switch (FROM_TO_FLAG) {
    case (TO_2000) :
      ra_m  = ra0  - 0.5 * (M + N * sin (ra0) * tan (dec0));
      dec_m = dec0 - 0.5 * N * cos (ra_m);
      *ra1_ptr  = ra0  - M - N * sin (ra_m) * tan (dec_m);
      *dec1_ptr = dec0 - N * cos (ra_m);
      break;
    case (FROM_2000) : 
      ra_m  = ra0  + 0.5 * (M + N * sin (ra0) * tan (dec0));
      dec_m = dec0 + 0.5 * N * cos (ra_m);
      *ra1_ptr  = ra0  + M + N * sin (ra_m) * tan (dec_m);
      *dec1_ptr = dec0 + N * cos (ra_m);
      break;
    default : 
	fprintf (stderr, "precess_2000: illegal value: FROM_TO_FLAG = %d!\n",
		FROM_TO_FLAG);
	exit (1);
  }
  limit_celestial_coords (ra1_ptr, dec1_ptr);
}

/*****************************************************************************/
/* 	Calling shell for new & improved precession routine.		     */
/*****************************************************************************/
void precess (ra0, dec0, ra1_ptr, dec1_ptr, from_equinox, to_equinox)
  int from_equinox, to_equinox;
  double ra0, dec0, *ra1_ptr, *dec1_ptr;
{
  void precess_2000();
  double ra_2000, dec_2000;
  precess_2000 (ra0, dec0, &ra_2000, &dec_2000, from_equinox, TO_2000);
  precess_2000 (ra_2000, dec_2000, ra1_ptr, dec1_ptr, to_equinox, FROM_2000);
}

/*****************************************************************************/
/* 	Precess 2000.0 RA & DEC --> 1950.0.				     */
/* 	All coordinates & angles input & output in RADIANS.		     */
/*****************************************************************************/
void precess_2000_1950 (ra0, dec0, ra1_ptr, dec1_ptr)
  double ra0, dec0, *ra1_ptr, *dec1_ptr;
{
  void precess();
  precess (ra0, dec0, ra1_ptr, dec1_ptr, 2000, 1950);
}

/*****************************************************************************/
/* 	Precess 1950.0 RA & DEC --> 2000.0.				     */
/* 	All coordinates & angles input & output in RADIANS.		     */
/*****************************************************************************/
void precess_1950_2000 (ra0, dec0, ra1_ptr, dec1_ptr)
    double ra0, dec0, *ra1_ptr, *dec1_ptr;
{
  void precess();
  precess (ra0, dec0, ra1_ptr, dec1_ptr, 1950, 2000);
}

/*****************************************************************************/
/* 	Precess 1950.0 RA & DEC --> 1900.0.				     */
/* 	All coordinates & angles input & output in RADIANS.		     */
/*****************************************************************************/
void precess_1950_1900 (ra0, dec0, ra1_ptr, dec1_ptr)
    double ra0, dec0, *ra1_ptr, *dec1_ptr;
{
  void precess();
  precess (ra0, dec0, ra1_ptr, dec1_ptr, 1950, 1900);
}

/*****************************************************************************/
/* 	Precess 1900.0 RA & DEC --> 1950.0.				     */
/* 	All coordinates & angles input & output in RADIANS.		     */
/*****************************************************************************/
void precess_1900_1950 (ra0, dec0, ra1_ptr, dec1_ptr)
    double ra0, dec0, *ra1_ptr, *dec1_ptr;
{
  void precess();
  precess (ra0, dec0, ra1_ptr, dec1_ptr, 1900, 1950);
}

/*****************************************************************************/
/** given decimal 1950.0 equatorial coordinates, return decimal Galactic **/
/** all angles in RADIANS **/
/*****************************************************************************/
void eq1950_to_gal (ra1950, dec1950, glon_ptr, glat_ptr)
  double ra1950, dec1950, *glon_ptr, *glat_ptr;
{
  void limit_celestial_coords();
  double smart_atan(), ngp_ra1950, ngp_dec1950, asc_node1950, glon, glat;

  ngp_ra1950   = 192.25 * PI / 180.0;  /* North Galactic Pole RA */
  ngp_dec1950  =  27.4  * PI / 180.0;  /* North Galactic Pole DEC */
  asc_node1950 =  33.0  * PI / 180.0;  /* Asc. Node of Gal. Pl on Eq. (L=33) */

  glat = asin (cos (dec1950) * cos (ngp_dec1950) * cos (ra1950 - ngp_ra1950)
		+ sin (dec1950) * sin (ngp_dec1950));
  glon = smart_atan (  (sin (dec1950) - sin (glat) * sin (ngp_dec1950)), 
		       (cos (dec1950) * sin (ra1950 - ngp_ra1950) 
				      * cos (ngp_dec1950)))
	 + asc_node1950;
  limit_celestial_coords (&glon, &glat);

  *glon_ptr = glon;
  *glat_ptr = glat;
}

/*****************************************************************************/
/** given decimal 2000.0 equatorial coordinates, return decimal Galactic **/
/** all angles in RADIANS **/
/*****************************************************************************/
void eq2000_to_gal (ra2000, dec2000, glon_ptr, glat_ptr)
  double ra2000, dec2000, *glon_ptr, *glat_ptr;
{
  void precess_2000_1950(), eq1950_to_gal();
  double ra1950, dec1950, glon, glat;

  precess_2000_1950 (ra2000, dec2000, &ra1950, &dec1950);
  eq1950_to_gal (ra1950, dec1950, &glon, &glat);

  *glon_ptr = glon;
  *glat_ptr = glat;
}

/*****************************************************************************/
/** given decimal Galactic coordinates, return decimal 1950.0 equatorial **/
/** all angles in RADIANS **/
/*****************************************************************************/
void gal_to_eq1950 (glon, glat, ra1950_ptr, dec1950_ptr)
  double glon, glat, *ra1950_ptr, *dec1950_ptr;
{
  void limit_celestial_coords();
  double smart_atan(), ngp_ra1950, ngp_dec1950, asc_node1950, ra1950, dec1950;

  ngp_ra1950   = 192.25 * PI / 180.0;  /* North Galactic Pole RA */
  ngp_dec1950  =  27.4  * PI / 180.0;  /* North Galactic Pole DEC */
  asc_node1950 =  33.0  * PI / 180.0;  /* Asc. Node of Gal. Pl on Eq. (L=33) */

  dec1950 = asin (  cos (glat) * cos (ngp_dec1950) * sin (glon - asc_node1950)
		  + sin (glat) * sin (ngp_dec1950) );
  ra1950 = smart_atan ( cos (glat) * cos (glon - asc_node1950), 
		        (  sin (glat) * cos (ngp_dec1950) 
			 -   cos (glat) * sin (ngp_dec1950)
			   * sin (glon - asc_node1950)))
	   + ngp_ra1950;
  limit_celestial_coords (&ra1950, &dec1950);

  *ra1950_ptr  = ra1950;
  *dec1950_ptr = dec1950;
}

/*****************************************************************************/
/** given decimal Galactic coordinates, return decimal 2000.0 equatorial **/
/** all angles in RADIANS **/
/*****************************************************************************/
void gal_to_eq2000 (glon, glat, ra2000_ptr, dec2000_ptr)
  double glon, glat, *ra2000_ptr, *dec2000_ptr;
{
  void precess_1950_2000(), gal_to_eq1950(), limit_celestial_coords();
  double ra1950, dec1950, ra2000, dec2000;

  gal_to_eq1950 (glon, glat, &ra1950, &dec1950);
  precess_1950_2000 (ra1950, dec1950, &ra2000, &dec2000);
  limit_celestial_coords (&ra2000, &dec2000);

  *ra2000_ptr  = ra2000;
  *dec2000_ptr = dec2000;
}

/*****************************************************************************/
/** given decimal ecliptic coordinates & equinox, return decimal equatorial **/
/** all angles in RADIANS **/
/*****************************************************************************/
void ecliptic_to_equatorial (lambda, beta, ra_ptr, dec_ptr, equinox)
  double lambda, beta, *ra_ptr, *dec_ptr, equinox;
{
  void limit_celestial_coords();
  double smart_atan(), ecliptic_obliquity(), epsilon, ra, dec;

  epsilon = ecliptic_obliquity (equinox);
  ra  = smart_atan ((sin(lambda)*cos(epsilon) - tan(beta)*sin(epsilon)), 
			cos(lambda));
  dec = asin (sin(beta)*cos(epsilon) + cos(beta)*sin(epsilon)*sin(lambda));
  limit_celestial_coords (&ra, &dec);

  *ra_ptr  = ra;
  *dec_ptr = dec;
}

/*****************************************************************************/
/** given decimal equatorial coordinates & equinox, return decimal ecliptic **/
/** all angles in RADIANS **/
/*****************************************************************************/
void equatorial_to_ecliptic (ra, dec, lambda_ptr, beta_ptr, equinox)
  double ra, dec, *lambda_ptr, *beta_ptr, equinox;
{
  void limit_celestial_coords();
  double smart_atan(), ecliptic_obliquity(), epsilon, lambda, beta;

  epsilon = ecliptic_obliquity (equinox);
  lambda = smart_atan ((sin(ra)*cos(epsilon) + tan(dec)*sin(epsilon)),cos(ra));
  beta   = asin (sin(dec)*cos(epsilon) - cos(dec)*sin(epsilon)*sin(ra));
  limit_celestial_coords (&lambda, &beta);

  *lambda_ptr = lambda;
  *beta_ptr   = beta;
}

/*****************************************************************************/
/** given equinox in years, return obliquity of ecliptic in radians **/
/*****************************************************************************/
double ecliptic_obliquity (equinox)
  double equinox;
{
  double t, epsilon;

  /* Strictly speaking, t is supposed to be the number of Julian centuries
     since 1900 January 0.5.  This should be close enough however. */
  t = (equinox - 1900.0) / 100.0;

  /* epsilon is the mean obliquity of the ecliptic, in radians */
  epsilon = 0.409319755 - 2.27111e-04*t - 2.86e-08*t*t - 8.775e-09*t*t*t;

  return epsilon;
}

/*****************************************************************************/
/** read RA & DEC in sexigesimal format and return radians 		    **/
/*****************************************************************************/
int read_coords (instr, clon_ptr, clat_ptr, FORMTYPE)
  char instr[];
  int FORMTYPE;
  double *clon_ptr, *clat_ptr;
{
  void exit(), h2hms(), hms2h();
  
  char coordstr[1024];
  int i, h_int1, h_int2, m_int1, m_int2, PARSE_OK;
  double s_real1, s_real2, h_real, h_real2;
  
  /* clean up string before parsing */
  strcpy (coordstr, instr);
  for (i=0; i<strlen(coordstr); i++) {
    if (   coordstr[i] == '+' || coordstr[i] == ',' 
	|| coordstr[i] == ';' || coordstr[i] == ':'
	|| coordstr[i] == '(' || coordstr[i] == ')'
	|| coordstr[i] == '[' || coordstr[i] == ']'
	|| coordstr[i] == '{' || coordstr[i] == '}'
	|| coordstr[i] == 'H' || coordstr[i] == 'h'
	|| coordstr[i] == 'D' || coordstr[i] == 'd'
	|| coordstr[i] == 'M' || coordstr[i] == 'm'
	|| coordstr[i] == 'S' || coordstr[i] == 's'
	|| coordstr[i] == '\'' || coordstr[i] == '\"') {
      coordstr[i] = ' ';
    }
  }
  *clon_ptr = 0.0;
  *clat_ptr = 0.0;
  PARSE_OK = FALSE;
  switch (FORMTYPE) {
    case SEX_HRS :
      if (sscanf (coordstr, "%d %d %lf %d %d %lf", 
  	  &h_int1, &m_int1, &s_real1, &h_int2, &m_int2, &s_real2) == 6) {
        if (s_real1 >= 60.0) { s_real1 -= 60.0; m_int1 += 1; }  /** kludge **/
        if (m_int1  >= 60  ) { m_int1  -= 60  ; h_int1 += 1; }
        if (h_int1  >= 24  ) { h_int1  -= 24;                }
        if (s_real2 >= 60.0) { s_real2 -= 60.0; m_int2 += 1; }
        if (m_int2  >= 60  ) { m_int2  -= 60  ; h_int2 += 1; }
        if (h_int1 >= 24 || h_int1 < 0 || m_int1 >= 60 || m_int1 < 0 
  	    || s_real1 >= 60.0 || s_real1 < 0.0) {
  	  fprintf (stderr, "*** Illegal CLON value: |%s| ***\n", coordstr);
  	  exit(1);
        }
        hms2h (h_int1, m_int1, s_real1, &h_real);
        if (h_real >= 24.0 || h_real < 0.0) {
  	  fprintf (stderr, "*** Illegal CLON value: |%s| ***\n", coordstr);
  	  exit(1);
        }
        *clon_ptr = h_real * 15.0 * PI / 180.0;
        if (h_int2 > 90 || h_int2 < -90 || m_int2 >= 60 || m_int2 < 0 
  	    || s_real2 >= 60.0 || s_real2 < 0.0) {
  	  fprintf (stderr, "*** Illegal CLAT value: |%s| ***\n", coordstr);
  	  exit(1);
        }
        hms2h (h_int2, m_int2, s_real2, &h_real);
        if (h_real >= 90.0 || h_real < -90.0) {
  	  fprintf (stderr, "*** Illegal CLAT value: |%s| ***\n", coordstr);
  	  exit(1);
        }
        *clat_ptr = h_real * PI / 180.0;
	PARSE_OK = TRUE;
      }
      break;
    case SEX_DEG : 
      if (sscanf (coordstr, "%d %d %lf %d %d %lf", 
  	  &h_int1, &m_int1, &s_real1, &h_int2, &m_int2, &s_real2) == 6) {
	if (s_real1 >= 60.0) { s_real1 -= 60.0; m_int1 += 1; }  /** kludge **/
	if (m_int1  >= 60  ) { m_int1  -= 60  ; h_int1 += 1; }
	if (h_int1  >= 360 ) { h_int1  -= 360 ;              }
	if (s_real2 >= 60.0) { s_real2 -= 60.0; m_int2 += 1; }
	if (m_int2  >= 60  ) { m_int2  -= 60  ; h_int2 += 1; }
	if (h_int1 >= 360 || h_int1 < 0 || m_int1 >= 60 || m_int1 < 0 
	    || s_real1 >= 60.0 || s_real1 < 0.0) {
	  fprintf (stderr, "*** Illegal CLON value: |%s| ***\n", coordstr);
	  exit(1);
	}
	hms2h (h_int1, m_int1, s_real1, &h_real);
	if (h_real >= 360.0 || h_real < 0.0) {
	  fprintf (stderr, "*** Illegal CLON value: |%s| ***\n", coordstr);
	  exit(1);
	}
	*clon_ptr = h_real * PI / 180.0;
	if (h_int2 > 90 || h_int2 < -90 || m_int2 >= 60 || m_int2 < 0 
	    || s_real2 >= 60.0 || s_real2 < 0.0) {
	  fprintf (stderr, "*** Illegal CLAT value: |%s| ***\n", coordstr);
	  exit(1);
	}
	hms2h (h_int2, m_int2, s_real2, &h_real);
	if (h_real >= 90.0 || h_real < -90.0) {
	  fprintf (stderr, "*** Illegal CLAT value: |%s| ***\n", coordstr);
	  exit(1);
	}
	*clat_ptr = h_real * PI / 180.0;
	PARSE_OK = TRUE;
      }
      break;
    case DEC_HRS : 
      if (sscanf (coordstr, "%le %le", &h_real, &h_real2) == 2) {
	if (h_real >= 24.0 || h_real < 0.0) {
	  fprintf (stderr, "*** Illegal CLON value: |%s| ***\n", coordstr);
	  exit(1);
	}
	*clon_ptr = h_real * 15.0 * PI / 180.0;
	if (h_real2 >= 90.0 || h_real2 < -90.0) {
	  fprintf (stderr, "*** Illegal CLAT value: |%s| ***\n", coordstr);
	  exit(1);
	}
	*clat_ptr = h_real2 * PI / 180.0;
	PARSE_OK = TRUE;
      }
      break;
    case DEC_DEG : 
      if (sscanf (coordstr, "%le %le", &h_real, &h_real2) == 2) {
	if (h_real >= 360.0 || h_real < 0.0) {
	  fprintf (stderr, "*** Illegal CLON value: |%s| ***\n", coordstr);
	  exit(1);
	}
	*clon_ptr = h_real * PI / 180.0;
	if (h_real2 >= 90.0 || h_real2 < -90.0) {
	  fprintf (stderr, "*** Illegal CLAT value: |%s| ***\n", coordstr);
	  exit(1);
	}
	*clat_ptr = h_real2 * PI / 180.0;
	PARSE_OK = TRUE;
      }
      break;
    default :
      fprintf (stderr, "*** error: illegal data format type flag ***\n");
      exit (1);
  }
  return PARSE_OK;
}

/*****************************************************************************/
/** given decimal h, return h m s (int int real) **/
/*****************************************************************************/
void h2hms (h_real, h_int_ptr, m_int_ptr, s_real_ptr)
    double h_real, *s_real_ptr;
    int *h_int_ptr, *m_int_ptr;
{
    int h_int, m_int, sign = 1;
    double s_real, m_real;

    if (h_real < 0.0) h_real *= (sign = -1);
    h_int  = sign * (int) floor (h_real);
    m_real = 60.0 * (h_real - floor (h_real));
    m_int  = (int) floor (m_real);
    s_real = 60.0 * (m_real - floor (m_real));

    *h_int_ptr  = h_int;
    *m_int_ptr  = m_int;
    *s_real_ptr = s_real;
}

/*****************************************************************************/
/*** given h m s (int int real), return decimal h ***/
/*****************************************************************************/
void hms2h (h_int, m_int, s_real, h_real_ptr)
    int h_int, m_int;
    double s_real, *h_real_ptr;
{
    int sign = 1;

    if (h_int < 0) h_int *= (sign = -1);
    *h_real_ptr = sign * (1.0*h_int + (1.0/60.0)*m_int + (1.0/3600.0)*s_real);
}

/*****************************************************************************/
/** given Y and X, compute angle (in radians) with sign ambiguities removed **/
/*****************************************************************************/
double smart_atan (y, x)
  double y, x;
{
  double theta;

  if (x >= 0.0 && y == 0.0) theta = 0.0*PI;
  if (x >  0.0 && y >  0.0) theta = 0.0*PI + atan (y / x);
  if (x == 0.0 && y >  0.0) theta = 0.5*PI;
  if (x <  0.0 && y >  0.0) theta = 1.0*PI + atan (y / x);
  if (x <  0.0 && y == 0.0) theta = 1.0*PI;
  if (x <  0.0 && y <  0.0) theta = 1.0*PI + atan (y / x);
  if (x == 0.0 && y <  0.0) theta = 1.5*PI;
  if (x >  0.0 && y <  0.0) theta = 2.0*PI + atan (y / x);

  return theta;
}

/*****************************************************************************/
/** given a celestial longitude & latitude pair (RA+DEC, L+B, etc.),	    **/
/** ensure that 0 <= CLON < 2*PI and -PI/2 <= CLAT <= +PI/2		    **/
/** all angles in RADIANS 						    **/
/*****************************************************************************/
void limit_celestial_coords (clon_ptr, clat_ptr)
  double *clon_ptr, *clat_ptr;
{
  while (*clon_ptr >= 2.0*PI) *clon_ptr -= 2.0*PI;
  while (*clon_ptr <  0.0   ) *clon_ptr += 2.0*PI;

  if (*clat_ptr > 0.5*PI || *clat_ptr < -0.5*PI) {
    *clat_ptr = asin (sin (*clat_ptr));
  }
}

/*****************************************************************************/
/** given a celestial longitude & latitude pair (RA+DEC, L+B, etc.),	    **/
/** ensure that -PI < CLON <= +PI and -PI/2 <= CLAT <= +PI/2		    **/
/** all angles in RADIANS 						    **/
/*****************************************************************************/
void limit_celestial_coords_2 (clon_ptr, clat_ptr)
  double *clon_ptr, *clat_ptr;
{
  while (*clon_ptr >     PI  ) *clon_ptr -= 2.0*PI;
  while (*clon_ptr <= -1.0*PI) *clon_ptr += 2.0*PI;

  if (*clat_ptr > 0.5*PI || *clat_ptr < -0.5*PI) {
    *clat_ptr = asin (sin (*clat_ptr));
  }
}

/*****************************************************************************/
/** force RA value to lie in the range [0.0,2.0*PI) **/
/*****************************************************************************/
void limit_ra (ra_ptr)
    double *ra_ptr;
{
    while (*ra_ptr <     0.0) *ra_ptr += 2.0*PI;
    while (*ra_ptr >= 2.0*PI) *ra_ptr -= 2.0*PI;
}

	/* Mathematical symbols and constants */

#define FALSE		0
#define TRUE		1

#ifndef PI
#define PI		3.14159265358979323846264
#endif

	/* Physical constants [cgs] */

#define	c_light		2.99792458e10	/* cm sec^-1	 */
#define	h_Planck	6.6260755e-27	/* erg sec 	 */
#define k_Boltzmann	1.380658e-16	/* erg K^-1 	 */
#define G_grav		6.67259e-08	/* dyn cm^2 g^-2 */
					/* 5.67051e-05 erg cm^-2 K^-4 s^-1 */
#define sigma_SB	(2.0 * PI*PI*PI*PI*PI * k_Boltzmann*k_Boltzmann*k_Boltzmann*k_Boltzmann / (15 * h_Planck*h_Planck*h_Planck * c_light*c_light))

	/* Angular unit conversions */

#define RAD2D		(       180.0/PI)	/* radians -> degrees */
#define RAD2AM		(  60.0*180.0/PI)	/* radians -> arcmin  */
#define RAD2AS		(3600.0*180.0/PI)	/* radians -> arcsec  */
#define RAD2H		(        12.0/PI)	/* radians -> hours   */
#define RAD2M		(  60.0* 12.0/PI)	/* radians -> minutes */
#define RAD2S		(3600.0* 12.0/PI)	/* radians -> seconds */

#define  D2RAD		(1.0 / RAD2D)      	/* degrees -> radians */
#define AM2RAD		(1.0 / RAD2AM)   	/* arcmin  -> radians */
#define AS2RAD		(1.0 / RAD2AS)     	/* arcsec  -> radians */
#define  H2RAD		(1.0 / RAD2H)      	/* hours   -> radians */
#define  M2RAD		(1.0 / RAD2M)     	/* minutes -> radians */
#define  S2RAD		(1.0 / RAD2S)      	/* seconds -> radians */

	/* Other unit conversions */

#define SIGMA_TO_FWHM	2.354820044		/* (2.0*sqrt(2.0*ln(2.0)) */
#define FWHM_TO_SIGMA	(1.0 / SIGMA_TO_FWHM)

#define MAG_TO_DEPTH	0.9210340372		/* 0.4 * ln(2.0) */
#define DEPTH_TO_MAG	(1.0 / MAG_TO_DEPTH)

#define PC2CM	3.085678e+18	/* number of centimeters in one parsec */
#define CM2PC	(1.0/PC2CM)	/* number of parsecs in one centimeter */

#define HMASS	1.6735e-24	/* mass of Hydrogen atom in grams */
#define SMASS	1.9891e+33	/* mass of Sun in grams */
#define R_Sun	6.9599e+10	/* radius of Sun in cm */

	/* HI column density conversion factors; C_AREA assumes Gaussian */

#define C_AREA  5.2e-19		/* N_HI = (tau_max * fwhm * T_s) / C_AREA */
#define HI_COL  1.823e+18	/* N_HI(tau<<1) = HI_COL * Integral{Tb * dv} */
#define H2_COL	2.3e+20		/* N_H2 = H2_COL * Integral{Tb(12CO) * dv} */
				   /* from Strong et al. (1988 A&A 207, 1) */

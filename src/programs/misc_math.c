/*############################################################################
 #############################################################################
 ##									    ##
 ##			 M I S C _ M A T H . C				    ##
 ##									    ##
 ##	Source code for sundry miscellaneous math utility routines, to be   ##
 ##	#included in various data reduction codes. 			    ##
 ##									    ##
 #############################################################################
 ############################################################################*/

#include "mathdefs.h"
#include "misc_math.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/*****************************************************************************
 **									    **
 **		               F S Q R	        			    **
 **									    **
 ** 		Return square of a double-precision float. 		    **
 **									    **
 *****************************************************************************/
double fsqr (x)
  double x;
{
  return x * x;
}

/*****************************************************************************
 **									    **
 **		               R O U N D        			    **
 **									    **
 ** 		Round up if frac >= 0.5, else down 			    **
 **									    **
 *****************************************************************************/
double round (x)
    double x;
{
    double roundval;

    if ((x - floor(x)) < 0.5) {
        roundval = floor (x);
    } else {
        roundval = ceil (x);
    }
    return roundval;
}

/*****************************************************************************
 **									    **
 **		               L E S S E R        			    **
 **									    **
 ** 		Returns the lesser of two integers.   			    **
 **									    **
 *****************************************************************************/
int lesser (a, b)
    int a, b;
{
    if (a < b) {
	return a;
    } else {
	return b;
    }
}

/*****************************************************************************
 **									    **
 **		               G R E A T E R        			    **
 **									    **
 ** 		Returns the greater of two integers.   			    **
 **									    **
 *****************************************************************************/
int greater (a, b)
    int a, b;
{
    if (a > b) {
	return a;
    } else {
	return b;
    }
}

/*****************************************************************************
 **									    **
 **		             F L E S S E R        			    **
 **									    **
 ** 		Returns the lesser of two doubles.   			    **
 **									    **
 *****************************************************************************/
double flesser (a, b)
    double a, b;
{
    if (a < b) {
	return a;
    } else {
	return b;
    }
}

/*****************************************************************************
 **									    **
 **		              F G R E A T E R        			    **
 **									    **
 ** 		Returns the greater of two doubles.   			    **
 **									    **
 *****************************************************************************/
double fgreater (a, b)
    double a, b;
{
    if (a > b) {
	return a;
    } else {
	return b;
    }
}

/*****************************************************************************
 **									    **
 **		              I N T E R P            			    **
 **									    **
 ** 	Interpolates linearly the Y value of a point between two others.    **
 **									    **
 *****************************************************************************/
double interp (x, x1, x2, y1, y2)
  double x1, x2, y1, y2, x;
{
  double f, return_val;

  if ((f=(x2-x1)) != 0.0) {
    return_val = (y1 + (y2 - y1) * ((x - x1) / f));
  } else {
    fprintf (stderr, "*** ERROR in INTERP: x=%e x1=%e x2=%e y1=%e y2=%e ***\n",
							x, x1, x2, y1, y2);
    exit(1);
  }
  return return_val;
}

/*****************************************************************************
 **									    **
 **		           L O G _ I N T E R P            		    **
 **									    **
 ** Interpolates LOGARITHMICALLY the Y value of a point between two others. **
 **									    **
 *****************************************************************************/
double log_interp (x, x1, x2, y1, y2)
  double x, x1, x2, y1, y2;
{
  double f, return_val;

  if (x1!=0.0 && x2!=0.0 && y1!=0.0 && y2!=0.0 && (f=log(x2/x1))!=0.0) {
    return_val = (y1 * pow (x/x1, log(y2/y1) / f));
  } else {
    fprintf (stderr, "*** ERROR in LOG_INTERP ***\n");
    exit(1);
  }
  return return_val;
}

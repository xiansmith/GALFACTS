/*############################################################################
 #############################################################################
 ##									    ##
 ##			       I O . C					    ##
 ##									    ##
 ##	Source code for a few primitive I/O routines, to be	 	    ##
 ##	#included where needed, e.g., in FITSIO.C. 			    ##
 ##									    ##
 #############################################################################
 ############################################################################*/

#include "io.h"
#include <stdlib.h>


/*****************************************************************************
 **									    **
 **			     O P E N _ R E A D  			    **
 **									    **
 **	Open file of given name for read operations; die if error.	    **
 **									    **
 *****************************************************************************/
void open_read (fileptr_ptr, fname)
  const char fname[];
  FILE **fileptr_ptr;
{

  if ((*fileptr_ptr = fopen (fname, "r")) == NULL) {
    fprintf (stderr, "*** error: cannot open %s for reading ***\n", fname);
    exit(1);
  }
}

/*****************************************************************************
 **									    **
 **			     O P E N _ W R I T E			    **
 **									    **
 **	Open file of given name for write operations; die if error.	    **
 **									    **
 *****************************************************************************/
void open_write (fileptr_ptr, fname)
  const char fname[];
  FILE **fileptr_ptr;
{

  if ((*fileptr_ptr = fopen (fname, "w")) == NULL) {
    fprintf (stderr, "*** error: cannot open %s for writing ***\n", fname);
    exit(1);
  }
}


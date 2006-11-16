/* ARITH_FITS.C -- perform 2-operand (+ - * /) arithmetic on specified FITS
   images and output result; if a filename is readable as a number, treat is
   as such.

	COMMAND LINE ARGUMENTS: 1 input image file name
				2 operand (one of A,+,-,*,/,N,X,V)
				3 second image file name or numeric constant
				4 output image file name
				5 output image object string (optional)
				6 output image bunit string (optional)
				7 zeroes for blanks? (opt; 0=no=default; 1=yes)

   adapted  2002 July 22 from sub_fits.c */

#define _LARGEFILE64_SOURCE	/* Allow output > 2 GB; see io.c for details.*/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <errno.h>
#include "/u/gibson/src/util/chardefs.h"
#include "/u/gibson/src/util/mathdefs.h"
#include "/u/gibson/src/util/misc.c"
#include "/u/gibson/src/util/misc_math.c"
#include "/u/gibson/src/fits/fitsio.c"
#include "/u/gibson/src/util/io_big.c"

		/* operand flag values */

#define ABS_OP  0	/* output = |op1| (op2 ignored) */
#define ADD_OP	1	/* output = op1 + op2 */
#define SUB_OP	2	/* output = op1 - op2 */
#define MUL_OP	3	/* output = op1 * op2 */
#define DIV_OP	4	/* output = op1 / op2 */
#define MIN_OP  5	/* output = minimum of op1 and op2 */
#define MAX_OP  6	/* output = maximum of op1 and op2 */
#define AVG_OP  7	/* output = 0.5 * (op1 + op2) */

/*****************************************************************************
 **									    **
 **				M A I N					    **
 **									    **
 *****************************************************************************/
main (argc, argv)
  int argc;		/* # command line args, INCLUDING command itself */
  char *argv[];		/* arg strings; argv[0] is command string */
{
  void exit();

  char opchar, inname1[256], inname2[256], outname[256];
  int i, j, k, kk, n_planes, ndim, m[MAX_NUM_AXES], 
	NUM_OP, UNARY_OP, OPFLAG, KLUDGE_MINMAX, ZERO_BLANKS, CUBE_AND_PLANE;
  size_t n_total, n[MAX_NUM_AXES];
  float *inplane1, *inplane2, *outplane;
  double op1, op2, outval, kludge_min, kludge_max;
  header_param_list hpar_in1, hpar_in2, hpar_out;
  FILE *infile1, *infile2, *outfile;

  KLUDGE_MINMAX = TRUE;
  if (argc-1 != 4 && argc-1 != 5 && argc-1 != 6 && argc-1 != 7) {
    fprintf (stderr, "calling syntax:\n");
    fprintf (stderr, 
	     "  arith_fits <infl1> <op(A+-*/NXV)> <infl2|const> <outfl>\n");
    fprintf (stderr, "      [<obj>] [<bunit>] [<zero_blanks(1=yes)>]\n");
    exit (0);
  }
  strcpy (inname1, argv[1]);
  sscanf (argv[2], "%c", &opchar);
  strcpy (inname2, argv[3]);
  strcpy (outname, argv[4]);

  if (sscanf (inname2, "%le", &op2) == 1) {
    NUM_OP = TRUE;
  } else {
    NUM_OP = FALSE;
  }
  switch (opchar) {
    case 'A' : OPFLAG = ABS_OP;  UNARY_OP = TRUE;   break;
    case '+' : OPFLAG = ADD_OP;  UNARY_OP = FALSE;  break;
    case '-' : OPFLAG = SUB_OP;  UNARY_OP = FALSE;  break;
    case '*' : OPFLAG = MUL_OP;  UNARY_OP = FALSE;  break;
    case '/' : OPFLAG = DIV_OP;  UNARY_OP = FALSE;  break;
    case 'N' : OPFLAG = MIN_OP;  UNARY_OP = FALSE;  break;
    case 'X' : OPFLAG = MAX_OP;  UNARY_OP = FALSE;  break;
    case 'V' : OPFLAG = AVG_OP;  UNARY_OP = FALSE;  break;
    default : fprintf (stderr, "*** illegal operator %c ***\n",opchar);exit(1);
  }
  if (NUM_OP == TRUE && op2 == 0.0 && OPFLAG == DIV_OP) {
    fprintf (stderr, "*** specified operation would divide image by 0! ***\n");
    exit (1);
  }
  ZERO_BLANKS = FALSE;
  if (argc-1 >= 7) {
    sscanf (argv[7], "%d", &ZERO_BLANKS);
  }
  switch (ZERO_BLANKS) {
    case 0 : printf (" BLANK pixels not computed.\n"); 
	     ZERO_BLANKS = FALSE; break;
    case 1 : printf (" BLANK pixels treated as zeros in computation.\n"); 
	     ZERO_BLANKS = TRUE; break;
    default : fprintf (stderr, "*** illegal value %d for ZERO_BLANKS! ***\n");
	      exit (1);
  }
  open_read (&infile1, inname1);
  readfits_header (infile1, &hpar_in1);
  ndim = hpar_in1.num_axes;
  if (ndim < 2 || ndim > MAX_NUM_AXES) {
    fprintf (stderr, 
	"*** error: number of dimensions must be between 2 and %d ***\n",
							MAX_NUM_AXES, argv[1]);
    exit(1);
  }
  printf (" Image 1 number of axes = %d\n", ndim);
  n_total = 1;
  for (i=0; i<ndim; i++) {
    n[i] = hpar_in1.naxis[i];
    n_total *= n[i];
    printf ("  naxis%d=%d", i+1, n[i]);
  }
  printf ("\n");
  n_planes = (int) (n_total / (n[0] * n[1]));
  printf (" n_total = %ld;  n_planes = %d\n", n_total, n_planes);
  if (NUM_OP == FALSE && UNARY_OP == FALSE) {
    open_read (&infile2, inname2);
    readfits_header (infile2, &hpar_in2);
    printf (" Image 2 number of axes = %d\n", hpar_in2.num_axes);
    for (i=0; i<hpar_in2.num_axes; i++) {
      printf ("  naxis%d=%d", i+1, hpar_in2.naxis[i]);
    }
    printf ("\n");
    CUBE_AND_PLANE = FALSE;
    for (i=0; i<ndim; i++) {
	/** Check to see if any axis lengths disagree.  Note this is a more
	    forgiving test than checking whether the numbers of axes disagree,
	    since it ignores degenerate axes.  However, the only reason this
	    works is because READFITS_HEADER calls INIT_HEADER_PARAM_LIST,
	    which initializes all axis lengths to 1, including those whose rank
	    exceeds NUM_AXES.  If the initialization is changed, this code will
	    not perform properly. 
	    One exception is allowed: if the first image is a cube and the
	    second is a plane. **/	    
      if (n[i] != hpar_in2.naxis[i]) {
	if (i == 2 && n[i] > 1 && hpar_in2.naxis[i] == 1) {
	  CUBE_AND_PLANE = TRUE;
	} else {
	  fprintf (stderr,"*** error: input image dimensions disagree! ***\n");
	  exit (1);
	}
      }
    }
    if (CUBE_AND_PLANE == TRUE) {
      /* This presumes that no other axis disagreements were found above,
	 since these would have caused the code to exit. */
      fprintf (stderr, "%%% Warning: n3_1 = %d and n3_2 = %d! %%%\n",
						n[2], hpar_in2.naxis[2]);
      fprintf (stderr, "%%% Will perform (cube) (op) (plane) = (cube). %%%\n");
    }
  }
  copy_header_param_list (&hpar_out, &hpar_in1);
  if (argc-1 >= 5) {
    strcpy (hpar_out.object, argv[5]);
  }
  if (argc-1 == 6) {
    strcpy (hpar_out.bunit, argv[6]);
  }
  if (hpar_out.bitpix==16) {
    if (KLUDGE_MINMAX == TRUE) {
      kludge_min = -300.0;
      kludge_max =  300.0;
      fprintf (stderr, 
	      "%%% Warning: kludge min & max 16-bit scaling in effect! %%%\n");
      fprintf (stderr, "    Using kludge min = %f, max = %f\n", 
						kludge_min, kludge_max);
      set_bzero_bscale (&hpar_out, kludge_min, kludge_max);
    } else {
      fprintf (stderr,
		"*** error: no scaling protocol for 16-bit output! ***\n");
      exit (1);
    }
  }
  open_write_big (&outfile, outname);
  writefits_header (outfile, &hpar_out);
  inplane1 = (float *) malloc (n[0] * n[1] * sizeof (float));
  if (NUM_OP == FALSE && UNARY_OP == FALSE) {
    inplane2 = (float *) malloc (n[0] * n[1] * sizeof (float));
    if (CUBE_AND_PLANE == TRUE) {
      readfits_plane (infile2, inplane2, &hpar_in2);
    }
  }
  outplane = (float *) malloc (n[0] * n[1] * sizeof (float));
  if (UNARY_OP == FALSE) {
    printf ("Computing (%s) %c (%s) = (%s).\n",inname1,opchar,inname2,outname);
  } else {
    printf ("Computing %c (%s) = (%s).\n",opchar,inname1,outname);
  }
  for (k=0; k<n_planes; k++) {
    kk = k;
    for (i=ndim-1; i>=3; i--) {
      m[i] = kk / ((int) n[i-1]);
      kk -= m[i] * ((int) n[i-1]);
    }
    m[2] = kk;
    printf ("  Reading plane %3d of %3d; by dimension:", k+1, n_planes);
    for (i=2; i<ndim; i++) {
      printf ("  %3d of %3d", m[i]+1, n[i]);
      if (ndim-2 > 1 && i < ndim-1) printf (",");
    }
    printf ("\n%c[1A", 27);
    readfits_plane (infile1, inplane1, &hpar_in1);
    if (NUM_OP == FALSE && CUBE_AND_PLANE == FALSE && UNARY_OP == FALSE) {
      readfits_plane (infile2, inplane2, &hpar_in2);
    }
    for (j=0; j<n[1]; j++) {
      for (i=0; i<n[0]; i++) {
	outval = BLANK_PIXEL;
	op1 = inplane1[i+j*n[0]];
	if (NUM_OP == FALSE && UNARY_OP == FALSE) op2 = inplane2[i+j*n[0]];
	if (ZERO_BLANKS == TRUE) {
	  if (op1 >= BLANK_THRESH) op1 = 0.0;
	  if (op2 >= BLANK_THRESH && NUM_OP == FALSE) op2 = 0.0;
	}
	if (   op1 < BLANK_THRESH 
	    && (UNARY_OP == TRUE || NUM_OP == TRUE || op2 < BLANK_THRESH)) {
	  switch (OPFLAG) {
	    case ABS_OP : outval = fabs (op1);  break;
	    case ADD_OP : outval = op1 + op2;  break;
	    case SUB_OP : outval = op1 - op2;  break;
	    case MUL_OP : outval = op1 * op2;  break;
	    case DIV_OP : if (op2 != 0.0) outval = op1 / op2;  break;
	    case MIN_OP : if (op1<op2) outval=op1; else outval=op2; break;
	    case MAX_OP : if (op1>op2) outval=op1; else outval=op2; break;
	    case AVG_OP : outval = 0.5 * (op1 + op2);
	  }
	}
	outplane[i+j*n[0]] = (float) outval;
      }
    }
    writefits_plane (outfile, outplane, &hpar_out);
  }
  writefits_pad_end (outfile, &hpar_out);
  fclose (outfile);
  fclose (infile1);
  if (NUM_OP == FALSE && UNARY_OP == FALSE) fclose (infile2);
  free (outplane);
  free (inplane1);
  if (NUM_OP == FALSE && UNARY_OP == FALSE) free (inplane2);
  printf ("\nDone.\n");
}

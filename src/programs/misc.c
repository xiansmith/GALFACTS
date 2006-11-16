/*############################################################################
 #############################################################################
 ##									    ##
 ##				M I S C . C				    ##
 ##									    ##
 ##	Source code for sundry miscellaneous utility routines, to be 	    ##
 ##	#included in various data reduction codes. 			    ##
 ##									    ##
 ##	consolidated Monday 29 May 1995 sjg				    ##
 ##									    ##
 #############################################################################
 ############################################################################*/

/** The following definition files must be included prior to MISC.C: 

	chardefs.h

**/

#define MAX_PARSE_LENGTH 	1024	/* maximum line length for parsing */

/*****************************************************************************
 **									    **
 **			     W A I T _ F O R _ C R			    **
 **									    **
 **	Print message and wait for carriage return from standard input.	    **
 **									    **
 *****************************************************************************/
void wait_for_CR (message)
    char message[];
{
    char c;

    c = NULL_CHAR;
    fprintf (stderr, "%s", message);
    do {
	c=fgetc(stdin);
/*
    } while (!feof(stdin) && c!=LINE_FEED && c!=END_OF_FILE);
*/
    } while (!feof(stdin) && c!=LINE_FEED);
/*
    printf ("last c=%d\n", c);
*/
}

/*****************************************************************************
 **									    **
 **			     I S I N S T R				    **
 **									    **
 **		Return TRUE if second string is in first. 		    **
 **									    **
 *****************************************************************************/
int isinstr (s, sub)
  char s[], sub[];
{
  int i, j, slen, sublen, return_val;

  return_val = FALSE;
  if ((slen = strlen (s)) >= (sublen = strlen (sub))) {
    for (i=0; i<slen; i++) {
/*
      if (strncmp(&s[i], sub, sublen) == 0) return_val = TRUE;
*/
      j = i;
      for (j=i; j<i+sublen; j++) {
	if (j==slen || sub[j-i] != s[j]) break;
	if (j==i+sublen-1) return_val = TRUE;
      }
    }
  }
  return return_val;
}

/*****************************************************************************
 **									    **
 **			     S U B S T R P O S				    **
 **									    **
 **     Return 1st-char position of second string inside first string. 	    **
 **	If 1st does not contain second, return -1.			    **
 **									    **
 *****************************************************************************/
int substrpos (s, sub)
  char s[], sub[];
{
  int i, j, slen, sublen, return_val;

  return_val = -1;
  if ((slen = strlen (s)) >= (sublen = strlen (sub))) {
    i = 0;
    while (i < slen && return_val == -1) {
/*
      if (strncmp(&s[i], sub, sublen) == 0) return_val = i;
*/
      j = i;
      for (j=i; j<i+sublen; j++) {
	if (j==slen || sub[j-i] != s[j]) break;
	if (j==i+sublen-1) return_val = i;
      }
      i++;
    }
  }
  return return_val;
}

/*****************************************************************************
 **									    **
 **			     C O M M A _ I N T    			    **
 **									    **
 **     Write integer output to a string, formatted with commas every	    **
 **	three digits.                                       		    **
 **									    **
 *****************************************************************************/
void comma_int (com_str, num)
  char com_str[];
  int num;
{
  char uncom_str[MAX_PARSE_LENGTH+1];
  int i, j, k, comma_flag[MAX_PARSE_LENGTH+1];

  sprintf (uncom_str, "%d", num);
  i = strlen (uncom_str) - 1;
  j = 0;
  while (uncom_str[i] != ' ' && i >= 0) {
    if (j > 0 && j % 3 == 0) {
      comma_flag[i] = TRUE;
    } else {
      comma_flag[i] = FALSE;
    }
    j++;
    i--;
  }
  i = strlen (uncom_str);
  k = 0;
  for (j=0; j<i; j++) {
    com_str[k++] = uncom_str[j];
    if (comma_flag[j] == TRUE) com_str[k++] = ',';
  }
  com_str[k++] = END_OF_STRING;
}

/*****************************************************************************
 **									    **
 **			     S G E T _ I N T      			    **
 **									    **
 **     Return value of int in string.  If sscanf gives an error,   	    **
 **	print an informative error message and exit program.		    **
 **									    **
 *****************************************************************************/
int sget_int (buffer, varname)
  char buffer[], varname[];
{
  void exit();
  int return_val;

  if (sscanf (buffer, "%d", &return_val) != 1) {
    fprintf (stderr, "*** error: improper value for %s ***\n", varname);
    exit (1);
  }
  return return_val;
}

/*****************************************************************************
 **									    **
 **			     S G E T _ D O U B L E			    **
 **									    **
 **     Return value of double in string.  If sscanf gives an error,   	    **
 **	print an informative error message and exit program.		    **
 **									    **
 *****************************************************************************/
double sget_double (buffer, varname)
  char buffer[], varname[];
{
  void exit();
  double return_val;

  if (sscanf (buffer, "%le", &return_val) != 1) {
    fprintf (stderr, "*** error: improper value for %s ***\n", varname);
    exit (1);
  }
  return return_val;
}

/*****************************************************************************
 **									    **
 **		              G E T _ L I N E                               **
 **									    **
 **	Read all characters from file pointer into input buffer until	    **
 **	either a line feed or end of file is encountered; the line feed	    **
 **	itself is read and discarded, so that a subsequent call will not    **
 **	find it still there.  If either an EOF is encountered or the        **
 **	buffer size is exceeded, the function will return a 0 value;        **
 **	otherwise, the operation is deemed successful, and a 1 is returned. **
 **									    **
 *****************************************************************************/
int get_line (file_ptr, buffer, maxchar)
  FILE *file_ptr;
  char buffer[];
  int maxchar;
{
  char c;
  int i, return_val;
  for (i=0; i<maxchar; i++) buffer[i] = END_OF_STRING;
  i = 0;
  c = 'a';
/*
  while (!feof(file_ptr) && c!=END_OF_FILE && c!=LINE_FEED && i<maxchar) {
    c = fgetc (file_ptr);
    if (c!=END_OF_FILE && c!=LINE_FEED) {
      buffer[i] = c;
      i++;
    }
  }
  buffer[i] = END_OF_STRING;
  if (feof(file_ptr) || c==END_OF_FILE || i==maxchar) return_val = 0;
*/
  while (!feof(file_ptr) && c!=LINE_FEED && i<maxchar) {
    c = fgetc (file_ptr);
    if (c!=LINE_FEED) {
      buffer[i] = c;
      i++;
    }
  }
  buffer[i] = END_OF_STRING;
  if (feof(file_ptr) || i==maxchar) return_val = 0;
  else return_val = 1;
  return return_val;
}

/*****************************************************************************
 **									    **
 **		         T R I M _ W H I T E S P A C E                      **
 **									    **
 **	Trim a string of leading and trailing whitespace characters.	    **
 **									    **
 *****************************************************************************/
void trim_whitespace (str) 
  char str[];
{
  int iswhite();
  int i, j, len, left, right;

  len = strlen (str);
  i = 0; while (str[i] != END_OF_STRING) i++;
  len = i;
  left = 0; while (iswhite(str[left])==TRUE && left < len) left++;
  for (i=0; i<len-left; i++) str[i] = str[i+left];
  for (j=i; j<len; j++) str[j] = SPACE;
  right = len-1; while (iswhite(str[right])==TRUE) right--;
  str[right+1] = END_OF_STRING;
}

/*****************************************************************************
 **									    **
 **		         K I L L _ W H I T E S P A C E                      **
 **									    **
 **	Remove all whitespace characters in a string by compression.	    **
 **									    **
 *****************************************************************************/
void kill_whitespace (str) 
  char str[];
{
  int iswhite();
  int i, j;

  i = 0; 
  while (str[i] != END_OF_STRING) {
    if (iswhite (str[i]) == TRUE) {
      j = i;
      while (str[j] != END_OF_STRING) {
	str[j] = str[j+1];
	j++;
      }
    }
    i++;
  }
}

/*****************************************************************************
 **									    **
 **		             I S W H I T E                                  **
 **									    **
 **	Return TRUE if a character is a whitespace character, else FALSE.   **
 **									    **
 *****************************************************************************/
int iswhite (c)
  char c;
{
  int return_value;

/*
  if (c==NULL_CHAR || c==END_OF_STRING || c==TAB || c==LINE_FEED 
	|| c==VERTICAL_TAB || c==FORM_FEED || c==CARRIAGE_RETURN
	|| c==SPACE || c==END_OF_FILE) {
*/
  if (c==NULL_CHAR || c==END_OF_STRING || c==TAB || c==LINE_FEED 
	|| c==VERTICAL_TAB || c==FORM_FEED || c==CARRIAGE_RETURN || c==SPACE) {
    return_value = TRUE;
  } else {
    return_value = FALSE;
  }
  return return_value;
}

/*****************************************************************************
 **									    **
 **		S M A L L O C -- "smart malloc" - check for malloc error    **
 **									    **
 *****************************************************************************/
char *smalloc (size)
  unsigned size;
{
  void exit();
  char *malloc_return;

  malloc_return = (char *) malloc (size);
  if (malloc_return == NULL) {
    fprintf (stderr, "*** error: malloc call unsuccessful! ***\n");
    exit(1);
  } 
  return malloc_return;
}

/*****************************************************************************
 **									    **
 **		E R R O R  --  print an error message and die		    **
 **									    **
 *****************************************************************************/
void error (message)
    char message[];
{
    void exit();

    fprintf (stderr, "%s\n", message);
    exit(1);
}

/*****************************************************************************/
void reverse_video ()
{
  printf ("%c[7m", 27);
}

/*****************************************************************************/
void normal_video ()
{
  printf ("%c[0m", 27);
}

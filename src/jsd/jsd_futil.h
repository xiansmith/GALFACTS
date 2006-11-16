/* 
 * Jeff's C file utilities.
 * 
 * Jeff Dever (jsd)
 */

#ifndef _JSD_FUTIL_H
#define _JSD_FUTIL_H

#include <stdio.h>
#include <stdlib.h>

/* 
 * Counts the number of lines in the file from the current position.
 * The current position is not changed.
 * 
 * @param file must be open for reading
 * @return -1 if file is null, else the number of lines
*/
int jsd_line_count(FILE * file);

#endif //_JSD_FUTIL_H


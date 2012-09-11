#include <stdio.h>

#ifndef _IO_H
#define _IO_H

void open_read (FILE **fileptr_ptr, const char *fname);
void open_write (FILE **fileptr_ptr, const char *fname);

#endif


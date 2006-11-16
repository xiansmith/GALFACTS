#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chardefs.h"
#include "mathdefs.h"
#include "misc.c"
#include "misc_math.c"
#include "fitsio.c"
//----------------------------------------------------------------------

#define TRUE      1                      // Define some handy constants  
#define FALSE     0                      // Define some handy constants  

int main(void)
{

  int  i, j, k, n1, n2;
  float gamma, omega, theta, sigma;

  for(j=0; j<100; j++) {
   theta = 0.5 + (float)j/20.0;

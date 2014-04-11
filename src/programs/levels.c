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
 
int main(void) {

  int i; 
  float TCMB = 2.73;                   // CMB temperature
  float alpha = -0.7;                  // flux spectral index of synchrotron emission
  float TSYN = 25.0;                   // temperature of synchrotron brightness at 408 MHz
  float TFF = 2.0;                     // free-free brightness temperature at 408 MHz
  float TBCMB, TBSYN, TBFF;
  float startnu = 1.0e08;
  float step = 1.1;
  FILE *outfile;
  float B, nu, x;

  if( (outfile = fopen("levels.out","w") ) == NULL )
     { printf("can't open output file \n\r"); }

  nu = startnu;  
  for (i=0;i<260;i++) {
    nu = nu+startnu*step*(float)i;
    x = h_Planck*nu/(k_Boltzmann*TCMB);
    B = (2.0*h_Planck*nu*nu*nu)/(c_light*c_light)*(1.0/(exp(x)-1.0));
    TBCMB = B * (c_light*c_light)/(2.0*k_Boltzmann*nu*nu);
    TBSYN = TSYN * pow((nu/408.0e06),alpha-2.0);
    TBFF = TFF * pow((nu/408.0e06),-2.1);

  printf("%e %e %e %e %e \n", nu, B, TBCMB, TBSYN, TBFF);
    
  }
}
  


 


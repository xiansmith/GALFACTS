#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/chardefs.h"
#include "programs/mathdefs.h"
#include "programs/misc.c"
#include "programs/misc_math.c"
#include "programs/coord_utils.c"
#include "programs/fitsio.c"
//----------------------------------------------------------------------
//  Program to calibrate the Stokes parameters from the WAPP data
//

#define TRUE      1                      // Define some handy constants  
#define FALSE     0                      // Define some handy constants 

int main(void)
{

  void exit(), readfits_map(), writefits_map();
  int  i, j, k, kt, n, npoints, numscans;
  int RFI = TRUE;        // Set to TRUE to apply RFI flagging
  double RA[60000], DEC[60000], AST[60000];
  int rah, ram, dd, dm;
  int lowchan = 40;
  int highchan = 220;
  int smooth = 25;         //smoothing width in frequency domain for cal signal
  int tsmooth = 5;         //smoothing width in the time domain for sky signal
  float difflimit = 0.1;   // limit on difference signal to flag as RFI
  float diffxx, diffyy;
  float fcen = 1170.0;
  float df = 0.390625;
  float freq[256];
  int flag[256];
  float numchan, count, calsum, goodchan;
  int firstime;
  float ras, ds;
  float Tcal = 1.0;
  float gx, gy;
  float I[60000],Q[60000],U[60000],V[60000];
  float dumI, dumQ, dumU, dumV;
  double dumRA, dumDEC, dumt;
  char epoch[8], dummy[480];
  char line[80];
  FILE  *datafile, *calfile, *skyfile, *fluxfile, *rfifile, *gainfile;

  struct gainset {
    float x[256];
    float y[256];
    float phi[256];
  };
  struct PolSet {
    float xx[256];
    float yy[256];
    float xy[256];
    float yx[256];
  };
  struct StokeSet {
    float I[256];
    float Q[256];
    float U[256];
    float V[256];
 };
 
float calxx, calyy, calxy, calyx, calU, calV;
struct PolSet calon, caloff, cal, corcal;
struct StokeSet ObsCal, CalCal, ObsSky, TrueSky;
struct gainset gain;

  numchan = highchan - lowchan + 1;
  k = (smooth - 1) / 2;
  kt = (tsmooth - 1) /2;
  printf("\n %7.2f %4d %4d %4d %4d \n", numchan, smooth, k, tsmooth, kt);

  for(i=0;i<=255;i++) {
    freq[i] = fcen + ((float)(i - 127))*df;
  }
    

 if( (datafile = fopen("A1863.za_scan.wapp1.53108.spec","rb") ) == NULL )
     { printf("can't open input file \n\r"); }

 if( (calfile = fopen("calspec.dat","w") ) == NULL )
     { printf("can't open input file \n\r"); }

  if( (skyfile = fopen("skyspec.dat","w") ) == NULL )
     { printf("can't open output file \n\r"); }

  if( (fluxfile = fopen("fluxtime.dat","w") ) == NULL )
     { printf("can't open output file \n\r"); }

  if( (gainfile = fopen("gain.dat","w") ) == NULL )
     { printf("can't open output file \n\r"); }

  if( (rfifile = fopen("rfi.dat","w") ) == NULL )
     { printf("can't open output file \n\r"); }

firstime = TRUE;
n = 0;

   fread(&RA[n],sizeof(double),1,datafile);
   fread(&DEC[n],sizeof(double),1,datafile);
   fread(&epoch,8,1,datafile);
   fread(&AST[n],sizeof(double),1,datafile);
//  printf("RA = %7.4f, DEC = %7.4f, EPOCH = %s, AST = %9.3f \n", RAHRS, DECDEG, epoch, AST);
   fread(&dummy,480,1,datafile);
   fread(&calon,  sizeof(struct PolSet),1,datafile);
   fread(&caloff, sizeof(struct PolSet),1,datafile);


 while(!feof(datafile)) {  
 


//  calculate observed Stokes parameters over the band for this time step;
   for(i=lowchan-k;i<=highchan+k;i++){
     cal.xx[i] = calon.xx[i] - caloff.xx[i];
     cal.yy[i] = calon.yy[i] - caloff.yy[i];
     cal.xy[i] = calon.xy[i] - caloff.xy[i];
     cal.yx[i] = calon.yx[i] - caloff.yx[i];
     ObsCal.I[i] = cal.xx[i] + cal.yy[i];
     ObsCal.Q[i] = cal.xx[i] - cal.yy[i];      
     ObsCal.U[i] = cal.xy[i] + cal.yx[i];
     ObsCal.V[i] = cal.xy[i] - cal.yx[i];
     ObsSky.I[i] = caloff.xx[i] + caloff.yy[i];
     ObsSky.Q[i] = caloff.xx[i] - caloff.yy[i];
     ObsSky.U[i] = caloff.xy[i] + caloff.yx[i];
     ObsSky.V[i] = caloff.xy[i] - caloff.yx[i];
   }

// find RFI by looking at first difference in frequency domain. Check both cal-off
// and cal-on signals
   for(i=lowchan;i<=highchan;i++){ 
     flag[i] = 1;
     diffxx = caloff.xx[i+1] - caloff.xx[i];
     diffyy = caloff.yy[i+1] - caloff.yy[i];
     if( (fabs(diffxx) > difflimit) || (fabs(diffyy) > difflimit) ) { 
       flag[i] = 0;
       if(AST[n] > 0.0) {
          fprintf(rfifile,"%4d %9.3f %8.2f %8.3f %8.3f\n", i, freq[i], AST[n], diffxx, diffyy);
       }
     }
     diffxx = calon.xx[i+1] - calon.xx[i];
     diffyy = calon.yy[i+1] - calon.yy[i];
     if( (fabs(diffxx) > difflimit) || (fabs(diffyy) > difflimit) ) { 
       if(flag[i] > 0) {
         flag[i] = 0;
         if(AST[n] > 0.0) {
           fprintf(rfifile,"%4d %9.3f %8.2f %8.3f %8.3f\n", i, freq[i], AST[n], diffxx, diffyy);
         }
       }
     }    
   }


// calculate gains
   for(i=lowchan;i<=highchan;i++){    
     calxx = 0.0; calyy= 0.0; calxy = 0.0; calyx = 0.0;
     count = 0;
     for(j=0;j<=smooth;j++) {
       if(RFI) {
         if (flag[i+j-k] > 0) {
           calxx = calxx + cal.xx[i+j-k];
           calyy = calyy + cal.yy[i+j-k];
           calxy = calxy + cal.xy[i+j-k];
           calyx = calyx + cal.yx[i+j-k];
           count = count + 1.0;
         }
       } else {
         calxx = calxx + cal.xx[i+j-k];
         calyy = calyy + cal.yy[i+j-k];
         calxy = calxy + cal.xy[i+j-k];
         calyx = calyx + cal.yx[i+j-k];
         count = count + 1.0;    
       }
     }
     calxx = calxx/count;
     calyy = calyy/count;
     calxy = calxy/count;
     calyx = calyx/count;
     gain.x[i] = sqrt(Tcal/calxx);
     gain.y[i] = sqrt(Tcal/calyy);
     calU = calxy + calyx;
     calV = calxy - calyx;
     gain.phi[i] = atan2(calV,calU); 
   }

// if(firstime) {
// for(i=lowchan;i<=highchan;i++){
//    printf("%4d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f \n", 
//   i+1,calon.xx[i],calon.yy[i],calon.xy[i],calon.yx[i],caloff.xx[i],caloff.yy[i],caloff.xy[i],
//   caloff.yx[i]);
//  }
//}

// calibrate the cal signal
  for(i=lowchan;i<=highchan;i++){
    CalCal.I[i] =  ((gain.x[i]*gain.x[i] + gain.y[i]*gain.y[i])/4.0)*ObsCal.I[i] 
                 + ((gain.x[i]*gain.x[i] - gain.y[i]*gain.y[i])/4.0)*ObsCal.Q[i];
    CalCal.Q[i] =  ((gain.x[i]*gain.x[i] - gain.y[i]*gain.y[i])/4.0)*ObsCal.I[i] 
                 + ((gain.x[i]*gain.x[i] + gain.y[i]*gain.y[i])/4.0)*ObsCal.Q[i];
    CalCal.U[i] = (gain.x[i]*gain.y[i]) * (cos(gain.phi[i])*ObsCal.U[i]
                    + sin(gain.phi[i])*ObsCal.V[i]);
    CalCal.V[i] = (gain.x[i]*gain.y[i]) * (-sin(gain.phi[i])*ObsCal.U[i]
                    +cos(gain.phi[i])*ObsCal.V[i]);
  }

  if(firstime) {
     for(i=lowchan;i<=highchan;i++){
       fprintf(calfile,"%4d %7.4f %7.4f %7.4f %7.4f %9.4f %7.4f %7.4f %7.4f %8.3f %6.3f %6.3f\n",
       i+1,      
       ObsCal.I[i],ObsCal.Q[i],ObsCal.U[i],ObsCal.V[i],CalCal.I[i],CalCal.Q[i],CalCal.U[i],
       CalCal.V[i],gain.x[i],gain.y[i],gain.phi[i]);
     }
    printf("\n");
   }   

// calibrate the Sky signal

  for(i=lowchan;i<=highchan;i++){
    TrueSky.I[i] = ((gain.x[i]*gain.x[i] + gain.y[i]*gain.y[i])/4.0)*ObsSky.I[i] 
                 + ((gain.x[i]*gain.x[i] - gain.y[i]*gain.y[i])/4.0)*ObsSky.Q[i];
    TrueSky.Q[i] = ((gain.x[i]*gain.x[i] - gain.y[i]*gain.y[i])/4.0)*ObsSky.I[i] 
                 + ((gain.x[i]*gain.x[i] + gain.y[i]*gain.y[i])/4.0)*ObsSky.Q[i];
    TrueSky.U[i] = (gain.x[i]*gain.y[i]) * (cos(gain.phi[i])*ObsSky.U[i]
                    + sin(gain.phi[i])*ObsSky.V[i]);
    TrueSky.V[i] = (gain.x[i]*gain.y[i]) * (-sin(gain.phi[i])*ObsSky.U[i]
                    + cos(gain.phi[i])*ObsSky.V[i]);
  }

  if(firstime) {
     for(i=lowchan;i<=highchan;i++){
       fprintf(skyfile,"%4d %7.4f %7.4f %7.4f %7.4f %9.4f %7.4f %7.4f %7.4f\n", i+1,      
       ObsSky.I[i],ObsSky.Q[i],ObsSky.U[i],ObsSky.V[i],TrueSky.I[i],TrueSky.Q[i],
       TrueSky.U[i],TrueSky.V[i]);
      }
   }  
 
   firstime = FALSE;
//  calculate average calibrated Stokes parameters over the band and
//  average x and y gains.
   I[n]=0;
   Q[n]=0;
   U[n]=0;
   V[n]=0;
   gx = 0.0;
   gy = 0.0; 
   goodchan = 0.0;
   for (i=lowchan;i<=highchan;i++){
     if(RFI) {
       if (flag[i] > 0) {
         I[n] = I[n] + TrueSky.I[i];
         Q[n] = Q[n] + TrueSky.Q[i];
         U[n] = U[n] + TrueSky.U[i];
         V[n] = V[n] + TrueSky.V[i];
         gx = gx + gain.x[i]*gain.x[i];
         gy = gy + gain.y[i]*gain.y[i];
         goodchan = goodchan + 1.0;
       }
     } else {
       I[n] = I[n] + TrueSky.I[i];
       Q[n] = Q[n] + TrueSky.Q[i];
       U[n] = U[n] + TrueSky.U[i];
       V[n] = V[n] + TrueSky.V[i]; 
       gx = gx + gain.x[i]*gain.x[i]/2.0;
       gy = gy + gain.y[i]*gain.y[i]/2.0;
       goodchan = goodchan + 1.0;   
     }   
   }
   I[n] = I[n]/goodchan;
   Q[n] = Q[n]/goodchan;
   U[n] = U[n]/goodchan;
   V[n] = V[n]/goodchan; 
   gx = gx/goodchan;
   gy = gy/goodchan;
   if(AST[n] > 0.0) {
     fprintf(gainfile,"%9.5f %9.4f %5.0f %9.4f %9.4f\n", RA[n], AST[n], goodchan, gx, gy); 
   }

   n++; 
   fread(&RA[n],sizeof(double),1,datafile);
   fread(&DEC[n],sizeof(double),1,datafile);
   fread(&epoch,8,1,datafile);
   fread(&AST[n],sizeof(double),1,datafile);
//  printf("RA = %7.4f, DEC = %7.4f, EPOCH = %s, AST = %9.3f \n", RAHRS, DECDEG, epoch, AST);
   fread(&dummy,480,1,datafile);
   fread(&calon,  sizeof(struct PolSet),1,datafile);
   fread(&caloff, sizeof(struct PolSet),1,datafile);
 
  }

 npoints = n-1;
 printf("\n Read %8d integrations \n", npoints);

// smooth the band averaged data in time
 for(n=kt; n< npoints-kt;n+=tsmooth){ 
   count = 0.0;
   dumI = 0.0;
   dumQ = 0.0;
   dumU = 0.0;
   dumV = 0.0;
   dumRA = 0.0;
   dumDEC = 0.0;
   dumt = 0.0;
   for(j=0;j<=tsmooth;j++) {
     if(!isnan(I[n+j-kt])) {
       dumI = dumI + I[n+j-kt];
       dumQ = dumQ + Q[n+j-kt];
       dumU = dumU + U[n+j-kt];
       dumV = dumV + V[n+j-kt];
       dumRA = dumRA + RA[n+j-kt];
       dumDEC = dumDEC + DEC[n+j-kt];
       dumt = dumt + AST[n+j-kt];
       count = count + 1.0;
     }
   }
   dumI = dumI/count;
   dumQ = dumQ/count;
   dumU = dumU/count;
   dumV = dumV/count;
   dumRA = dumRA/count;
   dumDEC = dumDEC/count;
   dumt = dumt/count;
   fprintf(fluxfile,"%7.4f %7.4f %9.3f %9.4f %8.4f %8.4f %8.4f \n", 
         dumRA, dumDEC, dumt, dumI, dumQ, dumU, dumV);
 }
 fclose(datafile);
 fclose(calfile);
 fclose(skyfile);
 fclose(fluxfile);
 fclose(rfifile);
 fclose(gainfile);

}















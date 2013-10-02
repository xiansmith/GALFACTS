/*
  Chebyshev Approximation Library
  Copyright (C) 2011, M. Andrecut, mandrecu@ucalgary.ca
  Institute for Space Imaging Science
  University of Calgary, Alberta, Canada

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY.  See the GNU General Public License for 
  more details: <http://www.gnu.org/licenses/>
*/

#ifndef _CHEBYSHEV_H
#define _CHEBYSHEV_H

#define cmax( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define cmin( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define csign(x) ( ((x>=0)-(x<0)) )
#define CNORMALIZE(x,min,max) ( ((x*2.0)-min-max) / (max-min) )
#define CDENORMALIZE(x,min,max) ( (x*(max-min) + max + min) / 2.0 )
#define id(n, m, N) (((m) * (N) + (n)))

void chebyshev_minmax(float x[], int N, float *min, float *max);

void chebyshev_normalize(float x[], int N, float min, float max);

void chebyshev_denormalize(float x[], int N, float min, float max);

float chebyshev_eval(float x, float c[], int order);

float chebyshev_eval_surface(float x, float y, float c[], int order);

float chebyshev(float x, int n);

void chebyshev_cholesky(int n, float c[], int N, float x[], float y[]);

void chebyshev_qrgsr(int n, float c[], int N, float x[], float y[]);

void chebyshev_qrgsm(int n, float c[], int N, float x[], float y[]);

void chebyshev_qrgsr_surface(int n, float c[], int N, float x[], float y[], float z[]);

void chebyshev_fit_bw(float X[], float Y[], int size, float nsigma, float C[], int order);

void chebyshev_fit_dec(float X[], float Y[], int size, float nsigma, float C[], int order);
void chebyshev_fit_sat(float X[], float Y[], int size, float nsigma, float C[], int order, int RFI[], float RA);


#endif

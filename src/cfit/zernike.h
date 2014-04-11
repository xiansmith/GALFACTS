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

#ifndef _ZERNIKE_H
#define _ZERNIKE_H

#define CNORMALIZE(x,min,max) ( ((x*2.0)-min-max) / (max-min) )
#define CDENORMALIZE(x,min,max) ( (x*(max-min) + max + min) / 2.0 )
#define id(n, m, N) (((m) * (N) + (n)))

void chebyshev_minmax(float x[], int N, float *min, float *max);

void chebyshev_normalize(float x[], int N, float min, float max);

void chebyshev_denormalize(float x[], int N, float min, float max);

float chebyshev_eval(float x, float c[], int order);

float zernike_eval_surface(float x, float y, float c[], int order);

float chebyshev(float x, int n);

float zernike(float x, float y, int n, int m);

float chebyshev_surface(float x, float y, int n, int m);

void chebyshev_qrgsr(int n, float c[], int N, float x[], float y[], float pfit_lambda);

void zernike_qrgsr_surface(int n, float c[], int N, float x[], float y[], float z[], float pfit_lambda);

#endif

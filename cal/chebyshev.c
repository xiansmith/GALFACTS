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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "chebyshev.h"


void chebyshev_minmax(float *x, int n, float *xmin, float *xmax)
{
/*
 Finds min (xmin) and max (xmax) for data vector x.
*/

	int i;
	*xmin = *xmax = x[0];
	for(i=1; i<n; i++) 
		{
		if(x[i] < *xmin) *xmin = x[i];
		else if(x[i] > *xmax) *xmax = x[i];
		}
}


void chebyshev_normalize(float *x, int N, float min, float max)
{
/*
 Transforms vector x in interval [min, max] to interval [-1, +1]
*/

	int n;
	for(n=0; n<N; n++) 
		{
		x[n] = CNORMALIZE(x[n],min,max);	
		}
}


void chebyshev_denormalize(float *x, int N, float min, float max)
{
/*
 Transforms vector x in interval [-1, 1] to interval [min, max]
*/

	int n;
	for(n=0; n<N; n++) 
		{
		x[n] = CDENORMALIZE(x[n],min,max);
		}
}


float chebyshev_eval(float x, float *c, int n)
{
/*
 Evaluates 1d Chebyshev expansion at point x using coefficients c 
 and polynomials of degree 0,1,...,n
*/

	int i;
	float y, p[n+3];

	p[n+2] = p[n+1] = 0;
	for(i=n; i>=0; i--) 
		{
		p[i] = c[i] + 2*x*p[i+1] - p[i+2];	
		}
	y = p[0] - x*p[1];
	return y;
}


float chebyshev(float x, int n)
{
/*
 Evaluates 1d Chebyshev polynomial of degree n at point x
*/

	int i;
	float p[n+1];

	p[0] = 1; 
	p[1] = x;
	for(i=2; i<=n; i++) 
		{
		p[i] = 2*x*p[i-1] - p[i-2]; 
		}
	return p[n];
}


float chebyshev_eval_surface(float x, float y, float *c, int n)
{
/*
 Evaluates 2d Chebyshev expansion (surface) at point (x,y) using 
 coefficients c and surface polynomials of degree 0,1,...,n
*/

	int t = 0, j, k;
	float p = 0;

	for(j=0; j<=n; j++) 
		{
		for(k=0; k<=n; k++)
			{
			p += c[t]*chebyshev(x, j)*chebyshev(y, k);
			t++;
			}
		}

	return p;
}


float chebyshev_surface(float x, float y, int n, int m)
{
/*
 Evaluates 2d Chebyshev surface of degree (n,m) at point (x,y)
*/

	return chebyshev(x, n)*chebyshev(y, m);
}


void chebyshev_qrgsr_surface(int n, float *c, int N, float *x, float *y, float *z)
{
/*
 2d Least Squares Fit using QR factorization with Gramm-Schmidt reorthogonalization.
 Fits data (x,y,z) using a 2d Chebyshev expansion of order n.
 Returns the coefficients c.
*/

	int i, j, k, nn = (n+1)*(n+1) - 1, t, nnn = nn + 1;
	float s;
	float *q = (float*)malloc(N*(nn+1)*sizeof(float));
	float *r = (float*)malloc((nn+1)*(nn+1) * sizeof(float));

	t = 0;
	for(j=0; j<=n; j++) 
		{
		for(k=0; k<=n; k++)
			{
			for(i=0; i<N; i++) 
				{
				q[i+t*N] = chebyshev_surface(x[i], y[i], j, k);
				}
			t++;
			}
		}

	for(j=0; j<=nn; j++)
		{	
		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[i+j*N]*q[i+k*N]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[i+j*N] -= s*q[i+k*N];
				}
			r[k+j*nnn] = s;
			}

		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[i+j*N]*q[i+k*N]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[i+j*N] -= s*q[i+k*N];
				}
			r[k+j*nnn] += s;
			}

		r[j+j*nnn] = 0; 
		for(i=0; i<N; i++) 
			{
			r[j+j*nnn] += q[i+j*N]*q[i+j*N]; 
			}
		s = r[j+j*nnn] = sqrt(r[j+j*nnn]);
		if(s > 0) s = 1.0/s;
		for(i=0; i<N; i++) 
			{
			q[i+j*N] = s*q[i+j*N];
			}
		}	
	for(k=0; k<=nn; k++)
		{
		c[k] = 0; 
		for(i=0; i<N; i++) 
			{
			c[k] += q[i+k*N]*z[i];
			}
		}
	if(r[nn+nn*nnn] > 0)
		{
		c[nn] = c[nn]/r[nn+nn*nnn];
		}
	else
		{
		c[nn] = 0;
		}
	for(k=nn-1; k>=0; k--)
		{
		s = 0; 
		for(i=k+1; i<=nn; i++) 
			{
			s += r[k+i*nnn]*c[i];
			}
		if(r[k+k*nnn] > 0)
			{
			c[k] = (c[k] - s)/r[k+k*nnn];
			}
		else
			{
			c[k] = 0;
			}
		}
	free(q); 
	free(r);
}


void chebyshev_qrgsr(int n, float *c, int N, float *x, float *y)
{
/*
 1d Least Squares Fit using QR factorization with Gramm-Schmidt reorthogonalization.
 Fits data (x,y) using a 1d Chebyshev expansion of order n.
 Returns the coefficients c.
*/

	float s;
	int i, j, k, nn=n+1;
	float *q = (float*)malloc(N*nn*sizeof(float));
	float *r = (float*)malloc(nn*nn*sizeof(float));

	for(i=0; i<N; i++) 
		{
		for(j=0; j<=n; j++) 
			{
			q[i+j*N] = chebyshev(x[i], j);
			}
		}
	for(j=0; j<=n; j++)
		{	
		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[i+j*N]*q[i+k*N]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[i+j*N] -= s*q[i+k*N];
				}
			r[k+j*nn] = s;
			}			
		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[i+j*N]*q[i+k*N]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[i+j*N] -= s*q[i+k*N];
				}
			r[k+j*nn] += s;
			}			
		r[j+j*nn] = 0;
		for(i=0; i<N; i++) 
			{
			r[j+j*nn] += q[i+j*N]*q[i+j*N]; 
			}
		s = r[j+j*nn] = sqrt(r[j+j*nn]);
		if(s > 0) s = 1.0/s;
		for(i=0; i<N; i++) 
			{
			q[i+j*N] = s*q[i+j*N];
			}
		}	
	for(k=0; k<=n; k++)
		{
		c[k] = 0; 
		for(i=0; i<N; i++) 
			{
			c[k] += q[i+k*N]*y[i];
			}
		}
	if(r[n+n*nn]>0)
		{
		c[n] = c[n]/r[n+n*nn];
		}
	else
		{
		c[n] = 0;
		}
	for(k=n-1; k>=0; k--)
		{
		s = 0; 
		for(i=k+1; i<=n; i++) 
			{
			s += r[k+i*nn]*c[i];
			}
		if(r[k+k*nn]>0)
			{
			c[k] = (c[k] - s)/r[k+k*nn];
			}
		else
			{
			c[k] = 0;
			}
		}
	free(q);
	free(r);
}


void chebyshev_qrgsm(int n, float *c, int N, float *x, float *y)
{
/*
 1d Least Squares Fit using QR factorization with Gramm-Schmidt reorthogonalization.
 Fits data (x,y) using a 1d Chebyshev expansion of order n.
 Returns the coefficients c.
*/

	float s;
	int i, j, k, nn=n+1;
	float *q = (float*)malloc(N*nn*sizeof(float));
	float *r = (float*)malloc(nn*nn*sizeof(float));

	if(q == NULL || r == NULL)
		printf("ERROR: malloc failed in chebyshev_cholesky()\n");

	for(i=0; i<N; i++) 
		{
		for(j=0; j<=n; j++) 
			{
			q[i+j*N] = chebyshev(x[i], j);
			}
		}
	for(j=0; j<=n; j++)
		{	
		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[i+j*N]*q[i+k*N]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[i+j*N] -= s*q[i+k*N];
				}
			r[k+j*nn] = s;
			}			
		r[j+j*nn] = 0;
		for(i=0; i<N; i++) 
			{
			r[j+j*nn] += q[i+j*N]*q[i+j*N]; 
			}
		s = r[j+j*nn] = sqrt(r[j+j*nn]);
		if(s > 0) s = 1.0/s;
		for(i=0; i<N; i++) 
			{
			q[i+j*N] = s*q[i+j*N];
			}
		}	
	for(k=0; k<=n; k++)
		{
		c[k] = 0; 
		for(i=0; i<N; i++) 
			{
			c[k] += q[i+k*N]*y[i];
			}
		}
	if(r[n+n*nn]>0)
		{
		c[n] = c[n]/r[n+n*nn];
		}
	else
		{
		c[n] = 0;
		}
	for(k=n-1; k>=0; k--)
		{
		s = 0; 
		for(i=k+1; i<=n; i++) 
			{
			s += r[k+i*nn]*c[i];
			}
		if(r[k+k*nn]>0)
			{
			c[k] = (c[k] - s)/r[k+k*nn];
			}
		else
			{
			c[k] = 0;
			}
		}

	free(q);
	free(r);
}


void chebyshev_cholesky(int n, float *c, int N, float *x, float *y)
{
/*
 1d Least Squares Fit using Cholesky factorization.
 Fits data (x,y) using a 1d Chebyshev expansion of order n.
 Returns the coefficients c.
*/

	float s;
	int i, j, k, nn = n + 1;
	float *q = (float*)malloc(N*nn*sizeof(float));
	float *a = (float*)malloc(nn*nn*sizeof(float));


	for(i=0; i<N; i++) 
		{
		for(j=0; j<=n; j++) 
			{
			q[i+j*N] = chebyshev(x[i], j);
			}
		}

	for(i=0; i<=n; i++)
		{
		for(j=i; j<=n; j++)
			{
			a[j+i*nn] = 0.0; for(k=0; k<N; k++) a[j+i*nn] += q[k+i*N]*q[k+j*N];
			}
		c[i] = 0.0; for(k=0; k<N; k++) c[i] += y[k]*q[k+i*N];
		}

	//for(i=0; i<=n; i++) a[i+i*nn] += 0.001; // regularization	just for very bad data

	for(i=0; i<=n; i++)
		{
		for(j=0; j<i; j++)
			{
			s = 0.0; 
			for(k=0; k<j; k++) 
				{
				s += a[i+k*nn]*a[j+k*nn]*a[k+k*nn];
				}
			a[i+j*nn] = (a[i+j*nn] - s)/a[j+j*nn];
			}
		s = 0.0; 
		for(k=0; k<j; k++) 
			{
			s += a[i+k*nn]*a[i+k*nn]*a[k+k*nn];
			}
		a[i+i*nn] = a[i+i*nn] - s;
		}
	
	for(i=1; i<=n; i++) 
		{
		for(j=0; j<i; j++) 
			{
			c[i] -= a[i+j*nn]*c[j];
			}
		}
	for(i=0; i<=n; i++) 
		{
		c[i] /= a[i+i*nn];
		}
	for(i=n-1; i>=0; i--) 
		{
		for(j=i+1; j<=n; j++) 
			{
			c[i] -= a[j+i*nn]*c[j];
			}
		}

	free(q);
	free(a);
}



void chebyshev_fit_bw(float *X, float *Y, int size, float nsigma, float *C, int order)
{
/*
 Basketweaving fitting function: fit + outlier removal
*/

	int i, outlier;
	float mean, sigma;

	float *x = (float*)malloc(size * sizeof(float));
	float *y = (float*)malloc(size * sizeof(float));
	float *p = (float*)malloc(size * sizeof(float));
	float *diff = (float*)malloc(size * sizeof(float));

	if(x == NULL || y == NULL || p == NULL || diff == NULL)
		printf("ERROR: malloc failed in chebyshev_fit_bw()\n");

	for(i=0; i<size; i++){x[i] = X[i]; y[i] = Y[i];}

	int iter=0;
	do 
	{
		//chebyshev_cholesky(order, C, size, x, y);
		chebyshev_qrgsm(order, C, size, x, y);
		mean = 0.0;
		for(i=0; i<size; i++) 
		{
			p[i] = chebyshev_eval(x[i], C, order);
			diff[i] = y[i] - p[i];
			mean += diff[i];
		}
		mean = mean/size;
		sigma = 0.0;
		for(i=0; i<size; i++) 
		{
			diff[i] = diff[i] - mean;
			sigma += diff[i]*diff[i];
		}
		sigma = nsigma*sqrt(sigma/size);
		outlier = 0;
		for(i=0; i<size; i++) 
		{
			if(fabs(diff[i]) > sigma) 
			{
				outlier = 1;
				size --; 
				x[i] = x[size]; 
				y[i] = y[size];
			}
		}
		iter++;
	}while(outlier && size > order);

	if(size <= order)
		printf("WARNING: size is %d order is %d\n",size,order);

	free(x); 
	free(y); 
	free(diff); 
	free(p);
}


void chebyshev_fit_dec(float *X, float *Y, int size, float nsigma, float *C, int order)
{
/*
 DEC background fitting function: fit + outlier removal
*/

	int i, M, outlier;
	float mean, sigma;

	float *xx = (float*)malloc(size * sizeof(float));
	float *yy = (float*)malloc(size * sizeof(float));
	float *x = (float*)malloc(size * sizeof(float));
	float *y = (float*)malloc(size * sizeof(float));
	float *p = (float*)malloc(size * sizeof(float));
	float *diff = (float*)malloc(size * sizeof(float));

	/* Extract lower envelope of data */
	//for(i=0; i<size; i++){x[i] = X[i]; y[i] = Y[i];} 
	xx[0] = X[0]; yy[0] = Y[0]; M = 1; for(i=1; i<size-1; i++) if(Y[i]-Y[i-1]<0 && Y[i+1]-Y[i]>0){xx[M] = X[i]; yy[M] = Y[i]; M++;} xx[M] = X[size-1]; yy[M] = Y[size-1]; size = M + 1;
	x[0] = xx[0]; y[0] = yy[0]; M = 1; for(i=1; i<size-1; i++) if(yy[i]-yy[i-1]>0 && yy[i+1]-yy[i]<0){x[M] = xx[i]; y[M] = yy[i]; M++;} x[M] = xx[size-1]; y[M] = yy[size-1]; size = M + 1;
	
	do 
		{
		//chebyshev_cholesky(order, C, size, x, y);
		chebyshev_qrgsm(order, C, size, x, y); 
		mean = 0.0;
		for(i=0; i<size; i++) 
			{ 
			p[i] = chebyshev_eval(x[i], C, order);
			diff[i] = y[i] - p[i];
			mean += diff[i];
			}
		mean = mean/size;

		sigma = 0.0;
		for(i=0; i<size; i++) 
			{
			diff[i] = diff[i] - mean;
			sigma += diff[i]*diff[i];
			}
		sigma = nsigma*sqrt(sigma/size);

		outlier = 0;	
		for(i=0; i<size; i++) 
			{
			if(fabs(diff[i]) > sigma) 
				{
				outlier = 1;
				size --; 
				x[i] = x[size]; 
				y[i] = y[size];
				}
			}
		}while(outlier && size > order);
	free(xx); 
	free(yy); 
	free(x); 
	free(y); 
	free(diff); 
	free(p);
}


void chebyshev_fit_sat(float *X, float *Y, int size, float nsigma, float *C, int order, int *RFI, float RA)
{
/*
 satellite RFI fitting function: fit + outlier removal
*/

	int i, M, outlier, outliercount = 0;
	float mean, sigma;

	float *x = (float*)malloc(size * sizeof(float));
	float *y = (float*)malloc(size * sizeof(float));
	float *p = (float*)malloc(size * sizeof(float));
	float *diff = (float*)malloc(size * sizeof(float));
	for(i=0; i<size; i++){x[i] = X[i]; y[i] = Y[i];}

	do
		{
		//chebyshev_cholesky(order, C, size, x, y);
		chebyshev_qrgsm(order, C, size, x, y);
		mean = 0.0;
		for(i=0; i<size; i++)
			{
			p[i] = chebyshev_eval(x[i], C, order);
			diff[i] = y[i] - p[i];
			mean += diff[i];
			}
		mean = mean/size;

		sigma = 0.0;
		for(i=0; i<size; i++)
			{
			diff[i] = diff[i] - mean;
			sigma += diff[i]*diff[i];
			}
		sigma = nsigma*sqrt(sigma/size);

		//printf("size = %d ", size );
		/*if ( fabs(RA - 71.390589) < 0.001 )
		{
		char filename[50];
		sprintf(filename, "fit_%f.dat", RA);
		FILE *fit = fopen( filename , "w" );
		for(i=0; i<size; i++)
		{
			fprintf(fit, "%f", diff[i] );
			fprintf(fit, "\n");
		}
		fclose( fit );
		}
		*/

		outlier = 0;
		for(i=0; i<size; i++)
			{
			if(fabs(diff[i]) > sigma)
				{
				outlier = 1;
				outliercount++;
				RFI[i] = 1.0;
				size --;
				x[i] = x[size];
				y[i] = y[size];
				}
			}
		}while(outlier && size > order);
	free(x);
	free(y);
	free(diff);
	free(p);



	//printf("Chevy outlier = %d\n", outliercount );
}

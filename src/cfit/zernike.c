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
#include <float.h>
//#include "chebyshev.h"
#include "zernike.h"


void chebyshev_minmax(float x[], int n, float *xmin, float *xmax)
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


void chebyshev_normalize(float x[], int N, float min, float max)
{
/*
 Transforms vector x in interval [min, max] to interval [-1, +1]
*/

	int n;
	for(n=0; n<N; n++) 
		{
		x[n] = CNORMALIZE(x[n], min, max);
		}
}


void chebyshev_denormalize(float x[], int N, float min, float max)
{
/*
 Transforms vector x in interval [-1, 1] to interval [min, max]
*/

	int n;
	for(n=0; n<N; n++) 
		{
		x[n] = CDENORMALIZE(x[n], min, max);		
		}
}


float chebyshev_eval(float x, float c[], int n)
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


float chebyshev_eval_surface(float x, float y, float c[], int n)
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

static int factorial(int n)
{
        int i,fact;
        if(n==0)
                return 1;
        fact = 1;
        for(i = 1;i<=n;i++)
                fact*=i;
        return(fact);
}

//void zernike_pt(int n, int m, float x, float y, float *z,float radius)
float zernike(float x, float y, int n, int m)
{
	//float radius = 1.0;
	float radius = 2.0;
	float z;
        int j,jmax;
        float r, theta, norm,max;
        jmax = 0.5*(n-abs(m));
        norm=sqrt(2*(n+1)/(1+(m==0)));

        max = 0.0;
        z = 0;
        if(n == 0)
                z = 1;
        else if((n-m)%2)
                z = 0;
        else
        {
                r = sqrt(x*x+y*y);
                if((x>=0 && y>=0) || (x>=0 && y<0))
                        theta = atan(y/x);
                else
                        theta = M_PI + atan(y/x);
                for(j = 0;j <= jmax;j++)
                {
                        z+=pow(-1,j)*factorial(n-j)*pow(r/radius,n-2*j)/(factorial(j)*factorial(0.5*(n+abs(m))-j)*factorial(0.5*(n-abs(m))-j));
                }
                z=norm*(z)*((m>=0)*cos(m*theta)-(m<0)*sin(m*theta));
        }
	//printf("%f\n",z);
	return z;
}

float zernike_eval_surface(float x, float y, float c[], int n)
{
	int j;

	int a=0,b=0;
	float p=0;	
	for(j=0; j<=n; j++) 
		{

		p+=c[j]*zernike(x,y,a,b);
                if(a==b)
			{
                        a++;
                        b = -a;
                	}
                else
                        b+=2;

		}
	return p;
}


void zernike_qrgsr_surface(int n, float c[], int N, float x[], float y[], float z[], float pfit_lambda)
{
/*
 2d Least Squares Fit using QR factorization with Gramm-Schmidt reorthogonalization.
 Fits data (x,y,z) using a 2d Chebyshev expansion of order n.
 Returns the coefficients c.
*/

	//int i, j, k, nn = (n+1)*(n+1)-1, t;
	int i, j, k, t;
	float s;

	float *q = (float*)malloc(N*(n+1)*sizeof(float));

	//float *r = (float*)malloc((nn+1)*(nn+1) * sizeof(float));
	float *r = (float*)malloc((n+1)*(n+1) * sizeof(float));

	t = 0;

	int a=0,b=0;
	for(j=0; j<=n; j++) 
		{

//		for(k=0; k<=n; k++)
//			{
		for(i=0; i<N; i++) 
			{
			q[id(i,t,N)] = zernike(x[i], y[i], a, b);
			}
		t++;
//			}

                if(a==b)
			{
                        a++;
                        b = -a;
                	}
                else
                        b+=2;

		}

	for(j=0; j<=n; j++)
		{	
		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[id(i,j,N)]*q[id(i,k,N)]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[id(i,j,N)] -= s*q[id(i,k,N)];
				}
			r[id(k,j,n+1)] = s;
			}
		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[id(i,j,N)]*q[id(i,k,N)]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[id(i,j,N)] -= s*q[id(i,k,N)];
				}
			r[id(k,j,n+1)] += s;
			}
		r[id(j,j,n+1)] = 0; 
		for(i=0; i<N; i++) 
			{
			r[id(j,j,n+1)] += q[id(i,j,N)]*q[id(i,j,N)]; 
			}
		r[id(j,j,n+1)] = sqrt(r[id(j,j,n+1)]) + pfit_lambda; //regularization
		for(i=0; i<N; i++) 
			{
			q[id(i,j,N)] = q[id(i,j,N)]/r[id(j,j,n+1)];
			}
		}	
	for(k=0; k<=n; k++)
		{
		c[k] = 0; 
		for(i=0; i<N; i++) 
			{
			c[k] += q[id(i,k,N)]*z[i];
			}
		}
	c[n] = c[n]/r[id(n,n,n+1)];
	for(k=n-1; k>=0; k--)
		{
		s = 0; 
		for(i=k+1; i<=n; i++) 
			{
			s += r[id(k,i,n+1)]*c[i];
			}
		c[k] = (c[k] - s)/r[id(k,k,n+1)];
		}
	free(q); free(r);
}


void chebyshev_qrgsr(int n, float c[], int N, float x[], float y[], float pfit_lambda)
{
/*
 1d Least Squares Fit using QR factorization with Gramm-Schmidt reorthogonalization.
 Fits data (x,y) using a 1d Chebyshev expansion of order n.
 Returns the coefficients c.
*/

	int i, j, k;
	float *q = (float*)malloc(N*(n+1)*sizeof(float));
	float r[n+1][n+1], s;

	for(i=0; i<N; i++) 
		{
		for(j=0; j<=n; j++) 
			{
			q[id(i,j,N)] = chebyshev(x[i], j);
			}
		}
	for(j=0; j<=n; j++)
		{	
		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[id(i,j,N)]*q[id(i,k,N)]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[id(i,j,N)] -= s*q[id(i,k,N)];
				}
			r[k][j] = s;
			}
		for(k=0; k<j; k++) 
			{
			s = 0; 
			for(i=0; i<N; i++) 
				{
				s += q[id(i,j,N)]*q[id(i,k,N)]; 			
				}
			for(i=0; i<N; i++) 
				{
				q[id(i,j,N)] -= s*q[id(i,k,N)];
				}
			r[k][j] += s;
			}
		r[j][j] =r[j][j] = 0; 
		for(i=0; i<N; i++) 
			{
			r[j][j] += q[id(i,j,N)]*q[id(i,j,N)]; 
			}
		r[j][j] = sqrt(r[j][j]) + pfit_lambda; //regularization
		for(i=0; i<N; i++) 
			{
			q[id(i,j,N)] = q[id(i,j,N)]/r[j][j];
			}
		}	
	for(k=0; k<=n; k++)
		{
		c[k] = 0; 
		for(i=0; i<N; i++) 
			{
			c[k] += q[id(i,k,N)]*y[i];
			}
		}
	c[n] = c[n]/r[n][n];
	for(k=n-1; k>=0; k--)
		{
		s = 0; for(i=k+1; i<=n; i++) s += r[k][i]*c[i];
		c[k] = (c[k] - s)/r[k][k];
		}
	free(q);
}

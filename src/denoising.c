#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "denoising.h"

//----------------------------------------------------------------------------------------------------------------------------------------
void andtv_filter(float *x, int N, float tau, float lambda)
{
// Denoising filter:
// Adaptive nonlinear diffusion with total variation regularization
// x = the signal to be filtered, returns filtered signal
// N = length of the signal
// tau = time step (recommended values: 0.001 <= tau <= 0.1)
// lambda = regularization (recommended values: 0<=tau<=3)

int n, i, i_max=1000;
float ss, ss_max, q=tau, hh, hr, hl, D, R;

float *y; y = (float*)malloc(N*sizeof(float));
float *z; z = (float*)malloc(N*sizeof(float));

for(n=0; n<N; n++) z[n] = x[n];
ss_max = 0; 
for(n=0; n<N-1; n++) 
	{
	ss_max += x[n]*x[n] + 0.25*fabs(x[n+1] - x[n]); 
	}
ss_max += x[N-1]*x[N-1];
i = 0;
while(i < i_max)
	{
	i++;				
	for(n=2; n<N-2; n++) 
		{
		hr = 1.0/(0.5*fabs(x[n+2] - x[n]) + q); 
		hh = 1.0/(0.5*fabs(x[n+1] - x[n-1]) + q);
		hl = 1.0/(0.5*fabs(x[n] - x[n-2]) + q);
		D = 0.25*(hr + hh)*(x[n+1] - x[n]) - 0.25*(hh + hl)*(x[n] - x[n-1]);
		R = lambda*(x[n] - z[n]);
		y[n] = x[n] + tau*(D - R);
		}
	hr = 1.0/(0.5*fabs(x[2] - x[1]) + q); 
	hh = 1.0/(0.5*fabs(x[1] - x[0]) + q);
	hl = hh;  
	D = 0.25*(hr + hh)*(x[2] - x[1]) - 0.25*(hh + hl)*(x[1] - x[0]);
	R = lambda*(x[1] - z[1]);
	y[1] = x[1] + tau*(D - R);
	hh = 1.0/(0.5*fabs(x[N-1] - x[N-3]) + q);
	hr = hh; 
	hl = 1.0/(0.5*fabs(x[N-2] - x[N-4]) + q);  ;  
	D = 0.25*(hr + hh)*(x[N-1] - x[N-2]) - 0.25*(hh + hl)*(x[N-2] - x[N-3]);
	R = lambda*(x[N-2] - z[N-2]);
	y[N-2] = x[N-2] + tau*(D - R);
	D = 0.5*(x[1] - x[0]);
	R = lambda*(x[0] - z[0]);
	y[0] = x[0] + tau*(D - R);
	D = 0.5*(x[N-2] - x[N-1]);
	R = lambda*(x[N-1] - z[N-1]);
	y[N-1] = x[N-1] + tau*(D - R);							
	ss = 0; 
	for(n=0; n<N-1; n++) 
		{
		ss += (z[n] - y[n])*(z[n] - y[n]) + 0.25*fabs(y[n+1] - y[n]);
		} 
	ss += (z[N-1] - y[N-1])*(z[N-1] - y[N-1]);
	if(ss<ss_max)
		{
		ss_max = ss; 
		for(n=0; n<N; n++) x[n] = y[n]; 
		} 
	else break;
	}

free(y); free(z);
}
//----------------------------------------------------------------------------------------------------------------------------------------
void aldtv_filter(float *x, int N, float tau, float lambda)
{
// Denoising filter:
// Adaptive linear diffusion with total variation regularization
// x = the signal to be filtered, returns filtered signal
// N = length of the signal
// tau = time step (recommended values: 0.01 <= tau <= 1)
// lambda = regularization (recommended values: 0<=tau<=5)

int n, i, i_max=100;
float ss, ss_max;

float *y; y = (float*)malloc(N*sizeof(float));
float *z; z = (float*)malloc(N*sizeof(float));

for(n=0; n<N; n++) z[n] = x[n];
ss_max = 0; 
for(n=0; n<N-1; n++) 
	{
	ss_max += x[n]*x[n] + 0.25*fabs(x[n+1] - x[n]); 
	}
ss_max += x[N-1]*x[N-1];

i = 0;
while(i < i_max)
	{
	i++;	
	for(n=1; n<N-1; n++) 
		{
		y[n] = x[n] + 0.25*tau*(x[n-1] + x[n+1] - 2.0*x[n] - lambda*(x[n] - z[n]));
		}
	y[0] = x[0] + tau*(0.5*(x[1] - x[0]) - lambda*(x[0] - z[0]));
	y[N-1] = x[N-1] + tau*(0.5*(x[N-2] - x[N-1]) - lambda*(x[N-1] - z[N-1]));			
	ss = 0; 
	for(n=0; n<N-1; n++) 
		{
		ss += (z[n] - y[n])*(z[n] - y[n]) + 0.25*fabs(y[n+1] - y[n]); 
		}
	ss += (z[N-1] - y[N-1])*(z[N-1] - y[N-1]);
	if(ss<ss_max)
		{
		ss_max = ss; 
		for(n=0; n<N; n++) x[n] = y[n]; 
		} 
	else 
		{
		break;
		}
	}
free(y); free(z);
}
//----------------------------------------------------------------------------------------------------------------------------------------
void diffusion_filter(float *x, int N, int T)
{
// Denoising filter:
// Adaptive linear diffusion with total variation regularization
// x = the signal to be filtered, returns filtered signal
// N = length of the signal

int n, i, i_max = T;
float ss, ss_max;

float *y; y = (float*)malloc(N*sizeof(float));
float *z; z = (float*)malloc(N*sizeof(float));

for(n=0; n<N; n++) z[n] = x[n];

ss_max = 0; 
for(n=0; n<N-1; n++) 
	{
	ss_max += x[n]*x[n] + 0.25*fabs(x[n+1] - x[n]); 
	}
ss_max += x[N-1]*x[N-1];

i = 0;
while(i < i_max)
	{
	i++;	
	for(n=1; n<N-1; n++) 
		{
		y[n] = 0.5*x[n] + 0.25*(x[n-1] + x[n+1]);
		}
	y[0] = 0.5*(x[1] + x[0]);
	y[N-1] = 0.5*(x[N-2] + x[N-1]);			
	ss = 0; 
	for(n=0; n<N-1; n++) 
		{
		ss += (z[n] - y[n])*(z[n] - y[n]) + 0.25*fabs(y[n+1] - y[n]); 
		}
	ss += (z[N-1] - y[N-1])*(z[N-1] - y[N-1]);
	if(ss<ss_max)
		{
		ss_max = ss; 
		for(n=0; n<N; n++) x[n] = y[n]; 
		} 
	else 
		{
		//printf("%d\n", i);
		break;
		}
	}
free(y); 
free(z);
}
//----------------------------------------------------------------------------------------------------------------------------------------
void moving_average_filter(float *x, int N, int L)
{
// Denoising filter:
// L = window size
// x = the signal to be filtered, returns filtered signal
// N = length of the signal

int n, m, i, j, t, LL = 2*L+1;
float mean;

float *y; y = (float*)malloc(N*sizeof(float));


for(n=0; n<=L; n++)
	{
	mean = 0; for(i=-n; i<=L; i++) mean += x[n+i]; mean /= (L+n+1);
	y[n] = 0;
	for(i=-L; i<=L; i++)
		{
		if(n+i<0) y[n] += mean; else y[n] += x[n+i];
		}
	y[n] /= LL;
	}
for(n=L+1; n<N-1-L; n++)
	{
	y[n] = y[n-1] + (x[n+L] - x[n-L-1])/LL;
	}
for(n=N-1-L; n<N; n++)
	{
	mean = 0; for(i=-L; i<=N-1-n; i++) mean += x[n+i]; mean /=(L+1+N-1-n);
	y[n] = 0;
	for(i=-L; i<=L; i++)
		{
		if(n+i>N-1) y[n] += mean; else y[n] += x[n+i];
		}
	y[n] /= LL;
	}	
for(n=0; n<N; n++) x[n] = y[n];


free(y); 
}
//----------------------------------------------------------------------------------------------------------------------------------------
void gaussian_filter(float *x, int N, float *w, int L)
{
// Denoising filter:
// L = window size
// x = the signal to be filtered, returns filtered signal
// N = length of the signal

int n, m, i, j, t, LL = 2*L+1;
float meanL, meanR;

float *y; y = (float*)malloc(N*sizeof(float));


meanL = 0; for(i=0; i<=L; i++) meanL += x[i]; meanL /= (L+1);
meanR = 0; for(i=N-1-L; i<=N-1; i++) meanR += x[i]; meanR /=(L+1);
	
for(n=0; n<N; n++)
	{	
	y[n] = 0;
	for(i=-L; i<=L; i++)
		{
		if(n+i>=0 && n+i<=N-1) y[n] += x[n+i]*w[abs(i)]; else if(n+i<0) y[n] += meanL*w[abs(i)]; else y[n] += meanR*w[abs(i)]; 	
		}
	}

for(n=0; n<N; n++) x[n] = y[n];


free(y); 
}
//----------------------------------------------------------------------------------------------------------------------------------------

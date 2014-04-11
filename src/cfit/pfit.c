//old+new
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "pfit.h"
//-----------------------------------------------------------------------
void pfit_minmax(float x[], int n, float *xmin, float *xmax)
{
	int i;
	*xmin = *xmax = x[0];
	for(i=1; i<n; i++) 
		{
		if(x[i] < *xmin) *xmin = x[i];
		else if(x[i] > *xmax) *xmax = x[i];
		}
}
//-----------------------------------------------------------------------
void pfit_normalize(float x[], int N, float min, float max)
{
	int n;
	for(n=0; n<N; n++) x[n] = NORMALIZE(x[n], min, max);
}
//-----------------------------------------------------------------------

void pfit_denormalize(float x[], int N, float min, float max)
{
	int n;
	for (n=0; n<N; n++) x[n] = DENORMALIZE(x[n], min, max);		
}
//-----------------------------------------------------------------------
float pfit_poly_eval(float x, float c[], int n)
{
// n - the order of the polynomial
// c - the coefficients (must be of size n+1) [0..n]
// x - the value to evaluate for
	int i;
	float y;
	y = c[n]; for(i=n-1; i>=0; i--) y = y*x + c[i];
	return y;
}
//-----------------------------------------------------------------------
float pw(float x, int n)
{
	int i;
	float y = 1;
	for(i=1; i<=n; i++) y = y*x; 
	return y;
}
//-----------------------------------------------------------------------
void pfit_qrgsr(int n, float c[], int N, float x[], float y[], float pfit_lambda)
{
// n = degree of polynomial, n<N
// c = coefficients of polynomial
// N = nr. of data points
// x, y = data
int i, j, k;
float *q = (float*)malloc(N*(n+1)*sizeof(float));
float r[n+1][n+1], s;

for(i=0; i<N; i++) for(j=0; j<=n; j++) q[id(i,j,N)] = pw(x[i],j);

for(j=0; j<=n; j++)
	{	
	for(k=0; k<j; k++) 
		{
		s = 0; for(i=0; i<N; i++) s += q[id(i,j,N)]*q[id(i,k,N)]; 			
		for(i=0; i<N; i++) q[id(i,j,N)] -= s*q[id(i,k,N)];
		r[k][j] = s;
		}
	for(k=0; k<j; k++) 
		{
		s = 0; for(i=0; i<N; i++) s += q[id(i,j,N)]*q[id(i,k,N)]; 			
		for(i=0; i<N; i++) q[id(i,j,N)] -= s*q[id(i,k,N)];
		r[k][j] += s;
		}
	r[j][j] =r[j][j] = 0; for(i=0; i<N; i++) r[j][j] += q[id(i,j,N)]*q[id(i,j,N)]; 
	r[j][j] = sqrt(r[j][j]) + pfit_lambda;
	for(i=0; i<N; i++) q[id(i,j,N)] = q[id(i,j,N)]/r[j][j];
	}	

for(k=0; k<=n; k++){c[k] = 0; for(i=0; i<N; i++) c[k] += q[id(i,k,N)]*y[i];}

c[n] = c[n]/r[n][n];
for(k=n-1; k>=0; k--)
	{
	s = 0; for(i=k+1; i<=n; i++) s += r[k][i]*c[i];
	c[k] = (c[k] - s)/r[k][k];
	}
free(q);
}
//-----------------------------------------------------------------------
void pfit_qrgsr_adaptive(int n, float c[], int N, float x[], float y[])
{
// n = degree of polynomial, n<N
// c = coefficients of polynomial
// N = nr. of data points
// x, y = data
int i, j, k, tmax = 21, t = 0;
float *a = (float*)malloc(N*(n+1)*sizeof(float));
float *q = (float*)malloc(N*(n+1)*sizeof(float));
float r[n+1][n+1], cc[n+1], s, lambda = 0.0000001, max_chisq = FLT_MAX, chisq;

for(i=0; i<N; i++) for(j=0; j<=n; j++) a[id(i,j,N)] = pw(x[i],j);
while(t<tmax)
	{
	t++;
	for(i=0; i<N; i++) for(j=0; j<=n; j++) q[id(i,j,N)] = a[id(i,j,N)];
	for(j=0; j<=n; j++)
		{	
		for(k=0; k<j; k++) 
			{
			s = 0; for(i=0; i<N; i++) s += q[id(i,j,N)]*q[id(i,k,N)]; 			
			for(i=0; i<N; i++) q[id(i,j,N)] -= s*q[id(i,k,N)];
			r[k][j] = s;
			}
		for(k=0; k<j; k++) 
			{
			s = 0; for(i=0; i<N; i++) s += q[id(i,j,N)]*q[id(i,k,N)]; 			
			for(i=0; i<N; i++) q[id(i,j,N)] -= s*q[id(i,k,N)];
			r[k][j] += s;
			}
		r[j][j] =r[j][j] = 0; for(i=0; i<N; i++) r[j][j] += q[id(i,j,N)]*q[id(i,j,N)]; 
		r[j][j] = sqrt(r[j][j]) + lambda;
		for(i=0; i<N; i++) q[id(i,j,N)] = q[id(i,j,N)]/r[j][j];
		}	
	for(k=0; k<=n; k++){cc[k] = 0; for(i=0; i<N; i++) cc[k] += q[id(i,k,N)]*y[i];}
	cc[n] = cc[n]/r[n][n];
	for(k=n-1; k>=0; k--)
		{
		s = 0; for(i=k+1; i<=n; i++) s += r[k][i]*cc[i];
		cc[k] = (cc[k] - s)/r[k][k];
		}
	chisq = 0;
	for(i=0; i<N; i++)
		{
		s = pfit_poly_eval(x[i], cc, n);
		s = y[i] - s;
		chisq += s*s;
		}
	if(chisq < max_chisq)
		{
		max_chisq = chisq;
		lambda = lambda*2;
		for(i=0; i<=n; i++) c[i] = cc[i];
		}
		else break;
	}
free(a); free(q);
}
//-----------------------------------------------------------------------
void pfit_ldl(int n, float c[], int N, float x[], float y[], float pfit_lambda)
{
// n = degree of polynomial, n<N
// c = coefficients of polynomial
// N = nr. of data points
// x, y = data
int i, j, k;
float F[n+1][n+1], ftmp;
for(i=0; i<=n; i++)
	{
	for(j=i; j<=n; j++)
		{
		F[j][i] = 0.0; for(k=0; k<N; k++) F[j][i] += pw(x[k], i+j);
		}
	c[i] = 0.0; for(k=0; k<N; k++) c[i] += y[k]*pw(x[k], i);
	}
for(i=0; i<=n; i++) F[i][i] += pfit_lambda; //regularization
for(i=0; i<=n; i++)
	{
	for(j=0; j<=i; j++)
		{
		ftmp = 0.0; for(k=0; k<j; k++) ftmp += F[i][k]*F[j][k]*F[k][k];
		if(i==j) F[j][j] = F[j][j] - ftmp; else F[i][j] = (F[i][j] - ftmp)/F[j][j];
		}
	}
for(i=1; i<=n; i++) for(j=0; j<i; j++) c[i] -= F[i][j]*c[j];
for(i=0; i<=n; i++) c[i] /= F[i][i];
for(i=n-1; i>=0; i--) for(j=i+1; j<=n; j++) c[i] -= F[j][i]*c[j];
}
//-----------------------------------------------------------------------
void pfit_ldl_adaptive(int n, float c[], int N, float x[], float y[])
{
// n = degree of polynomial, n<N
// c = coefficients of polynomial
// N = nr. of data points
// x, y = data
int i, j, k, t=0, tmax = 21;
float A[n+1][n+1], F[n+1][n+1], aa[n+1], cc[n+1], ftmp, chisq, max_chisq = FLT_MAX, lambda = 1.0;

for(i=0; i<=n; i++)
	{
	for(j=i; j<=n; j++)
		{
		A[j][i] = 0.0; for(k=0; k<N; k++) A[j][i] += pw(x[k], i+j);
		}
	aa[i] = 0.0; for(k=0; k<N; k++) aa[i] += y[k]*pw(x[k], i);
	}
while(t<tmax)
	{
	t++;
	for(i=0; i<=n; i++) for(j=0; j<=n; j++) F[i][j] = A[i][j];
	for(i=0; i<=n; i++) F[i][i] += lambda; //regularization
	for(i=0; i<=n; i++) cc[i] = aa[i];
	for(i=0; i<=n; i++)
		{
		for(j=0; j<=i; j++)
			{
			ftmp = 0.0; for(k=0; k<j; k++) ftmp += F[i][k]*F[j][k]*F[k][k];
			if(i==j) F[j][j] = F[j][j] - ftmp; else F[i][j] = (F[i][j] - ftmp)/F[j][j];
			}
		}
	for(i=1; i<=n; i++) for(j=0; j<i; j++) cc[i] -= F[i][j]*cc[j];
	for(i=0; i<=n; i++) cc[i] /= F[i][i];
	for(i=n-1; i>=0; i--) for(j=i+1; j<=n; j++) cc[i] -= F[j][i]*cc[j];
	chisq = 0;
	for(i=0; i<N; i++)
		{
		ftmp = pfit_poly_eval(x[i], cc, n);
		ftmp = y[i] - ftmp;
		chisq += ftmp*ftmp;
		}
	if(chisq < max_chisq)
		{
		max_chisq = chisq;
		lambda = lambda*0.5;
		for(i=0; i<=n; i++) c[i] = cc[i];
		}
		else break;
	}
}
//-----------------------------------------------------------------------
void pfit_lsq(int n, float c[], int N, float x[], float y[], float pfit_lambda)
{
// n = degree of polynomial, n<N
// c = coefficients of polynomial
// N = nr. of data points
// x, y = data
int row[n+1], itmp, i, j, k, p, q;
float F[n+1][n+2], ftmp;

for(i=0; i<=n; i++)
	{
	for(j=i; j<=n; j++)
		{
		F[i][j] = 0.0; for(k=0; k<N; k++) F[i][j] += pw(x[k], i+j);
		}
	F[i][n+1] = 0.0; for(k=0; k<N; k++) F[i][n+1] += y[k]*pw(x[k], i);
	}
for(i=0; i<n; i++)
	{
	for(j=i+1; j<=n; j++)
		{
		F[j][i] = F[i][j];
		}
	}
for(i=0; i<=n; i++) F[i][i] += pfit_lambda; //regularization

for(j=0; j<=n; j++) row[j] = j;
for(p=0; p<n; p++)
	{
	for(k=p+1; k<=n; k++)
		{
		if(fabs(F[row[k]][p]) > fabs(F[row[p]][p]))
			{
			itmp = row[p]; 
			row[p] = row[k]; 
			row[k] = itmp;
			}
		}
	for(k=p+1; k<=n; k++)
		{
		ftmp = F[row[k]][p]/F[row[p]][p];
	    for(q=p+1; q<=n+1; q++) F[row[k]][q] -= ftmp*F[row[p]][q];
		}
	}
c[n] = F[row[n]][n+1]/F[row[n]][n];
for(k=n-1; k>=0; k--)
	{
	ftmp=0; for(q=k+1; q<=n; q++) ftmp += F[row[k]][q]*c[q];
	c[k] = (F[row[k]][n+1] - ftmp)/F[row[k]][k];
	}
}
//-----------------------------------------------------------------------
void pfit_lsq_adaptive(int n, float c[], int N, float x[], float y[])
{
int row[n+1], itmp, t_max = 21, t=0;
float A[n+1][n+2], F[n+1][n+2], cc[n+1], ftmp, chisq, max_chisq = FLT_MAX, lambda = 1.0;
int i, j, k, p, q;

for(i=0; i<=n; i++)
	{
	for(j=i; j<=n; j++)
		{
		A[i][j] = 0.0; for(k=0; k<N; k++) A[i][j] += pw(x[k], i+j);
		}
	A[i][n+1] = 0.0; for(k=0; k<N; k++) A[i][n+1] += y[k]*pw(x[k], i);
	}
for(i=0; i<n; i++) for(j=i+1; j<=n; j++) A[j][i] = A[i][j];

while(t<t_max)
	{
	t++;
	for(i=0; i<=n; i++) for(j=0; j<=n+1; j++) F[i][j] = A[i][j];
	for(i=0; i<=n; i++) F[i][i] += lambda; //regularization
	for(j=0; j<=n; j++) row[j] = j;
	for(p=0; p<n; p++)
		{
		for(k=p+1; k<=n; k++)
			{
			if(fabs(F[row[k]][p]) > fabs(F[row[p]][p]))
				{
				itmp = row[p];
				row[p] = row[k];
				row[k] = itmp;
				}
			}
		for(k=p+1; k<=n; k++)
			{
			ftmp = F[row[k]][p]/F[row[p]][p];
			for(q=p+1; q<=n+1; q++) F[row[k]][q] -= ftmp*F[row[p]][q];
			}
		}
	cc[n] = F[row[n]][n+1]/F[row[n]][n];
	for(k=n-1; k>=0; k--)
		{
		ftmp=0; for(q=k+1; q<=n; q++) ftmp += F[row[k]][q]*cc[q];
		cc[k] = (F[row[k]][n+1] - ftmp)/F[row[k]][k];
		}
	chisq = 0;
	for(i=0; i<N; i++)
		{
		ftmp = pfit_poly_eval(x[i], cc, n);
		ftmp = y[i] - ftmp;
		chisq += ftmp*ftmp;
		}
	if(chisq < max_chisq)
		{
		max_chisq = chisq;
		lambda = lambda*0.5;
		for(i=0; i<=n; i++) c[i] = cc[i];
		}
		else break;
	}
}
//-----------------------------------------------------------------------
int pfit_poly_fit(float X[], float Y[], int size, float nsigma, float C[], int order, int pfit_type, float pfit_lambda)
{
	int i;
	float sigma;
	int outlier;
	float mean;
	int err = 0; //assume no error

	float *x = (float*)malloc(size * sizeof(float));
	float *y = (float*)malloc(size * sizeof(float));
	float *diff = (float*)malloc(size * sizeof(float));

	for(i=0; i<size; i++){x[i] = X[i]; y[i] = Y[i];}

	do 
		{
		switch(pfit_type)
			{
			case 0: pfit_qrgsr(order, C, size, x, y, pfit_lambda); break;
			case 1: pfit_qrgsr_adaptive(order, C, size, x, y); break;
			case 2: pfit_lsq(order, C, size, x, y, pfit_lambda); break;
			case 3: pfit_lsq_adaptive(order, C, size, x, y); break;
			case 4: pfit_ldl(order, C, size, x, y, pfit_lambda); break;
			case 5: pfit_ldl_adaptive(order, C, size, x, y); break;
			default: pfit_qrgsr(order, C, size, x, y, pfit_lambda); break;
			}

		/* do the evaluations */
		mean = 0.0;
		for(i=0; i<size; i++) 
			{ 
			diff[i] = y[i] - pfit_poly_eval(x[i], C, order);
			mean += diff[i];
			}
		mean = mean/size;

		/* calc sigma */
		sigma = 0.0;
		for(i=0; i<size; i++) 
			{
			diff[i] = diff[i] - mean;
			sigma += diff[i]*diff[i];
			}
		sigma = nsigma*sqrt(sigma/size);

		/* determine outliers */
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

	free(x); free(y); free(diff);
	return err;
}
//-----------------------------------------------------------------------

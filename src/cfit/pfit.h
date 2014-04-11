#ifndef _P_FIT_H
#define _P_FIT_H

#include <stdio.h>

#define NORMALIZE(x,min,max) ( ((x*2.0)-min-max) / (max-min) )
#define DENORMALIZE(x,min,max) ( (x*(max-min) + max + min) / 2.0 )
#define id(n, m, N) (((m) * (N) + (n)))

void pfit_minmax(float x[], int N, float *min, float *max);
void pfit_normalize(float x[], int N, float min, float max);
void pfit_denormalize(float x[], int N, float min, float max);
float pfit_poly_eval(float x, float c[], int order);
float pw(float x, int n);
void pfit_qrgsr(int n, float c[], int N, float x[], float y[], float pfit_lambda);
void pfit_qrgsr_adaptive(int n, float c[], int N, float x[], float y[]);
void pfit_ldl(int n, float c[], int N, float x[], float y[], float pfit_lambda);
void pfit_ldl_adaptive(int n, float c[], int N, float x[], float y[]);
void pfit_lsq(int n, float c[], int N, float x[], float y[], float pfit_lambda);
void pfit_lsq_adaptive(int n, float c[], int N, float x[], float y[]);
int pfit_poly_fit(float X[], float Y[], int size, float nsigma, float C[], int order, int pfit_type, float pfit_lambda);
#endif

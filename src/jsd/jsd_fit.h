#ifndef _JSD_FIT_H
#define _JSD_FIT_H

#include <stdio.h>
/*
#define SQR(X) ((X)*(X))
*/
#define NORMALIZE(x,min,max) ( ((x*2.0)-min-max) / (max-min) )
#define DENORMALIZE(N,min,max) ( (N*(max-min) + max + min) / 2.0 )

void jsd_minmax(double A[], int n, double *p_min, double *p_max);
void jsd_normalize(double A[], int n, double min, double max);
void jsd_denormalize(double A[], int n, double min, double max);
float jsd_fpoly_eval(float x, float C[], int order);
double jsd_poly_eval(double x, double C[], int order);
double jsd_linear_eval(double x, double C[]);

/*
void jsd_iterative_fpoly_fit(float X[], float Y[], float sig[], int size, float C[], int order, float *chisq);
*/
int jsd_poly_fit(double X[], double Y[], int size, float nsigma, double C[], int order, double *chisq);
//void jsd_linear_fit(double X[], double Y[], int size, float nsigma, double C[], double *chisq);
void jsd_print_poly(FILE *file, double C[], int order);

#endif


#ifndef _NRFIT_H
#define _NRFIT_H

/* return zero on success, nonzero on error
*/
int svdfit(double x[], double y[], double sig[], int ndata, double a[], int ma,
	double **u, double **v, double w[], double *chisq,
	void (*funcs)(double, double [], int));

void svdvar(double **v, int ma, double w[], double **cvm);

void fpoly(double x, double p[], int np);

#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include "jsd_fit.h"
#include "nrfit.h"
#include "nrutil.h"

#define LOOP_COUNT_MAX 20

void jsd_minmax(double A[], int n, double *pmin, double *pmax)
{
	int i;
	double min, max;

	if (n <= 0) {
		*pmin = *pmax = 0.0;
		return;
	}

	min = max = A[0];
	for (i=1; i<n; i++) {
		double val = A[i];
		if (isless(val, min))
			min = val;
		else if (isgreater(val, max))
			max = val;
	}
	*pmin = min;
	*pmax = max;
}

void jsd_normalize(double A[], int n, double min, double max)
{
	int i;

	for (i=0; i<n; i++) {
		A[i] = NORMALIZE(A[i], min, max);
	}
}
void jsd_denormalize(double A[], int n, double min, double max)
{
	int i;

	for (i=0; i<n; i++) {
		A[i] = DENORMALIZE(A[i], min, max);
	}
}


/* order - the order of the polynomial
 * C - the coefficients (must be of size order+1) [0..order]
 * x - the value to evaluate for
 */
double jsd_poly_eval(double x, double C[], int order)
{
	int i;
	float val;

	val = C[order];
	for (i=order-1; i>=0; i--) {
		val = val * x + C[i];
	}

	return val;
}


double jsd_linear_eval(double x, double C[])
{
	return C[0] + C[1]*x;
}


/* order - the order of the polynomial
 * C - the coefficients (must be of size order+1) [0..order]
 * x - the value to evaluate for
 */
float jsd_fpoly_eval(float x, float C[], int order)
{
	int i;
	float val;

	val = C[order];
	for (i=order-1; i>=0; i--) {
		val = val * x + C[i];
	}

	return val;
}

/*
void jsd_linear_fit(double X[], double Y[], int size, float nsigma, double C[], double *chisq)
{
	double c0, c1, cov00, cov01, cov11;
	double sigma;
	int outlier;
	double mean;
	double *eval;
	int i;

	eval = malloc(sizeof(double) * size);

	do
	{
		if (size <= 1) {
			C[0] = Y[0];
			C[1] = 0;
			return;
		}

		gsl_fit_linear (X, 1, Y, 1, size, &c0, &c1, &cov00, &cov01, &cov11, chisq);
		if (isnan(c0+c1)) {
			C[0] = 0.0;
			C[1] = 0.0;
			return;
		}

		// do the evaluations
		mean = 0.0;
		for (i=0; i<size; i++) {
			eval[i] = c0+c1*X[i];
			mean += Y[i];
		}
		mean /= size;

		// calc sigma
		sigma = 0.0;
		for (i=0; i<size; i++) {
			sigma += SQR(Y[i]-mean);
		}
		sigma = sqrt(sigma/size);


		// determine outliers
		outlier = 0;
		for (i=0; i<size; i++) {
			if (fabs(mean-Y[i]) > nsigma*sigma) {
				outlier = 1;
				//printf("outlier=%g: (size=%i mean=%g sigma=%g nsigma=%g)\n",
				//	Y[i], size, mean, sigma, nsigma);
				size = size-1;
				Y[i] = Y[size];
				X[i] = X[size];
			}
		}
	} while (outlier && size > 1);

	C[0] = c0;
	C[1] = c1;

	free(eval);
}
*/
void jsd_print_poly(FILE *file, double C[], int order)
{
	int i;

	fprintf(file, "f%i(x) = ", order);
	for (i=0; i<order; i++) {
		fprintf(file, "%g*x**%i + ", C[i], i);
	}
	fprintf(file, "%g*x**%i\n", C[i], i);
}

int jsd_poly_fit(double X[], double Y[], int size, float nsigma, double C[], int order, double *chisq)
{
	int i;
	int terms = order+1;
	double *x,*y,*sig,*a,*w,**cvm,**u,**v;
	double sigma;
	int outlier;
	double mean;
	double *diff, *eval;
	int err = 0; //assume no error
	if (nsigma < 1.0) {
		printf("ERROR: nsigma less than 1.0\n");
		return -2;
	}
	if(size < order+1)
	{
		for (i=1;i<=terms;i++)
			C[i-1] = 0.0;
		return -1;
	}
	diff = malloc(sizeof(double) * size+1);
	eval = malloc(sizeof(double) * size+1);

	x=dvector(1,size);
	y=dvector(1,size);
	sig=dvector(1,size);
	a=dvector(1,terms);
	w=dvector(1,terms);
	cvm=dmatrix(1,terms,1,terms);
	u=dmatrix(1,size,1,terms);
	v=dmatrix(1,terms,1,terms);
		
	for (i=1;i<=size;i++) { //nr fortran style indexing
		x[i]=X[i-1];
		y[i]=Y[i-1];
		sig[i]=1.0;
	}
	
	do {
/*		for (i=1;i<=size;i++) { //nr fortran style indexing
			x[i]=X[i-1];
			y[i]=Y[i-1];
			sig[i]=1.0;
		}
*/
		err = svdfit(x,y,sig,size,a,terms,u,v,w,chisq,fpoly);
		if (err || !isfinite(a[1])) {
			//printf("ERROR: svdfit curve fit failed! %d\n",size);
			for (i=1;i<=terms;i++)
			C[i-1] = 0.0;
			err = -1;
			break;
		}
		for (i=1;i<=terms;i++)
			C[i-1] = a[i];

		// do the evaluations 
		mean = 0.0;
//		for (i=0; i<size; i++) { //normal indexing
		for (i=1; i<=size; i++) { //nr indexing
			eval[i] = jsd_poly_eval(x[i], C, order);
			diff[i] = eval[i] - y[i];
			mean += diff[i];
		}
		mean /= size;

		// calc sigma 
		sigma = 0.0;
//<<<<<<< jsd_fit.c
//		for (i=0; i<size; i++) {
//		for (i=1; i<=size; i++) {
//			sigma += (diff[i]-mean)*(diff[i]-mean);
//=======
		for (i=0; i<size; i++) {
			sigma += SQR(diff[i]-mean);
//>>>>>>> 1.4
		}
		sigma = sqrt(sigma/size);

		// determine outliers 
		outlier = 0;
		for (i=1; i<=size; i++) {
		//for (i=0; i<size; i++) {
			if (fabs(diff[i] - mean) > nsigma*sigma) {
				outlier = 1;
				//printf("outlier=%g: (size=%i mean=%g sigma=%g nsigma=%g)\n",
				//	Y[i], size, mean, sigma, nsigma);
				//Y[i] = eval[i];
				//size = size-1;
				x[i] = x[size];
				y[i] = y[size];
				size = size-1;
			}
		}
	} while (outlier && size > order+1);
	free(eval);
	free(diff);
	free_dmatrix(v,1,terms,1,terms);
	free_dmatrix(u,1,size,1,terms);
	free_dmatrix(cvm,1,terms,1,terms);
	free_dvector(w,1,terms);
	free_dvector(a,1,terms);
	free_dvector(sig,1,size);
	free_dvector(y,1,size);
	free_dvector(x,1,size);
	return err;
}


/*int jsd_poly_fit(double X[], double Y[], int size, float nsigma,double C[], int order, double *chisq)
{
	int i, j;
	gsl_matrix *x, *cov;
	gsl_vector *y, *c;
	int terms = order + 1;

	//allocate memory
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (size, terms);
	x = gsl_matrix_alloc (size, terms);
	y = gsl_vector_alloc (size);
	c = gsl_vector_alloc (terms);
	cov = gsl_matrix_alloc (terms, terms);

	//prepare the fit start conditions
	for (i = 0; i < size; i++)
	{
		double mv = 1.0;
		for (j=0; j<terms; j++) {
			gsl_matrix_set (x, i, j, mv);
			mv *= X[i];
		}

		gsl_vector_set (y, i, Y[i]);
	}

	//perform the fit
	gsl_multifit_linear (x, y, c, cov, chisq, work);

	//store the results in the output array
	for (i=0; i<terms; i++) {
		C[i] = gsl_vector_get (c, i);
	}

	//cleanup memory
	gsl_multifit_linear_free (work);
	gsl_matrix_free(x);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
	return 1;
}
*/


/* a - array of coefficients 1..ma
 * X is normalized on return of this function
 */
/*
void jsd_iterative_fpoly_fit(float X[], float Y[], float sig[], int size, float C[], int order, float *chisq)
{
    int i;
    int outlierFound;
    double sigma;
    float *eval;
	int nsigma = 3;  //TODO: parameterize this or somthing
	int count;

	int ndata = size;
	int ma = order+1;
	float **U, **V, *W; //workspaces for svdfit

	U = matrix(1, ndata, 1, ma);
	V = matrix(1, ma, 1, ma);
	W = vector(1, ma);
    eval = vector(1, size);

	for (i=0; i<size; i++) {
		sig[i] = 1.0;
	}

	count = 0;
    do {
        sigma = 0;
        outlierFound = 0;
        //do the fit
		jsd_normalize(X, size);
		nr_svdfit(X-1, Y-1, sig-1, ndata, C-1, ma, U, V, W, chisq, fpoly);

        //calculate the difference between the function evaluation and Y
        //and also calculate sigma using the function evaluation as the mean
        for (i=0; i<size; i++) {
            eval[i] = jsd_fpoly_eval(X[i], C, order);
            sigma += SQR(Y[i] - eval[i]);
        }
        sigma = sqrt(sigma/(size-1));

        //reject outliers by making the value equal to the function
        for (i=0; i<size; i++) {
            if (fabs(eval[i] - Y[i]) > nsigma*sigma*sig[i]) {
                outlierFound = 1;
				sig[i] = INFINITY;
				//printf("outlier at (%f, %f)\n",  X[i], Y[i]);
            }
        }

		//fprintf(stdout, "%c(x) = ", 'f'+count);
		//for (i=0; i<=order; i++) fprintf(stdout, "%+.6g*x**%i ", C[i], i);
		//fprintf(stdout, "\n");

        count++;
        //repeat until there are no more outliers

    } while (outlierFound && count < LOOP_COUNT_MAX);

    printf("curve fit iterations: %i\n", count);

    free_vector(eval, 1, size);
	free_vector(W, 1, ma);
	free_matrix(V, 1, ma, 1, ma);
	free_matrix(U, 1, ndata, 1, ma);
}

*/

#include <complex.h>

typedef struct
{
	float x;
	float y;
	complex double z;
}dftdata;

void dft2d(dftdata input[], int numpoints, dftdata output[]);
void idft2d(dftdata input[], int numpoints, dftdata output[]);
void zernike(int n, int m, float x[], float y[], int numpoints, float z[],float radius);
void zernike_pt(int n, int m, float x, float y, float *z,float radius);
float eval_dft(dftdata input[], int numpoints, float x, float y);
//int factorial(int n);

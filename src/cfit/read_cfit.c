/*
  Least-Squares 2d Chebyshev Surface Approximation.
  Copyright (C) 2011, M. Andrecut

  Example.
  Construct the least-squares surface of degree n
  that fits the N data points: 
  (x(0), y(0), z(0)),...,(x(N-1), y(N-1), z(N-1)).

  Use the Chebyshev Approximation Library
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "chebyshev.h"


int main(int argc, char **argv)
{
int i,j,k,n;

float tmp_dRA, tmp_dDEC, tmp_AST, tmp_I, tmp_Q, tmp_U, tmp_V;

float xmin, xmax, ymin, ymax;

FILE *file;

int beam = atoi(argv[1]); /* beam */
int chan = atoi(argv[2]); /* channel */

char filename[64];
sprintf(filename,"beam%d_fit%04i.dat",beam,chan);
printf("%s\n",filename);
fflush(stdout);
file = fopen(filename,"r");

fscanf(file,"%d %f %f %f %f\n",&n,&xmin,&xmax,&ymin,&ymax);
printf("%d %2.9f %2.9f %2.9f %2.9f\n",n,xmin,xmax,ymin,ymax);
printf("\n");
fflush(stdout);
float *c = (float*)malloc((n+1)*(n+1) * sizeof(float)); // the coefficients of the Chebyshev expansion

for(i=0; i<=n; i++)
{ 
	for(j=0; j<=n; j++)
	{
		fscanf(file,"%f ",&c[i*(n+1)+j]);
//		printf("%2.9f ",c[i*(n+1)+j]);
	}
//	printf("\n");
}
//printf("\n");

for(i=0; i<=n; i++)
{ 
	for(j=0; j<=n; j++)
	{
		fscanf(file,"%f ",&c[i*(n+1)+j]);
//		printf("%2.9f ",c[i*(n+1)+j]);
	}
//	printf("\n");
}
//printf("\n");

for(i=0; i<=n; i++)
{ 
	for(j=0; j<=n; j++)
	{
		fscanf(file,"%f ",&c[i*(n+1)+j]);
//		printf("%2.9f ",c[i*(n+1)+j]);
	}
//	printf("\n");
}
//printf("\n");

for(i=0; i<=n; i++)
{ 
	for(j=0; j<=n; j++)
	{
		fscanf(file,"%f ",&c[i*(n+1)+j]);
//		printf("%2.9f ",c[i*(n+1)+j]);
	}
//	printf("\n");
}
//printf("\n");


fclose(file);


printf("poly_order = %d\n", n);
sprintf(filename,"beam%d_model%04i.dat",beam,chan);
sprintf(filename,"beam5_model%04i.dat",chan);

file=fopen(filename, "r");

int N; fscanf(file, "%d", &N); printf("N = %d\n", N);

float *x = (float*)malloc(N * sizeof(float)); 
float *y = (float*)malloc(N * sizeof(float)); 
float *z = (float*)malloc(N * sizeof(float));
float *f = (float*)malloc(N * sizeof(float));
int count = 0;
for(i=0; i<N; i++)
	{
	fscanf(file, "%f %f %f %f %f %f %f\n", &tmp_dRA, &tmp_dDEC, &tmp_AST, &tmp_I, &tmp_Q, &tmp_U, &tmp_V); 
	if(fabs(tmp_dRA) < 0.14 && fabs(tmp_dDEC)< 0.14)
	{
		x[i] = tmp_dRA;
		y[i] = tmp_dDEC;
	//z[i] = tmp_I;    // choose I, Q, U, V !!!
	//z[i] = tmp_Q;    // choose I, Q, U, V !!!
	//z[i] = tmp_U;    // choose I, Q, U, V !!!
		z[i] = tmp_V;    // choose I, Q, U, V !!!
		count++;
	}
	}
fclose(file);
N = count;


chebyshev_minmax(x, N, &xmin, &xmax);
chebyshev_minmax(y, N, &ymin, &ymax);

chebyshev_normalize(x, N, xmin, xmax);
chebyshev_normalize(y, N, ymin, ymax);

for(i=0; i<N; i++) 
	{
	f[i] = chebyshev_eval_surface(x[i], y[i], c, n); // Chebyshev surface evaluation at points (x, y)
	}

chebyshev_denormalize(x, N, xmin, xmax);
chebyshev_denormalize(y, N, ymin, ymax);

		
//printf("The coefficients of the least squares surface are:\n");
//for(k=0; k<=(n+1)*(n+1)-1; k++) {printf("c[%d]=%f\n", k, c[k]);}


file=fopen("results1.dat", "w");
for(i=0; i<N; i++)
	{
	fprintf(file, "%f\t%f\t%f\t%f\n", x[i], y[i], z[i], f[i]);
	}
fclose(file);

FILE *pipe = popen("gnuplot -persist","w");
fprintf(pipe, "splot 'results1.dat' u 1:2:3 title 'original', 'results1.dat' u 1:2:4 with points title 'fit'\n");
fclose(pipe);
free(x); 
free(y); 
free(z); 
free(f); 
free(c);

return EXIT_SUCCESS;
}


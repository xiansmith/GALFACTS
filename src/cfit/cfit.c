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
int i,j,k;

float tmp_dRA, tmp_dDEC, tmp_AST, tmp_I, tmp_Q, tmp_U, tmp_V;

int beam = atoi(argv[1]); /* beam */
int chan = atoi(argv[2]); /* channel */
int n = atoi(argv[3]); /* degree of Chebyshev surface */
printf("poly_order = %d\n", n);

int pfit_type = 0; /* keep it 0 */

float pfit_lambda = 0; /* keep it 0 */

FILE *file;
char filename[64];
sprintf(filename,"beam%d_model%04i.dat",beam,chan);
printf("%s\n",filename);
file=fopen(filename, "r");

int N; fscanf(file, "%d", &N); printf("N = %d\n", N);

float *x = (float*)malloc(N * sizeof(float)); 
float *y = (float*)malloc(N * sizeof(float)); 
float *I = (float*)malloc(N * sizeof(float));
float *Q = (float*)malloc(N * sizeof(float));
float *U = (float*)malloc(N * sizeof(float));
float *V = (float*)malloc(N * sizeof(float));
float *f = (float*)malloc(N * sizeof(float));

for(i=0; i<N; i++)
	{
	fscanf(file, "%f %f %f %f %f %f %f\n", &tmp_dRA, &tmp_dDEC, &tmp_AST, &tmp_I, &tmp_Q, &tmp_U, &tmp_V); 
	x[i] = tmp_dRA;
	y[i] = tmp_dDEC;
	I[i] = tmp_I;    // choose I, Q, U, V !!!
	Q[i] = tmp_Q;    // choose I, Q, U, V !!!
	U[i] = tmp_U;    // choose I, Q, U, V !!!
	V[i] = tmp_V;    // choose I, Q, U, V !!!
	}
fclose(file);

/* The fitting starts here */

float xmin, xmax, ymin, ymax;

chebyshev_minmax(x, N, &xmin, &xmax);
chebyshev_minmax(y, N, &ymin, &ymax);

chebyshev_normalize(x, N, xmin, xmax);
chebyshev_normalize(y, N, ymin, ymax);

clock_t start=clock(); 

float *c = (float*)malloc((n+1)*(n+1) * sizeof(float)); // the coefficients of the Chebyshev expansion

printf("Computation time: %f sec\n", ((float)clock()-start)/CLOCKS_PER_SEC);    

sprintf(filename,"beam%d_fit%04i.dat",beam,chan);
file = fopen(filename,"w");
fprintf(file,"%d %2.9f %2.9f %2.9f %2.9f\n",n,xmin,xmax,ymin,ymax);
fprintf(file,"\n");

chebyshev_qrgsr_surface(n, c, N, x, y, I, pfit_lambda); // Chebyshev surface fit, i.e. finds coefficients c

for(i=0; i<=n; i++)
{ 
	for(j=0; j<=n; j++)
	{
		fprintf(file,"%2.9f ",c[i*(n+1)+j]);
	}
	fprintf(file,"\n");
}
fprintf(file,"\n");

chebyshev_qrgsr_surface(n, c, N, x, y, Q, pfit_lambda); // Chebyshev surface fit, i.e. finds coefficients c

for(i=0; i<=n; i++)
{ 
	for(j=0; j<=n; j++)
	{
		fprintf(file,"%2.9f ",c[i*(n+1)+j]);
	}
	fprintf(file,"\n");
}
fprintf(file,"\n");

chebyshev_qrgsr_surface(n, c, N, x, y, U, pfit_lambda); // Chebyshev surface fit, i.e. finds coefficients c

for(i=0; i<=n; i++)
{ 
	for(j=0; j<=n; j++)
	{
		fprintf(file,"%2.9f ",c[i*(n+1)+j]);
	}
	fprintf(file,"\n");
}
fprintf(file,"\n");

chebyshev_qrgsr_surface(n, c, N, x, y, V, pfit_lambda); // Chebyshev surface fit, i.e. finds coefficients c

for(i=0; i<=n; i++)
{ 
	for(j=0; j<=n; j++)
	{
		fprintf(file,"%2.9f ",c[i*(n+1)+j]);
	}
	fprintf(file,"\n");
}


fclose(file);


/*file=fopen("results.dat", "w");
float xx=-0.9,yy=-0.9,ff;
for(i=0; i<=100; i++)
{ 
	xx = -0.9;
	for(j=0; j<=100; j++)
	{
		xx = xx+0.018;
		ff = chebyshev_eval_surface(xx, yy, c, n); // Chebyshev surface evaluation at points (x, y)
		float xxx,yyy;
		xxx = CDENORMALIZE(xx,xmin,xmax);
		yyy = CDENORMALIZE(yy,ymin,ymax);
		fprintf(file, "%f\t%f\t%f\n", xxx, yyy, ff);
	}
	yy = yy+0.018;
}
fclose(file);
*/

for(i=0; i<N; i++) 
	{
	f[i] = chebyshev_eval_surface(x[i], y[i], c, n); // Chebyshev surface evaluation at points (x, y)
	}

chebyshev_denormalize(x, N, xmin, xmax);
chebyshev_denormalize(y, N, ymin, ymax);

/* The fitting ends here */
		
//printf("The coefficients of the least squares surface are:\n");
//for(k=0; k<=(n+1)*(n+1)-1; k++) {printf("c[%d]=%f\n", k, c[k]);}

/* Plot results */

file=fopen("results.dat", "w");
for(i=0; i<N; i++)
	{
	fprintf(file, "%f\t%f\t%f\t%f\n", x[i], y[i], V[i], f[i]);
	//fprintf(file, "%f\t%f\t%f\n", x[i], y[i], f[i]);
	}
fclose(file);

/*FILE *pipe = popen("gnuplot -persist","w");
fprintf(pipe, "splot 'results.dat' u 1:2:3 title 'original', 'results.dat' u 1:2:4 with points title 'fit'\n");
fclose(pipe);*/
free(x); 
free(y); 
free(I); 
free(Q); 
free(U); 
free(V); 
free(f); 
free(c);
return EXIT_SUCCESS;
}

#include <stdlib.h>
#include <stdio.h>

#include "jsd_fit.h"
#include "jsd_futil.h"

#define BUF_SIZE 128

int main(int argc, char *argv[])
{
	int i;
	int count;
	FILE *datafile;
	double *X, *Y, *Z;
	double C[100];
	double chisq;
	double nsigma = 3;
	char *filename;
	int order;
	double min, max;
	char buf[BUF_SIZE];	
	
	if (argc != 3) {
		printf("Usage: %s <file> <order>\n", argv[0]);
		return EXIT_FAILURE;
	}
	filename = argv[1];
	order = atoi(argv[2]);

	datafile = fopen(filename, "r");
	if (!datafile) {
		printf("ERROR: unable to open '%s'\n", filename);
		return EXIT_FAILURE;
	}
	
	count = jsd_line_count(datafile);
	X = malloc(count * sizeof(double));
	Y = malloc(count * sizeof(double));
	Z = malloc(count * sizeof(double));
	for (i=0; i<count; i++) {
		fgets(buf, BUF_SIZE, datafile);
		if (buf[0] == '#') continue;
		sscanf(buf, "%lf %lf\n", &X[i], &Y[i]);
	}
	fclose(datafile);

	jsd_minmax(X, count, &min, &max);
	printf("normalizing: %f %f\n", min, max);
	printf("first/last %f %f\n", X[0], X[count-1]);
	jsd_normalize(X, count, min, max);

	{
		double m, n;
		jsd_minmax(X, count, &m, &n);
		printf("normalized: %f %f\n", m, n);
		printf("first/last %f %f\n", X[0], X[count-1]);
	}

	jsd_linear_fit(X, Y, count, nsigma, C, &chisq);
	printf("f1(x) = %f + %f*x\n", C[0], C[1]);

	for (i=0; i<=order; i++) {
		jsd_poly_fit(X, Y, count, C, i, &chisq);
		jsd_print_poly(stdout, C, i);
	}
	
	datafile = fopen("normal.dat", "w");
	for (i=0; i<count; i++) {
		fprintf(datafile, "%lf %lf\n",  X[i], Y[i]);
	}
	fclose(datafile);



	jsd_denormalize(X, count, min, max);
	jsd_minmax(X, count, &min, &max);
	printf("denormalized: %f %f\n", min, max);
	printf("first/last %f %f\n", X[0], X[count-1]);

	free(X);
	free(Y);
	free(Z);

	return EXIT_SUCCESS;
}

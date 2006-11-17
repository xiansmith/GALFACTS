#include "calibrate.h"
#include "rfi.h"
#include <stdlib.h>
#include "jsd/jsd_fit.h"



//pre-compute the raw cal for every channel
//void compute_raw_cal(SpecRecord dataset[], int size)
void compute_raw_cal(SpecRecord dataset[],int lowchan, int highchan, int size)
{
	int n, i;
	for (n=0; n<size; n++)
	{
		//for (i=0; i<MAX_CHANNELS; i++)
		for (i=lowchan; i<=highchan; i++) 
		{
			dataset[n].cal.xx[i] = dataset[n].calon.xx[i] - dataset[n].caloff.xx[i];
			dataset[n].cal.yy[i] = dataset[n].calon.yy[i] - dataset[n].caloff.yy[i];
			dataset[n].cal.xy[i] = dataset[n].calon.xy[i] - dataset[n].caloff.xy[i];
			dataset[n].cal.yx[i] = dataset[n].calon.yx[i] - dataset[n].caloff.yx[i];
		}
	}
}

void print_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI)
{	
	int n, i;
	FILE * file;

	file = fopen("cal.dat", "w");

	fprintf(file, "#AST chan, calXX, calYY, calXY, calYX\n");

	for (n=0; n<size; n++)
	{
		for (i=lowchan; i<highchan; i++)
		{
			if (dataset[n].flagBAD) continue;
			if (ignoreRFI || dataset[n].flagRFI[i] == RFI_NONE) {
				fprintf(file, "%8.2f %4i %8.6f %8.6f %8.6f %8.6f\n",
						dataset[n].AST, i, 
						dataset[n].cal.xx[i], dataset[n].cal.yy[i], 
						dataset[n].cal.xy[i], dataset[n].cal.yx[i]);
			}
		}
	}
	fclose(file);
}


void linear_fit_cal(SpecRecord dataset[], int size, int lowchan, int highchan, int ignoreRFI)
{
	int n, chan;
	double C[2]; //assuming linear fit 
	float nsigma = 5.0;
	double sumsq;
	double *Xast, *Yxx, *Yyy, *Yxy, *Yyx;
	FILE * eqfile;
	FILE * chifile;
	size_t count;

	eqfile = fopen("caleq.dat", "w");
	fprintf(eqfile, "#chan XXc0 XXc1 YYc0 YYc1 XYc0 XYc1 YXc0 YXc1\n");

	chifile = fopen("calchi.dat", "w");
	fprintf(chifile, "#chan XX YY XY YX\n");

	//y is the cal value
	//x is the time steps
	Xast = (double*) malloc(sizeof(double) * size);
	Yxx = (double*) malloc(sizeof(double) * size);
	Yyy = (double*) malloc(sizeof(double) * size);
	Yxy = (double*) malloc(sizeof(double) * size);
	Yyx = (double*) malloc(sizeof(double) * size);

	for (chan=lowchan; chan<=highchan; chan++) 
//ssg	for (chan=0; chan<=MAX_CHANNELS; chan++) 
	{
		fprintf(chifile, "%i ", chan);

		count = 0;
		for (n=0; n<size; n++) {
			if (dataset[n].flagBAD) continue;
			if (ignoreRFI || dataset[n].flagRFI[chan] == RFI_NONE) {
				Xast[count] = dataset[n].AST;
				Yxx[count] = dataset[n].cal.xx[chan];
				Yyy[count] = dataset[n].cal.yy[chan];
				Yxy[count] = dataset[n].cal.xy[chan];
				Yyx[count] = dataset[n].cal.yx[chan];
				count++;
			}
		}

		if (count == 0) {
			//NaN values will result in the dataset cal
		}

		/* normalize the x values for the curve fit */
		{
			double min, max;
			jsd_minmax(Xast, count, &min, &max);
			jsd_normalize(Xast, count, min, max);
		}

		//XX
		jsd_linear_fit(Xast, Yxx, count, nsigma, C, &sumsq);
		fprintf (eqfile, "XX %i %g %g\n", chan, C[0], C[1]);
		fprintf(chifile, "%g ", sumsq);
		for (n=0; n<size; n++) //eval function
			dataset[n].cal.xx[chan] = jsd_linear_eval(Xast[n], C);

		//YY
		jsd_linear_fit(Xast, Yyy, count, nsigma, C, &sumsq);
		fprintf (eqfile, "YY %i %g %g\n", chan, C[0], C[1]);
		fprintf(chifile, "%g ", sumsq);
		for (n=0; n<size; n++) //eval function
			dataset[n].cal.yy[chan] = jsd_linear_eval(Xast[n], C);

		//XY
		jsd_linear_fit(Xast, Yxy, count, nsigma, C, &sumsq);
		fprintf (eqfile, "XY %i %g %g\n", chan, C[0], C[1]);
		fprintf(chifile, "%g ", sumsq);
		for (n=0; n<size; n++) //eval function
			dataset[n].cal.xy[chan] = jsd_linear_eval(Xast[n], C);

		//YX
		jsd_linear_fit(Xast, Yyx, count, nsigma, C, &sumsq);
		fprintf (eqfile, "YX %i %g %g\n", chan, C[0], C[1]);
		fprintf(chifile, "%g ", sumsq);
		for (n=0; n<size; n++) //eval function
			dataset[n].cal.yx[chan] = jsd_linear_eval(Xast[n], C);

		fprintf(chifile, "\n");

	}

	fclose(eqfile);
	fclose(chifile);
	free(Xast);
	free(Yxx);
	free(Yyy);
	free(Yxy);
	free(Yyx);
}



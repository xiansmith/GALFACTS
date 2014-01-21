#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "decdep.h"
#include "chebyshev.h"
//----------------------------------------------------------------------------------------------------------------------------------------
static void day_dec_dependence(FluxWappData * wappdata, int day, int beam, int order, int chan, float *cIc, float *cQc, float *cUc, float *cVc, int avg)
{
	int Hchan = 2752; // hard wired!!!
	int r, N, R, n, navg;
	if (avg == 0)
		navg = 1;
	else navg = avg;
	float Hfreq = 1420.4057, Cfreq = 1450, df = 0.042, freq = Cfreq - ((float) (chan + 1 + navg - MAX_CHANNELS / 2.0)) * df;
	float min, max, DEC, nsigma = 2.5;
	float cI[order + 1], cQ[order + 1], cU[order + 1], cV[order + 1];
	char decrmfilename[50];
	//sprintf(decrmfilename, "decrm/decremoval%d.dat", day);

	FluxDayData * daydata = &wappdata->daydata[day];
	sprintf(decrmfilename, "%s/beam%d/decremoval.dat", daydata->mjd,beam);

	R = daydata->numRecords;

	float *x = (float*) malloc(R * sizeof(float));
	float *yI = (float*) malloc(R * sizeof(float));
	float *yQ = (float*) malloc(R * sizeof(float));
	float *yU = (float*) malloc(R * sizeof(float));
	float *yV = (float*) malloc(R * sizeof(float));

// not doing average image, read coefficients from a file
// and apply them
	if (chan) {
		//printf(" from file\n");
		FILE *decfile = fopen(decrmfilename, "r");
		char line[500];
		int i = 0, num = 0;
		char *result;

		if (decfile != NULL) {

			// get coefficients for I from file
			if (fgets(line, 500, decfile) == NULL) {
				printf("in decdep.c, fgets read error line is %s\n", line);
				exit(1);
			}

			i = 0;
			result = strtok(line, " ");
			while (result != NULL) {
				sscanf(result, "%f", &cI[i]);
				result = strtok(NULL, " ");
				i++;
			}

			if (i != (order + 1)) {
				printf("ERROR: read %d dec removal coeffs when expected %d\n", i, (order + 1));
			}


			// now Q
			if (fgets(line, 500, decfile) == NULL) {
				printf("in decdep.c, fgets read error");
				exit(1);
			}

			i = 0;
			result = strtok(line, " ");
			while (result != NULL) {
				sscanf(result, "%f", &cQ[i]);
				result = strtok(NULL, " ");
				i++;
			}

			if (i != (order + 1)) {
				printf("ERROR: read %d dec removal coeffs when expected %d\n", i, (order + 1));
			}

			// U
			if (fgets(line, 500, decfile) == NULL) {
				printf("in decdep.c, fgets read error line is %s\n", line);
				exit(1);
			}

			i = 0;
			result = strtok(line, " ");
			while (result != NULL) {
				sscanf(result, "%f", &cU[i]);
				result = strtok(NULL, " ");
				i++;
			}

			if (i != (order + 1)) {
				printf("ERROR: read %d dec removal coeffs when expected %d\n", i, (order + 1));
			}

			// V
			if (fgets(line, 500, decfile) == NULL) {
				printf("in decdep.c, fgets read error line is %s\n", line);
			}

			i = 0;
			result = strtok(line, " ");
			while (result != NULL) {
				sscanf(result, "%f", &cV[i]);
				result = strtok(NULL, " ");
				i++;
			}

			if (i != (order + 1)) {
				printf("ERROR: read %d dec removal coeffs when expected %d\n", i, (order + 1));
			}

			fclose(decfile);
		}
		else {
			printf("ERROR: got null in decdep trying to get the coefficients from %s\n", decrmfilename);
		}

		N = 0;
		for (r = 0; r < R; r++) {
			if (isfinite(daydata->records[r].stokes.I) && isfinite(daydata->records[r].stokes.Q) && isfinite(daydata->records[r].stokes.U)
					&& isfinite(daydata->records[r].stokes.V)) {
				x[N] = daydata->records[r].DEC;
				yI[N] = daydata->records[r].stokes.I;
				yQ[N] = daydata->records[r].stokes.Q;
				yU[N] = daydata->records[r].stokes.U;
				yV[N] = daydata->records[r].stokes.V;
				N++;
			}
		}

		chebyshev_minmax(x, N, &min, &max);
		chebyshev_normalize(x, N, min, max);

		if (freq > Hfreq - 0.5 && freq < Hfreq + 0.5) {
			for (n = 0; n < order + 1; n++) {
				cI[n] = cIc[n + day * (order + 1)];
				cQ[n] = cQc[n + day * (order + 1)];
				cU[n] = cUc[n + day * (order + 1)];
				cV[n] = cVc[n + day * (order + 1)];
			}
		}
		else {
			chebyshev_fit_dec(x, yI, N, nsigma, cI, order);
			chebyshev_fit_dec(x, yQ, N, nsigma, cQ, order);
			chebyshev_fit_dec(x, yU, N, nsigma, cU, order);
			chebyshev_fit_dec(x, yV, N, nsigma, cV, order);
			for (n = 0; n < order + 1; n++) {
				cIc[n + day * (order + 1)] = cI[n];
				cQc[n + day * (order + 1)] = cQ[n];
				cUc[n + day * (order + 1)] = cU[n];
				cVc[n + day * (order + 1)] = cV[n];
			}
		}

		// now apply it
		for (r = 0; r < R; r++) {
			if (isfinite(daydata->records[r].stokes.I) && isfinite(daydata->records[r].stokes.Q) && isfinite(daydata->records[r].stokes.U)
					&& isfinite(daydata->records[r].stokes.V)) {
				DEC = CNORMALIZE(daydata->records[r].DEC, min, max);
	            // XXX temporary magic number    
				daydata->records[r].stokes.I = daydata->records[r].stokes.I  - chebyshev_eval(DEC, cI, order) + 0.75;
				//daydata->records[r].stokes.I -= chebyshev_eval(DEC, cI, order);
				daydata->records[r].stokes.Q -= chebyshev_eval(DEC, cQ, order);
				daydata->records[r].stokes.U -= chebyshev_eval(DEC, cU, order);
				daydata->records[r].stokes.V -= chebyshev_eval(DEC, cV, order);
			}
		}

	}
	else {
		// we are processing average image, create coefficients and write to file
		//printf(" and write to file\n");

		N = 0;
		for (r = 0; r < R; r++) {
			if (isfinite(daydata->records[r].stokes.I) && isfinite(daydata->records[r].stokes.Q) && isfinite(daydata->records[r].stokes.U)
					&& isfinite(daydata->records[r].stokes.V)) {
				x[N] = daydata->records[r].DEC;
				yI[N] = daydata->records[r].stokes.I;
				yQ[N] = daydata->records[r].stokes.Q;
				yU[N] = daydata->records[r].stokes.U;
				yV[N] = daydata->records[r].stokes.V;
				N++;
			}
		}

		chebyshev_minmax(x, N, &min, &max);
		chebyshev_normalize(x, N, min, max);

		if (freq > Hfreq - 0.5 && freq < Hfreq + 0.5) {
			for (n = 0; n < order + 1; n++) {
				cI[n] = cIc[n + day * (order + 1)];
				cQ[n] = cQc[n + day * (order + 1)];
				cU[n] = cUc[n + day * (order + 1)];
				cV[n] = cVc[n + day * (order + 1)];
			}
		}
		else {
			chebyshev_fit_dec(x, yI, N, nsigma, cI, order);
			chebyshev_fit_dec(x, yQ, N, nsigma, cQ, order);
			chebyshev_fit_dec(x, yU, N, nsigma, cU, order);
			chebyshev_fit_dec(x, yV, N, nsigma, cV, order);
			for (n = 0; n < order + 1; n++) {
				cIc[n + day * (order + 1)] = cI[n];
				cQc[n + day * (order + 1)] = cQ[n];
				cUc[n + day * (order + 1)] = cU[n];
				cVc[n + day * (order + 1)] = cV[n];
			}
		}

		for (r = 0; r < R; r++) {
			if (isfinite(daydata->records[r].stokes.I) && isfinite(daydata->records[r].stokes.Q) && isfinite(daydata->records[r].stokes.U)
					&& isfinite(daydata->records[r].stokes.V)) {
				DEC = CNORMALIZE(daydata->records[r].DEC, min, max);
                // XXX temporary magic number 
				daydata->records[r].stokes.I = daydata->records[r].stokes.I  - chebyshev_eval(DEC, cI, order) + 0.75;
				//daydata->records[r].stokes.I -= chebyshev_eval(DEC, cI, order);
				daydata->records[r].stokes.Q -= chebyshev_eval(DEC, cQ, order);
				daydata->records[r].stokes.U -= chebyshev_eval(DEC, cU, order);
				daydata->records[r].stokes.V -= chebyshev_eval(DEC, cV, order);
			}
		}

		// finished applying correction, write the coefficients to file
		//mkdir("decrm", 0777);
		FILE *decfile = fopen(decrmfilename, "w"); // write files to decrm/decremovalN.dat
		if (decfile != NULL) {
			int i = 0;

			for (i = 0; i < order + 1; i++) {
				fprintf(decfile, "%.8f", cI[i]);
				if (i != order) fprintf(decfile, " ");
			}
			fprintf(decfile, "\n");

			for (i = 0; i < order + 1; i++) {
				fprintf(decfile, "%.8f", cQ[i]);
				if (i != order) fprintf(decfile, " ");
			}
			fprintf(decfile, "\n");

			for (i = 0; i < order + 1; i++) {
				fprintf(decfile, "%.8f", cU[i]);
				if (i != order) fprintf(decfile, " ");
			}
			fprintf(decfile, "\n");

			for (i = 0; i < order + 1; i++) {
				fprintf(decfile, "%.8f", cV[i]);
				if (i != order) fprintf(decfile, " ");
			}
			fprintf(decfile, "\n");

			fclose(decfile);

		}
		else {
			// can not recover
			printf("couldn't write decremoval.dat file\n");
			exit(1);
		}

	} // end if (avg == 0)

	free(x);
	free(yI);
	free(yQ);
	free(yU);
	free(yV);
}
//-------------------------------------------------------------------------------
void beam_gain_calibration(FluxWappData * wappdata)
{
int r, day;
for(day=0; day<wappdata->numDays; day++) 
	{
	FluxDayData * daydata = &wappdata->daydata[day];
	for(r=0; r<daydata->numRecords; r++)
		{
		if(isfinite(daydata->records[r].stokes.I))
			{
			switch(day%7)
				{
				case 0:	daydata->records[r].stokes.I /= ( 0.00122530 * daydata->records[r].DEC + 0.97985562); break;
				case 1:	daydata->records[r].stokes.I /= ( 0.00116780 * daydata->records[r].DEC + 0.75177204); break;
				case 2: daydata->records[r].stokes.I /= (-0.00025934 * daydata->records[r].DEC + 0.79175157); break; 
				case 3: daydata->records[r].stokes.I /= ( 0.00036153 * daydata->records[r].DEC + 0.75028658); break; 
				case 4: daydata->records[r].stokes.I /= (-0.00090732 * daydata->records[r].DEC + 0.74814666); break; 
				case 5: daydata->records[r].stokes.I /= (-0.00602525 * daydata->records[r].DEC + 0.96467256); break; 
				case 6: daydata->records[r].stokes.I /= ( 0.00116780 * daydata->records[r].DEC + 0.75177204); break; // like beam 1
				default: break;
				}
			}
		}
	}	
}
//-------------------------------------------------------------------------------
void beam_gain_calibration_table(FluxWappData * wappdata, int cal_low, int cal_high, float cal_table[][7], int chan)
{
int r, day;
/*
if(chan < cal_low || chan > cal_high)
	{
	for(day=0; day<wappdata->numDays; day++) 
		{
		FluxDayData * daydata = &wappdata->daydata[day];
		for(r=0; r<daydata->numRecords; r++)
			{
			if(isfinite(daydata->records[r].stokes.I))
				{
				switch(day%7)
					{
					case 0:	daydata->records[r].stokes.I /= ( cal_table[0][0] * daydata->records[r].DEC + cal_table[1][0]); break;
					case 1:	daydata->records[r].stokes.I /= ( cal_table[0][1] * daydata->records[r].DEC + cal_table[1][1]); break;
					case 2: daydata->records[r].stokes.I /= ( cal_table[0][2] * daydata->records[r].DEC + cal_table[1][2]); break;
					case 3: daydata->records[r].stokes.I /= ( cal_table[0][3] * daydata->records[r].DEC + cal_table[1][3]); break;
					case 4: daydata->records[r].stokes.I /= ( cal_table[0][4] * daydata->records[r].DEC + cal_table[1][4]); break;
					case 5: daydata->records[r].stokes.I /= ( cal_table[0][5] * daydata->records[r].DEC + cal_table[1][5]); break;
					case 6: daydata->records[r].stokes.I /= ( cal_table[0][6] * daydata->records[r].DEC + cal_table[1][6]); break; 
					default: break;
					}
				}
			}
		}	
	}
	else
		{
		for(day=0; day<wappdata->numDays; day++) 
			{
			FluxDayData * daydata = &wappdata->daydata[day];
			for(r=0; r<daydata->numRecords; r++)
				{
				if(isfinite(daydata->records[r].stokes.I))
					{
					switch(day%7)
						{
						case 0:	daydata->records[r].stokes.I /= ( cal_table[0][0] * daydata->records[r].DEC + cal_table[chan - cal_low +2][0]); break;
						case 1:	daydata->records[r].stokes.I /= ( cal_table[0][1] * daydata->records[r].DEC + cal_table[chan - cal_low +2][1]); break;
						case 2: daydata->records[r].stokes.I /= ( cal_table[0][2] * daydata->records[r].DEC + cal_table[chan - cal_low +2][2]); break;
						case 3: daydata->records[r].stokes.I /= ( cal_table[0][3] * daydata->records[r].DEC + cal_table[chan - cal_low +2][3]); break;
						case 4: daydata->records[r].stokes.I /= ( cal_table[0][4] * daydata->records[r].DEC + cal_table[chan - cal_low +2][4]); break;
						case 5: daydata->records[r].stokes.I /= ( cal_table[0][5] * daydata->records[r].DEC + cal_table[chan - cal_low +2][5]); break;
						case 6: daydata->records[r].stokes.I /= ( cal_table[0][6] * daydata->records[r].DEC + cal_table[chan - cal_low +2][6]); break; 
						default: break;
						}
					}
				}
			}				
		}
*/


if(chan == 0)
	{
	for(day=0; day<wappdata->numDays; day++) 
		{
		FluxDayData * daydata = &wappdata->daydata[day];
		for(r=0; r<daydata->numRecords; r++)
			{
			if(isfinite(daydata->records[r].stokes.I))
				{
				switch(day%7)
					{
					case 0:	daydata->records[r].stokes.I /= ( cal_table[0][0] * daydata->records[r].DEC + cal_table[1][0]); break;
					case 1:	daydata->records[r].stokes.I /= ( cal_table[0][1] * daydata->records[r].DEC + cal_table[1][1]); break;
					case 2: daydata->records[r].stokes.I /= ( cal_table[0][2] * daydata->records[r].DEC + cal_table[1][2]); break;
					case 3: daydata->records[r].stokes.I /= ( cal_table[0][3] * daydata->records[r].DEC + cal_table[1][3]); break;
					case 4: daydata->records[r].stokes.I /= ( cal_table[0][4] * daydata->records[r].DEC + cal_table[1][4]); break;
					case 5: daydata->records[r].stokes.I /= ( cal_table[0][5] * daydata->records[r].DEC + cal_table[1][5]); break;
					case 6: daydata->records[r].stokes.I /= ( cal_table[0][6] * daydata->records[r].DEC + cal_table[1][6]); break; 
					default: break;
					}
				}
			}
		}	
	}
else if(chan >= cal_low && chan <= cal_high)
		{
		for(day=0; day<wappdata->numDays; day++) 
			{
			FluxDayData * daydata = &wappdata->daydata[day];
			for(r=0; r<daydata->numRecords; r++)
				{
				if(isfinite(daydata->records[r].stokes.I))
					{
					switch(day%7)
						{
						case 0:	daydata->records[r].stokes.I /= ( cal_table[0][0] * daydata->records[r].DEC + cal_table[chan - cal_low +2][0]); break;
						case 1:	daydata->records[r].stokes.I /= ( cal_table[0][1] * daydata->records[r].DEC + cal_table[chan - cal_low +2][1]); break;
						case 2: daydata->records[r].stokes.I /= ( cal_table[0][2] * daydata->records[r].DEC + cal_table[chan - cal_low +2][2]); break;
						case 3: daydata->records[r].stokes.I /= ( cal_table[0][3] * daydata->records[r].DEC + cal_table[chan - cal_low +2][3]); break;
						case 4: daydata->records[r].stokes.I /= ( cal_table[0][4] * daydata->records[r].DEC + cal_table[chan - cal_low +2][4]); break;
						case 5: daydata->records[r].stokes.I /= ( cal_table[0][5] * daydata->records[r].DEC + cal_table[chan - cal_low +2][5]); break;
						case 6: daydata->records[r].stokes.I /= ( cal_table[0][6] * daydata->records[r].DEC + cal_table[chan - cal_low +2][6]); break; 
						default: break;
						}
					}
				}
			}				
		}

}
//-------------------------------------------------------------------------------
void calculate_dec_dependence(FluxWappData * wappdata, int order, int chan, float *cIc, float *cQc, float *cUc, float *cVc, int avg)
{
int d,beam;
for(d=0; d<wappdata->numDays; d++) 
	{
	//printf("Day %d of %d", d+1, wappdata->numDays);
        if(!strcmp(wappdata->wapp,"multibeam"))
	{
		beam = d%7;
	}
	else
		beam = atoi(&wappdata->wapp[4]);

	day_dec_dependence(wappdata, d, beam, order, chan, cIc, cQc, cUc, cVc, avg);
	}
}
//-------------------------------------------------------------------------------

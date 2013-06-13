#include "rfi.h"
#include <stdlib.h>
#include <math.h>

void mark_bad_channels(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[],
		int *badchannels) {
	int k, i;

	for (k = lowchan; k < highchan; k++) {
		if (badchannels[k]) {
			for (i = 0; i < size; i++) {
				dataset[i].flagRFI[k] |= RFI_CALOFF_XX;
				dataset[i].flagRFI[k] |= RFI_CALOFF_YY;
				dataset[i].flagRFI[k] |= RFI_CALON_XX;
				dataset[i].flagRFI[k] |= RFI_CALON_YY;
				dataset[i].flagRFI[k] |= RFI_CALOFF_XY;
				dataset[i].flagRFI[k] |= RFI_CALOFF_YX;
				dataset[i].flagRFI[k] |= RFI_CALON_XY;
				dataset[i].flagRFI[k] |= RFI_CALON_YX;
			}
		}
	}
}


void rfi_detection_frequency_domain(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband,float freq[])
{
	int i, j, k, N;
	float mean[8], sigma[8], diffs[size], diff, delta, minf = hidrogenfreq - hidrogenband, maxf = hidrogenfreq + hidrogenband;
	int outlierFound;
	char filename[50];
	FILE *diffoutput;
	int diffcounter = 0;

	for (i = 0; i < size; i++) {
		if (!dataset[i].flagBAD) {
			for (k = lowchan; k < highchan; k++)
				dataset[i].flagRFI[k] = RFI_NONE;

			float diffxxon[MAX_CHANNELS], diffyyon[MAX_CHANNELS],diffxyon[MAX_CHANNELS],diffyxon[MAX_CHANNELS];
			float diffxxoff[MAX_CHANNELS], diffyyoff[MAX_CHANNELS],diffxyoff[MAX_CHANNELS],diffyxoff[MAX_CHANNELS];

			for (k = lowchan; k < highchan - 1; k++) {
				if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE) {
					diffxxoff[k] = dataset[i].caloff.xx[k + 1] - dataset[i].caloff.xx[k];
					diffyyoff[k] = dataset[i].caloff.yy[k + 1] - dataset[i].caloff.yy[k];
					diffxxon[k] = dataset[i].calon.xx[k + 1] - dataset[i].calon.xx[k];
					diffyyon[k] = dataset[i].calon.yy[k + 1] - dataset[i].calon.yy[k];

					diffxyoff[k] = dataset[i].caloff.xy[k + 1] - dataset[i].caloff.xy[k];
					diffyxoff[k] = dataset[i].caloff.yx[k + 1] - dataset[i].caloff.yx[k];
					diffxyon[k] = dataset[i].calon.xy[k + 1] - dataset[i].calon.xy[k];
					diffyxon[k] = dataset[i].calon.yx[k + 1] - dataset[i].calon.yx[k];
				}
			}

			do {
				N = 0;
				for (j = 0; j < 8; j++) {
					mean[j] = 0.0;
					sigma[j] = 0.0;
				}

				for (k = lowchan; k < highchan - 1; k++) {
					if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE && (freq[k] < minf || freq[k] > maxf) )
					{
						mean[0] += diffxxoff[k];
						mean[1] += diffyyoff[k];
						mean[2] += diffxxon[k];
						mean[3] += diffyyon[k];
						mean[4] += diffxyoff[k];
						mean[5] += diffyxoff[k];
						mean[6] += diffxyon[k];
						mean[7] += diffyxon[k];

						N++;
					}
				}

				mean[0] /= N;
				mean[1] /= N;
				mean[2] /= N;
				mean[3] /= N;
				mean[4] /= N;
				mean[5] /= N;
				mean[6] /= N;
				mean[7] /= N;


				if (fabs(dataset[i].RA - 71.390589) < 0.001) {
					for (j = 0; j < 8; j++) {
						printf("mean %d = %f\n", j, mean[j]);
					}
				}

				for (k = lowchan; k < highchan - 1; k++) {
					if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE && (freq[k] < minf || freq[k] > maxf)) {
						sigma[0] += (diffxxoff[k] - mean[0]) * (diffxxoff[k] - mean[0]);
						sigma[1] += (diffyyoff[k] - mean[1]) * (diffyyoff[k] - mean[1]);
						sigma[2] += (diffxxon[k] - mean[2]) * (diffxxon[k] - mean[2]);
						sigma[3] += (diffyyon[k] - mean[3]) * (diffyyon[k] - mean[3]);

						sigma[4] += (diffxyoff[k] - mean[4]) * (diffxyoff[k] - mean[4]);
						sigma[5] += (diffyxoff[k] - mean[5]) * (diffyxoff[k] - mean[5]);
						sigma[6] += (diffxyon[k] - mean[6]) * (diffxyon[k] - mean[6]);
						sigma[7] += (diffyxon[k] - mean[7]) * (diffyxon[k] - mean[7]);

					}
				}

				for (j = 0; j < 8; j++) {
					if (N > 1)
						sigma[j] = sqrt(sigma[j] / (N - 1));
					else sigma[j] = sqrt(sigma[j]);

					if (fabs(dataset[i].RA - 71.390589) < 0.001) {
						printf("sigma %d = %f\n", j, sigma[j]);
					}
				}

				if (fabs(dataset[i].RA - 71.390589) < 0.001) {
					sprintf(filename, "diff_%f_%f.dat%d", dataset[i].RA, dataset[i].DEC, diffcounter++);
					diffoutput = fopen(filename, "w");

					if (diffoutput == NULL) {
						printf("diffoutput was null!!!!\n");
						exit(1);
					}

					for (k = lowchan; k < highchan - 1; k++) {
						fprintf(diffoutput, "%f %f %f %f %f %f %f %f ", diffxxoff[k], diffyyoff[k], diffxxon[k], diffyyon[k], diffxyoff[k], diffyxoff[k], diffxyon[k], diffyxon[k]);
						if (dataset[i].flagRFI[k] != RFI_NONE || dataset[i].flagRFI[k + 1] != RFI_NONE)
							fprintf(diffoutput, "%f\n", (sigma[0] * numSigma));
						else fprintf(diffoutput, "%f\n", 0.0);
					}
					fclose(diffoutput);
				}


				outlierFound = 0;
				for (k = lowchan; k < highchan - 1; k++) {
					if (freq[k] < minf || freq[k] > maxf) {
						//if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE) {

						if (fabs(diffxxoff[k] - mean[0]) > numSigma * sigma[0]) {
							if (dataset[i].flagRFI[k] == RFI_NONE || dataset[i].flagRFI[k + 1] == RFI_NONE) outlierFound = 1;
							dataset[i].flagRFI[k] |= RFI_CALOFF_XX;
							dataset[i].flagRFI[k + 1] |= RFI_CALOFF_XX;

						}
						if (fabs(diffyyoff[k] - mean[1]) > numSigma * sigma[1]) {
							if (dataset[i].flagRFI[k] == RFI_NONE || dataset[i].flagRFI[k + 1] == RFI_NONE) outlierFound = 1;
							dataset[i].flagRFI[k] |= RFI_CALOFF_YY;
							dataset[i].flagRFI[k + 1] |= RFI_CALOFF_YY;
						}
						if (fabs(diffxxon[k] - mean[2]) > numSigma * sigma[2]) {
							if (dataset[i].flagRFI[k] == RFI_NONE || dataset[i].flagRFI[k + 1] == RFI_NONE) outlierFound = 1;
							dataset[i].flagRFI[k] |= RFI_CALON_XX;
							dataset[i].flagRFI[k + 1] |= RFI_CALON_XX;
						}
						if (fabs(diffyyon[k] - mean[3]) > numSigma * sigma[3]) {
							if (dataset[i].flagRFI[k] == RFI_NONE || dataset[i].flagRFI[k + 1] == RFI_NONE) outlierFound = 1;
							dataset[i].flagRFI[k] |= RFI_CALON_YY;
							dataset[i].flagRFI[k + 1] |= RFI_CALON_YY;
						}

						// Cross correlations
						if (fabs(diffxyoff[k] - mean[4]) > numSigma * sigma[4]) {
							if (dataset[i].flagRFI[k] == RFI_NONE || dataset[i].flagRFI[k + 1] == RFI_NONE) outlierFound = 1;
							dataset[i].flagRFI[k] |= RFI_CALOFF_XX;
							dataset[i].flagRFI[k + 1] |= RFI_CALOFF_XX;

						}
						if (fabs(diffyxoff[k] - mean[5]) > numSigma * sigma[5]) {
							if (dataset[i].flagRFI[k] == RFI_NONE || dataset[i].flagRFI[k + 1] == RFI_NONE) outlierFound = 1;
							dataset[i].flagRFI[k] |= RFI_CALOFF_YY;
							dataset[i].flagRFI[k + 1] |= RFI_CALOFF_YY;
						}
						if (fabs(diffxyon[k] - mean[6]) > numSigma * sigma[6]) {
							if (dataset[i].flagRFI[k] == RFI_NONE || dataset[i].flagRFI[k + 1] == RFI_NONE) outlierFound = 1;
							dataset[i].flagRFI[k] |= RFI_CALON_XX;
							dataset[i].flagRFI[k + 1] |= RFI_CALON_XX;
						}
						if (fabs(diffyxon[k] - mean[7]) > numSigma * sigma[7]) {
							if (dataset[i].flagRFI[k] == RFI_NONE || dataset[i].flagRFI[k + 1] == RFI_NONE) outlierFound = 1;
							dataset[i].flagRFI[k] |= RFI_CALON_YY;
							dataset[i].flagRFI[k + 1] |= RFI_CALON_YY;
						}

						//}
					}
				}

				if (fabs(dataset[i].RA - 71.390589) < 0.001) {
					printf("outlierfound is %d at RA %f\n", outlierFound, dataset[i].RA);
				}
			}
			while (outlierFound);
		}
	}
}


// RFI detection in frequency domain 
// first order differences
// dataset - the channel data for a particular time step
// lowchan - the lowest channel number to do computations for
// highchan - the highest channel number to do computations for
// numSigma - the number of sigma away from the mean the difference should be rejected
// ignore hidrogenfeq(+-)hidrogenband
void rfi_detection_frequency_domain_gobble(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[]) {
	int i, j, k, N;
	float sum[4], mean[4], sigma[4], diffs[size], diff, delta, minf = hidrogenfreq - hidrogenband, maxf = hidrogenfreq + hidrogenband;
	char outlierFound;


	for (i = 0; i < size; i++) {
		if (!dataset[i].flagBAD) {
			for (k = lowchan; k < highchan; k++)
				dataset[i].flagRFI[k] = RFI_NONE;

			float diffxxon[MAX_CHANNELS], diffyyon[MAX_CHANNELS]; //,diffxyon[MAX_CHANNELS],diffyxon[MAX_CHANNELS];
			float diffxxoff[MAX_CHANNELS], diffyyoff[MAX_CHANNELS]; //,diffxyoff[MAX_CHANNELS],diffyxoff[MAX_CHANNELS];

			for (k = lowchan; k < highchan - 1; k++) {
				if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE) {
					diffxxoff[k] = dataset[i].caloff.xx[k + 1] - dataset[i].caloff.xx[k];
					diffyyoff[k] = dataset[i].caloff.yy[k + 1] - dataset[i].caloff.yy[k];

					diffxxon[k] = dataset[i].calon.xx[k + 1] - dataset[i].calon.xx[k];
					diffyyon[k] = dataset[i].calon.yy[k + 1] - dataset[i].calon.yy[k];
				}
			}

			char filename[50];
			FILE *diffoutput;


			int diffcounter =0;

			do {

				printf("RFI iteration %d\n", diffcounter);

				//if( i % 5000 == 0 ) {
				if( fabs(dataset[i].RA - 71.390589) < 0.001 )
				{
					//if( dataset[i].RA > 70.5 && dataset[i].RA < 71.5) {
				sprintf( filename, "diff_%f_%f.dat%d", dataset[i].RA, dataset[i].DEC, diffcounter );
				diffoutput = fopen( filename, "w" );

				if( diffoutput == NULL ) {
					printf("diffoutput was null!!!!\n");
					exit(1);
				}

				for (k = lowchan; k < highchan - 1; k++) {
					fprintf(diffoutput, "%f %f %f %f\n", diffxxoff[k], diffyyoff[k], diffxxon[k], diffyyon[k]);
				}

				fclose( diffoutput);
				}


				N = 0;
				for (j = 0; j < 4; j++) {
					mean[j] = 0;
					sigma[j] = 0;
				}

				for (k = lowchan; k < highchan - 1; k++) {
					if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE) {
						mean[0] += diffxxoff[k];
						mean[1] += diffyyoff[k];
						mean[2] += diffxxon[k];
						mean[3] += diffyyon[k];

						N++;
					}
				}

				mean[0] /= N;
				mean[1] /= N;
				mean[2] /= N;
				mean[3] /= N;

				for (k = lowchan; k < highchan - 1; k++) {
					if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE) {
						sigma[0] += (diffxxoff[k] - mean[0]) * (diffxxoff[k] - mean[0]);
						sigma[1] += (diffyyoff[k] - mean[1]) * (diffyyoff[k] - mean[1]);
						sigma[2] += (diffxxon[k] - mean[2]) * (diffxxon[k] - mean[2]);
						sigma[3] += (diffyyon[k] - mean[3]) * (diffyyon[k] - mean[3]);
					}
				}

				for (j = 0; j < 4; j++) {
					if (N > 1)
						sigma[j] = sqrt(sigma[j] / (N - 1));
					else sigma[j] = sqrt(sigma[j]);
				}

				outlierFound = 0;
				for (k = lowchan; k < highchan - 1; k++) {
					if (freq[k] < minf || freq[k] > maxf) {

						int RFIstart = k;
						int step = 1;


						while (dataset[i].flagRFI[RFIstart] != RFI_NONE && RFIstart < (highchan-1) ) {
							RFIstart++;
						}


						while (dataset[i].flagRFI[RFIstart+step] != RFI_NONE && RFIstart+step < highchan ) {
							step++;
						}

						///////  XX OFF
						if (fabs(dataset[i].caloff.xx[RFIstart + step] - dataset[i].caloff.xx[RFIstart] - mean[0]) > numSigma * sigma[0]) {
							dataset[i].flagRFI[RFIstart + step] |= RFI_CALOFF_XX;
							outlierFound = 1;
							if( RFIstart + step < highchan - 1 ) step++;
							else break;

							while (fabs(dataset[i].caloff.xx[RFIstart + step] - dataset[i].caloff.xx[RFIstart] - mean[0]) > numSigma * sigma[0]) {
								dataset[i].flagRFI[RFIstart + step] |= RFI_CALOFF_XX;
								if( RFIstart + step < highchan - 1 ) step++;
								else break;
							}
						}
						k = RFIstart;
					}
				}

				for (k = lowchan; k < highchan - 1; k++) {
					if (freq[k] < minf || freq[k] > maxf) {

						int RFIstart = k;
						int step = 1;

						while (dataset[i].flagRFI[RFIstart] == RFI_CALOFF_YY) {
							RFIstart++;
						}

						///////  YY OFF
						if (fabs(dataset[i].caloff.yy[RFIstart + step] - dataset[i].caloff.yy[RFIstart] - mean[1]) > numSigma * sigma[1]) {
							dataset[i].flagRFI[RFIstart + step] |= RFI_CALOFF_YY;
							outlierFound = 1;
							step++;

							while (fabs(dataset[i].caloff.yy[RFIstart + step] - dataset[i].caloff.yy[RFIstart] - mean[1]) > numSigma * sigma[1]) {
								dataset[i].flagRFI[RFIstart + step] |= RFI_CALOFF_YY;
								step++;
							}
						}
					}
				}

				for (k = lowchan; k < highchan - 1; k++) {
					if (freq[k] < minf || freq[k] > maxf) {

						int RFIstart = k;
						int step = 1;

						while (dataset[i].flagRFI[RFIstart] == RFI_CALON_XX) {
							RFIstart++;
						}

						///////  XX ON
						if (fabs(dataset[i].calon.xx[RFIstart + step] - dataset[i].calon.xx[RFIstart] - mean[2]) > numSigma * sigma[2]) {
							dataset[i].flagRFI[RFIstart + step] |= RFI_CALON_XX;
							outlierFound = 1;
							step++;

							while (fabs(dataset[i].calon.xx[RFIstart + step] - dataset[i].calon.xx[RFIstart] - mean[2]) > numSigma * sigma[2]) {
								dataset[i].flagRFI[RFIstart + step] |= RFI_CALON_XX;
								step++;
							}
						}
					}
				}

				for (k = lowchan; k < highchan - 1; k++) {
					if (freq[k] < minf || freq[k] > maxf) {

						int RFIstart = k;
						int step = 1;

						while (dataset[i].flagRFI[RFIstart] == RFI_CALON_YY) {
							RFIstart++;
						}

						///////  YY ON
						if (fabs(dataset[i].calon.yy[RFIstart + step] - dataset[i].calon.yy[RFIstart] - mean[3]) > numSigma * sigma[3]) {
							dataset[i].flagRFI[RFIstart + step] |= RFI_CALON_YY;
							outlierFound = 1;
							step++;

							while (fabs(dataset[i].calon.yy[RFIstart + step] - dataset[i].calon.yy[RFIstart] - mean[3]) > numSigma * sigma[3]) {
								dataset[i].flagRFI[RFIstart + step] |= RFI_CALON_YY;
								step++;
							}
						}
					}
				}

				diffcounter++;
			}
			while (outlierFound);
		}
	}
}

// RFI detection in frequency domain
// first order differences
// dataset - the channel data for a particular time step
// lowchan - the lowest channel number to do computations for
// highchan - the highest channel number to do computations for
// numSigma - the number of sigma away from the mean the difference should be rejected
// ignore hidrogenfeq(+-)hidrogenband
void rfi_detection_frequency_domain_old(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband,
		float freq[]) {
	int i, j, k, N;
	float mean[8], sigma[8], diff, delta, minf = hidrogenfreq - hidrogenband, maxf = hidrogenfreq + hidrogenband;
	char outlierFound;

	for (i = 0; i < size; i++) {
		if (!dataset[i].flagBAD) {
			for (k = lowchan; k < highchan; k++)
				dataset[i].flagRFI[k] = RFI_NONE;
			do {
				N = 0;
				for (j = 0; j < 8; j++) {
					mean[j] = 0;
					sigma[j] = 0;
				}
				for (k = lowchan; k < highchan - 1; k++) {
					if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE) {
						N++;
						diff = dataset[i].caloff.xx[k + 1] - dataset[i].caloff.xx[k];
						delta = diff - mean[0];
						mean[0] += delta / N;
						sigma[0] += delta * (diff - mean[0]);
						diff = dataset[i].caloff.yy[k + 1] - dataset[i].caloff.yy[k];
						delta = diff - mean[1];
						mean[1] += delta / N;
						sigma[1] += delta * (diff - mean[1]);
						diff = dataset[i].calon.xx[k + 1] - dataset[i].calon.xx[k];
						delta = diff - mean[2];
						mean[2] += delta / N;
						sigma[2] += delta * (diff - mean[2]);
						diff = dataset[i].calon.yy[k + 1] - dataset[i].calon.yy[k];
						delta = diff - mean[3];
						mean[3] += delta / N;
						sigma[3] += delta * (diff - mean[3]);
						diff = dataset[i].caloff.xy[k + 1] - dataset[i].caloff.xy[k];
						delta = diff - mean[4];
						mean[4] += delta / N;
						sigma[4] += delta * (diff - mean[4]);
						diff = dataset[i].caloff.yx[k + 1] - dataset[i].caloff.yx[k];
						delta = diff - mean[5];
						mean[5] += delta / N;
						sigma[5] += delta * (diff - mean[5]);
						diff = dataset[i].calon.xy[k + 1] - dataset[i].calon.xy[k];
						delta = diff - mean[6];
						mean[6] += delta / N;
						sigma[6] += delta * (diff - mean[6]);
						diff = dataset[i].calon.yx[k + 1] - dataset[i].calon.yx[k];
						delta = diff - mean[7];
						mean[7] += delta / N;
						sigma[7] += delta * (diff - mean[7]);
					}
				}
				for (j = 0; j < 8; j++) {
					if (N > 1)
						sigma[j] = sqrt(sigma[j] / (N - 1));
					else sigma[j] = sqrt(sigma[j]);
				}
				outlierFound = 0;
				for (k = lowchan; k < highchan - 1; k++) {
					if (freq[k] < minf || freq[k] > maxf) {
						if (dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k + 1] == RFI_NONE) {
							if (fabs(dataset[i].caloff.xx[k + 1] - dataset[i].caloff.xx[k] - mean[0]) > numSigma * sigma[0]) {
								dataset[i].flagRFI[k] |= RFI_CALOFF_XX;
								dataset[i].flagRFI[k + 1] |= RFI_CALOFF_XX;
								outlierFound = 1;
							}
							if (fabs(dataset[i].caloff.yy[k + 1] - dataset[i].caloff.yy[k] - mean[1]) > numSigma * sigma[1]) {
								dataset[i].flagRFI[k] |= RFI_CALOFF_YY;
								dataset[i].flagRFI[k + 1] |= RFI_CALOFF_YY;
								outlierFound = 1;
							}
							if (fabs(dataset[i].calon.xx[k + 1] - dataset[i].calon.xx[k] - mean[2]) > numSigma * sigma[2]) {
								dataset[i].flagRFI[k] |= RFI_CALON_XX;
								dataset[i].flagRFI[k + 1] |= RFI_CALON_XX;
								outlierFound = 1;
							}
							if (fabs(dataset[i].calon.yy[k + 1] - dataset[i].calon.yy[k] - mean[3]) > numSigma * sigma[3]) {
								dataset[i].flagRFI[k] |= RFI_CALON_YY;
								dataset[i].flagRFI[k + 1] |= RFI_CALON_YY;
								outlierFound = 1;
							}
							if (fabs(dataset[i].caloff.xy[k + 1] - dataset[i].caloff.xy[k] - mean[4]) > numSigma * sigma[4]) {
								dataset[i].flagRFI[k] |= RFI_CALOFF_XY;
								dataset[i].flagRFI[k + 1] |= RFI_CALOFF_XY;
								outlierFound = 1;
							}
							if (fabs(dataset[i].caloff.yx[k + 1] - dataset[i].caloff.yx[k] - mean[5]) > numSigma * sigma[5]) {
								dataset[i].flagRFI[k] |= RFI_CALOFF_YX;
								dataset[i].flagRFI[k + 1] |= RFI_CALOFF_YX;
								outlierFound = 1;
							}
							if (fabs(dataset[i].calon.xy[k + 1] - dataset[i].calon.xy[k] - mean[6]) > numSigma * sigma[6]) {
								dataset[i].flagRFI[k] |= RFI_CALON_XY;
								dataset[i].flagRFI[k + 1] |= RFI_CALON_XY;
								outlierFound = 1;
							}
							if (fabs(dataset[i].calon.yx[k + 1] - dataset[i].calon.yx[k] - mean[7]) > numSigma * sigma[7]) {
								dataset[i].flagRFI[k] |= RFI_CALON_YX;
								dataset[i].flagRFI[k + 1] |= RFI_CALON_YX;
								outlierFound = 1;
							}
						}
					}
				}
			}
			while (outlierFound);
		}
	}
}

void rfi_detection_time_domain1(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq,
		float hidrogenband, float freq[]) {
// sigma
	int i, j, k, N;
	float mean[8], sigma[8], diff, delta, minf = hidrogenfreq - hidrogenband, maxf = hidrogenfreq + hidrogenband, nSigma = numSigma;
	char outlierFound;
	float signal[size][8];
	int rfi_count = 0;

	for (i = 0; i < size; i++) {
		if (!dataset[i].flagBAD) {
			N = 0;
			for (j = 0; j < 8; j++) {
				mean[j] = 0;
				sigma[j] = 0;
			}
			for (k = lowchan; k < highchan - 1; k++) {
				N++;
				diff = dataset[i].caloff.xx[k + 1] - dataset[i].caloff.xx[k];
				delta = diff - mean[0];
				mean[0] += delta / N;
				sigma[0] += delta * (diff - mean[0]);
				diff = dataset[i].caloff.yy[k + 1] - dataset[i].caloff.yy[k];
				delta = diff - mean[1];
				mean[1] += delta / N;
				sigma[1] += delta * (diff - mean[1]);
				diff = dataset[i].calon.xx[k + 1] - dataset[i].calon.xx[k];
				delta = diff - mean[2];
				mean[2] += delta / N;
				sigma[2] += delta * (diff - mean[2]);
				diff = dataset[i].calon.yy[k + 1] - dataset[i].calon.yy[k];
				delta = diff - mean[3];
				mean[3] += delta / N;
				sigma[3] += delta * (diff - mean[3]);
				diff = dataset[i].caloff.xy[k + 1] - dataset[i].caloff.xy[k];
				delta = diff - mean[4];
				mean[4] += delta / N;
				sigma[4] += delta * (diff - mean[4]);
				diff = dataset[i].caloff.yx[k + 1] - dataset[i].caloff.yx[k];
				delta = diff - mean[5];
				mean[5] += delta / N;
				sigma[5] += delta * (diff - mean[5]);
				diff = dataset[i].calon.xy[k + 1] - dataset[i].calon.xy[k];
				delta = diff - mean[6];
				mean[6] += delta / N;
				sigma[6] += delta * (diff - mean[6]);
				diff = dataset[i].calon.yx[k + 1] - dataset[i].calon.yx[k];
				delta = diff - mean[7];
				mean[7] += delta / N;
				sigma[7] += delta * (diff - mean[7]);
			}
			for (j = 0; j < 8; j++) {
				if (N > 1)
					signal[i][j] = sqrt(sigma[j] / (N - 1));
				else signal[i][j] = sqrt(sigma[j]);
			}
		}
	}
	do {
		N = 0;
		for (j = 0; j < 8; j++) {
			mean[j] = 0;
			sigma[j] = 0;
		}
		for (i = 1; i < size - 1; i++) {
			if (dataset[i].flagBAD && dataset[i - 1].flagBAD && dataset[i + 1].flagBAD) continue;
			N++;
			for (j = 0; j < 8; j++) {
				diff = signal[i + 1][j] - 2 * signal[i][j] + signal[i - 1][j];
				delta = diff - mean[j];
				mean[j] += delta / N;
				sigma[j] += delta * (diff - mean[j]);
			}
		}
		for (j = 0; j < 8; j++) {
			if (N > 1)
				sigma[j] = sqrt(sigma[j] / (N - 1));
			else sigma[j] = sqrt(sigma[j]);
		}
		outlierFound = 0;
		for (i = 1; i < size - 1; i++) {
			if (field[0] == 'N' && field[1] == '1' && ((dataset[i].RA > 82.9333 && dataset[i].RA < 84.0666) || (dataset[i].RA > 68.95 && dataset[i].RA < 69.6)))
				continue; // Just for N1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (dataset[i].flagBAD && dataset[i - 1].flagBAD && dataset[i + 1].flagBAD) continue;
			if (fabs(signal[i + 1][0] - 2 * signal[i][0] + signal[i - 1][0] - mean[0]) > numSigma * sigma[0]
					|| fabs(signal[i + 1][1] - 2 * signal[i][1] + signal[i - 1][1] - mean[1]) > numSigma * sigma[1]
					|| fabs(signal[i + 1][2] - 2 * signal[i][2] + signal[i - 1][2] - mean[2]) > numSigma * sigma[2]
					|| fabs(signal[i + 1][3] - 2 * signal[i][3] + signal[i - 1][3] - mean[3]) > numSigma * sigma[3]
					|| fabs(signal[i + 1][4] - 2 * signal[i][4] + signal[i - 1][4] - mean[4]) > numSigma * sigma[4]
					|| fabs(signal[i + 1][5] - 2 * signal[i][5] + signal[i - 1][5] - mean[5]) > numSigma * sigma[5]
					|| fabs(signal[i + 1][6] - 2 * signal[i][6] + signal[i - 1][6] - mean[6]) > numSigma * sigma[6]
					|| fabs(signal[i + 1][7] - 2 * signal[i][7] + signal[i - 1][7] - mean[7]) > numSigma * sigma[7]) {
				dataset[i - 1].flagBAD = 1;
				dataset[i].flagBAD = 1;
				dataset[i + 1].flagBAD = 1;
				for (k = lowchan; k < highchan; k++) {
					{
						dataset[i - 1].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i + 1].flagRFI[k] |= RFI_OUTOFBAND;
						rfi_count++;
					}
				}
				outlierFound = 1;
			}
		}
	}
	while (outlierFound);

	printf("!!! Found %d RFI\n", rfi_count);
}

void rfi_detection_time_domain2(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq,
		float hidrogenband, float freq[]) {
// average
	int i, k, N, j, cur = 0;
	float mean[8], sigma[8], diff, delta, minf = hidrogenfreq - hidrogenband, maxf = hidrogenfreq + hidrogenband, nSigma = numSigma;
	char outlierFound;
	float signal[size][8];
	float diffs[size][8];
	int rfi_count = 0;
	int rfi_chancount = 0;
	int outliercount = 0;

	float RA, DEC, radius = 0.06666;
	float s[10000][3];
	int sourcecount;
	int exceptionCount = 0;
	int badcount = 0;

	FILE *StrongSourceFile = fopen("../../rfi.list", "r");
	FILE *RFIin = fopen("rfitimein.ann", "w");
	FILE *AnnotationFile = fopen("rfiexception.ann", "w");
	fprintf(AnnotationFile, "COLOUR YELLOW\n");

	fprintf(RFIin, "COLOUR RED\n");

	if (StrongSourceFile != NULL) {
		sourcecount = jsd_line_count(StrongSourceFile);
		rewind(StrongSourceFile);

		for (i = 0; i < sourcecount; i++) {
			int num, ret;
			float fileradius;
			char line[80];

			if (fgets(line, 80, StrongSourceFile) == NULL) {
				printf("in rfi.c, fgets read error");
			}

			num = sscanf(line, "%f %f %f", &RA, &DEC, &fileradius);

			s[cur][0] = RA;
			s[cur][1] = DEC;
			s[cur][2] = radius;

			if (num == 3) {
				s[cur][2] = fileradius;
			}
			cur++;
		}
	}
	else {
		printf("No time domain RFI exception list file found\n");
	}

	int M = 0;
	int count[8] = {0,0,0,0,0,0,0,0};

	for (i = 0; i < size; i++) {
		if (!dataset[i].flagBAD) {
			for (j = 0; j < 8; j++)
				signal[i][j] = 0;
			N = 0;
			M++;
			for (k = lowchan; k < highchan; k++) {
				if( freq[k] < minf || freq[k] > maxf )
				{
					if( dataset[i].flagRFI[k] != RFI_CALOFF_XX ){
						signal[i][0] += dataset[i].caloff.xx[k];
						count[0]++;
					}
					if( dataset[i].flagRFI[k] != RFI_CALOFF_YY ) {
						signal[i][1] += dataset[i].caloff.yy[k];
						count[1]++;
					}
					if( dataset[i].flagRFI[k] != RFI_CALON_XX ){
						signal[i][2] += dataset[i].calon.xx[k];
						count[2]++;
					}
					if( dataset[i].flagRFI[k] != RFI_CALON_YY ){
						signal[i][3] += dataset[i].calon.yy[k];
						count[3]++;
					}
					if( dataset[i].flagRFI[k] != RFI_CALOFF_XY ){
						signal[i][4] += dataset[i].caloff.xy[k];
						count[4]++;
					}
					if( dataset[i].flagRFI[k] != RFI_CALOFF_YX ){
						signal[i][5] += dataset[i].caloff.yx[k];
						count[5]++;
					}
					if( dataset[i].flagRFI[k] != RFI_CALON_XY ){
						signal[i][6] += dataset[i].calon.xy[k];
						count[6]++;
					}
					if( dataset[i].flagRFI[k] != RFI_CALON_YX ){
						signal[i][7] += dataset[i].calon.yx[k];
						count[7]++;
					}
					//N++;
				}
			}

			for (j = 0; j < 8; j++) {
				signal[i][j] = signal[i][j] / count[j];
			}

		}
	}

	//printf("Counts %d %d %d %d %d %d %d %d", count[0], count[1],count[2],count[3],count[4],count[5],count[6],count[7]);

	do {
		N = 0;
		M = 0;
		for (j = 0; j < 8; j++) {
			mean[j] = 0;
			sigma[j] = 0;
		}

		for (i = 1; i < size - 1; i++) {
			if (dataset[i].flagBAD && dataset[i - 1].flagBAD && dataset[i + 1].flagBAD) continue;

			for (j = 0; j < 8; j++) {
				diffs[i][j] = signal[i + 1][j] - 2 * signal[i][j] + signal[i - 1][j];
				mean[j] += diffs[i][j];
			}
			M++;
		}

		for (j = 0; j < 8; j++) {
			mean[j] /= M;
			//printf("Mean %d = %.10f\n", j, mean[j] );
		}

		for (i = 1; i < size - 1; i++) {
			if (dataset[i].flagBAD && dataset[i - 1].flagBAD && dataset[i + 1].flagBAD) continue;

			for (j = 0; j < 8; j++) {
				sigma[j] += (diffs[i][j] - mean[j]) * (diffs[i][j] - mean[j]);
			}
			N++;
		}

		for (j = 0; j < 8; j++) {
			if (N > 1)
				sigma[j] = sqrt(sigma[j] / (N - 1));
			else sigma[j] = sqrt(sigma[j]);

			//printf("Sigma %d = %.10f\n", j, sigma[j] );
		}

		outlierFound = 0;
		for (i = 1; i < size - 1; i++) {
			//if(field[0] == 'N' && field[1] == '1' && ((dataset[i].RA>82.9333 && dataset[i].RA<84.0666) || (dataset[i].RA>68.95 && dataset[i].RA<69.6))) continue; // Just for N1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (dataset[i].flagBAD && dataset[i - 1].flagBAD && dataset[i + 1].flagBAD) continue;
			if (fabs(signal[i + 1][0] - 2 * signal[i][0] + signal[i - 1][0] - mean[0]) > numSigma * sigma[0]
					|| fabs(signal[i + 1][1] - 2 * signal[i][1] + signal[i - 1][1] - mean[1]) > numSigma * sigma[1]
					|| fabs(signal[i + 1][2] - 2 * signal[i][2] + signal[i - 1][2] - mean[2]) > numSigma * sigma[2]
					|| fabs(signal[i + 1][3] - 2 * signal[i][3] + signal[i - 1][3] - mean[3]) > numSigma * sigma[3]
					|| fabs(signal[i + 1][4] - 2 * signal[i][4] + signal[i - 1][4] - mean[4]) > numSigma * sigma[4]
					|| fabs(signal[i + 1][5] - 2 * signal[i][5] + signal[i - 1][5] - mean[5]) > numSigma * sigma[5]
					|| fabs(signal[i + 1][6] - 2 * signal[i][6] + signal[i - 1][6] - mean[6]) > numSigma * sigma[6]
					|| fabs(signal[i + 1][7] - 2 * signal[i][7] + signal[i - 1][7] - mean[7]) > numSigma * sigma[7]) {

				fprintf(RFIin, "CROSS %f %f 0.2 0.2\n", dataset[i].RA, dataset[i].DEC);

				int exception = 0;

				if (StrongSourceFile != NULL) {
					// found an outlier, check if its a strong source
					for (k = 0; k < sourcecount; k++) {
						// radius must be greater than 0
						if (s[k][2] > 0.00) {
							//if( sqrt( powf(dataset[i].RA - s[k][0],2) + powf(dataset[i].DEC - s[k][1],2) ) < s[k][2]  )
							if ((fabs(dataset[i].RA - s[k][0]) < s[k][2]) && (fabs(dataset[i].DEC - s[k][1]) < s[k][2])) {
								// GOT A HIT, remove the RFI flag
								//printf("Exception hit at %f %f in dataset at %f %f\n",  s[k][0], s[k][1], dataset[i].RA, dataset[i].DEC  );
								fprintf(AnnotationFile, "CROSS %f %f 0.2 0.2\n", dataset[i].RA, dataset[i].DEC);
								exception = 1;
								exceptionCount++;
								break;
							}
						}
					}
				}

				if (!exception) {
					dataset[i - 1].flagBAD = 1;
					dataset[i].flagBAD = 1;
					dataset[i + 1].flagBAD = 1;
					for (k = lowchan; k < highchan; k++) {
						rfi_chancount++;
						{
							dataset[i - 1].flagRFI[k] |= RFI_OUTOFBAND;
							dataset[i].flagRFI[k] |= RFI_OUTOFBAND;
							dataset[i + 1].flagRFI[k] |= RFI_OUTOFBAND;
						}
					}
					outlierFound = 1;
					outliercount++;
				}
			}
		}
	}
	while (outlierFound);

	fclose(AnnotationFile);
	fclose(RFIin);

	printf("Time domain. Found %d outliers and %d exceptions\n", outliercount, exceptionCount);
	//printf("!!! Found %d RFI points in %d channels and %d outliers\n", rfi_count, rfi_chancount, outliercount);
}

void rfi_detection_time_domain3(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq,
		float hidrogenband, float freq[]) {
// difference
	int i, j, k, N;
	float mean[10], sigma[10], diff, delta, minf = hidrogenfreq - hidrogenband, maxf = hidrogenfreq + hidrogenband, nSigma = numSigma;
	char outlierFound;
	float signal[size][10];

	for (i = 0; i < size; i++) {
		if (!dataset[i].flagBAD) {
			for (j = 0; j < 8; j++)
				signal[i][j] = 0;
			N = 0;
			for (k = lowchan; k < highchan; k++) {
				N++;
				signal[i][0] += dataset[i].caloff.xx[k];
				signal[i][1] += dataset[i].caloff.yy[k];
				signal[i][2] += dataset[i].calon.xx[k];
				signal[i][3] += dataset[i].calon.yy[k];
				signal[i][4] += dataset[i].caloff.xy[k];
				signal[i][5] += dataset[i].caloff.yx[k];
				signal[i][6] += dataset[i].calon.xy[k];
				signal[i][7] += dataset[i].calon.yx[k];
			}
			for (j = 0; j < 8; j++)
				signal[i][j] = signal[i][j] / N;
			signal[i][8] = signal[i][0] - signal[i][1];
			signal[i][9] = signal[i][2] - signal[i][3];
		}
	}

	do {
		N = 0;
		for (j = 0; j < 10; j++) {
			mean[j] = 0;
			sigma[j] = 0;
		}
		for (i = 1; i < size - 1; i++) {
			if (dataset[i].flagBAD && dataset[i - 1].flagBAD && dataset[i + 1].flagBAD) continue;
			N++;
			for (j = 0; j < 10; j++) {
				diff = signal[i + 1][j] - 2 * signal[i][j] + signal[i - 1][j];
				delta = diff - mean[j];
				mean[j] += delta / N;
				sigma[j] += delta * (diff - mean[j]);
			}
		}
		for (j = 0; j < 10; j++) {
			if (N > 1)
				sigma[j] = sqrt(sigma[j] / (N - 1));
			else sigma[j] = sqrt(sigma[j]);
		}
		outlierFound = 0;
		for (i = 1; i < size - 1; i++) {
			if (field[0] == 'N' && field[1] == '1' && ((dataset[i].RA > 82.9333 && dataset[i].RA < 84.0666) || (dataset[i].RA > 68.95 && dataset[i].RA < 69.6)))
				continue; // Just for N1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (dataset[i].flagBAD && dataset[i - 1].flagBAD && dataset[i + 1].flagBAD) continue;
			if (fabs(signal[i + 1][8] - 2 * signal[i][8] + signal[i - 1][8] - mean[8]) > numSigma * sigma[8]
					|| fabs(signal[i + 1][9] - 2 * signal[i][9] + signal[i - 1][9] - mean[9]) > numSigma * sigma[9]
					|| fabs(signal[i + 1][4] - 2 * signal[i][4] + signal[i - 1][4] - mean[4]) > numSigma * sigma[4]
					|| fabs(signal[i + 1][5] - 2 * signal[i][5] + signal[i - 1][5] - mean[5]) > numSigma * sigma[5]
					|| fabs(signal[i + 1][6] - 2 * signal[i][6] + signal[i - 1][6] - mean[6]) > numSigma * sigma[6]
					|| fabs(signal[i + 1][7] - 2 * signal[i][7] + signal[i - 1][7] - mean[7]) > numSigma * sigma[7]) {
				dataset[i - 1].flagBAD = 1;
				dataset[i].flagBAD = 1;
				dataset[i + 1].flagBAD = 1;
				for (k = lowchan; k < highchan; k++) {
					{
						dataset[i - 1].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i + 1].flagRFI[k] |= RFI_OUTOFBAND;
					}
				}
				outlierFound = 1;
			}
		}
	}
	while (outlierFound);
}

void outofbandrfi_ann(SpecRecord dataset[], int size, int lowchan, int highchan) {
	int i;

	FILE * file;
	file = fopen("outofbandrfi.ann", "w");
	fprintf(file, "COLOUR RED\n");

	/*
	 for (i = 0; i < size; i++)
	 {
	 if(dataset[i].flagRFI[lowchan] && RFI_OUTOFBAND) fprintf(file, "LINE W %i %f %i %f\n", i, dataset[i].RA, 30, dataset[i].AST);
	 }
	 fclose(file);*/

	//file = fopen("rfi_rd.ann", "w");
	for (i = 0; i < size; i++) {
		if (dataset[i].flagRFI[lowchan] == RFI_OUTOFBAND) fprintf(file, "CROSS %f %f 0.2 0.2\n", dataset[i].RA, dataset[i].DEC);
	}
	fclose(file);

}

void rfi_ann(SpecRecord dataset[], int size, int lowchan, int highchan, float freq[]) {
	int i, k;

	FILE * file;
	file = fopen("rfi.ann", "w");

	for (i = 0; i < size; i++) {
		for (k = lowchan; k < highchan; k++) {
			if (dataset[i].flagRFI[k] != RFI_NONE) fprintf(file, "%f %f\n", freq[k], dataset[i].AST);
		}
	}

	fclose(file);
}

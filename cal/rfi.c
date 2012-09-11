#include "rfi.h"
#include <stdlib.h>
#include <math.h>

void mark_bad_channels(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[], int *badchannels)
{
int k, i;

for(k = lowchan; k < highchan; k++)
	{
	if(badchannels[k])
		{
		for(i = 0; i < size; i++)
			{
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

// RFI detection in frequency domain 
// first order differences
// dataset - the channel data for a particular time step
// lowchan - the lowest channel number to do computations for
// highchan - the highest channel number to do computations for
// numSigma - the number of sigma away from the mean the difference should be rejected
// ignore hidrogenfeq(+-)hidrogenband
void rfi_detection_frequency_domain(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[])
{
	int i, j, k, N;
	float mean[8], sigma[8], diff, delta, minf=hidrogenfreq-hidrogenband, maxf=hidrogenfreq+hidrogenband;
	char outlierFound;

	for(i = 0; i < size; i++)
		{
		if(!dataset[i].flagBAD)
			{
			for(k = lowchan; k < highchan; k++) dataset[i].flagRFI[k] = RFI_NONE;
			do
				{
				N = 0; 
				for(j=0; j<8; j++){mean[j] = 0; sigma[j] = 0;} 
				for(k = lowchan; k < highchan - 1; k++)
					{
					if(dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k+1] == RFI_NONE)
						{
						N++;
						diff = dataset[i].caloff.xx[k+1] - dataset[i].caloff.xx[k];
						delta = diff - mean[0];
						mean[0] += delta/N;
						sigma[0] += delta*(diff - mean[0]);
						diff = dataset[i].caloff.yy[k+1] - dataset[i].caloff.yy[k];
						delta = diff - mean[1];
						mean[1] += delta/N;
						sigma[1] += delta*(diff - mean[1]);
						diff = dataset[i].calon.xx[k+1] - dataset[i].calon.xx[k];
						delta = diff - mean[2];
						mean[2] += delta/N;
						sigma[2] += delta*(diff - mean[2]);
						diff = dataset[i].calon.yy[k+1] - dataset[i].calon.yy[k];
						delta = diff - mean[3];
						mean[3] += delta/N;
						sigma[3] += delta*(diff - mean[3]);
						diff = dataset[i].caloff.xy[k+1] - dataset[i].caloff.xy[k];
						delta = diff - mean[4];
						mean[4] += delta/N;
						sigma[4] += delta*(diff - mean[4]);
						diff = dataset[i].caloff.yx[k+1] - dataset[i].caloff.yx[k];
						delta = diff - mean[5];
						mean[5] += delta/N;
						sigma[5] += delta*(diff - mean[5]);
						diff = dataset[i].calon.xy[k+1] - dataset[i].calon.xy[k];
						delta = diff - mean[6];
						mean[6] += delta/N;
						sigma[6] += delta*(diff - mean[6]);
						diff = dataset[i].calon.yx[k+1] - dataset[i].calon.yx[k];
						delta = diff - mean[7];
						mean[7] += delta/N;
						sigma[7] += delta*(diff - mean[7]);
						}
					}
				for(j=0; j<8; j++) {if(N > 1) sigma[j] = sqrt(sigma[j]/(N-1)); else sigma[j] = sqrt(sigma[j]);}
				outlierFound = 0;
				for(k = lowchan; k < highchan - 1; k++)
					{
					if(freq[k] < minf || freq[k] > maxf)
						{
						if(dataset[i].flagRFI[k] == RFI_NONE && dataset[i].flagRFI[k+1] == RFI_NONE)
							{
							if(fabs(dataset[i].caloff.xx[k+1] - dataset[i].caloff.xx[k] - mean[0]) > numSigma*sigma[0])
								{
								dataset[i].flagRFI[k] |= RFI_CALOFF_XX;
								dataset[i].flagRFI[k+1] |= RFI_CALOFF_XX;
								outlierFound = 1;
								}
							if(fabs(dataset[i].caloff.yy[k+1] - dataset[i].caloff.yy[k] - mean[1]) > numSigma*sigma[1])
								{
								dataset[i].flagRFI[k] |= RFI_CALOFF_YY;
								dataset[i].flagRFI[k+1] |= RFI_CALOFF_YY;
								outlierFound = 1;								
								}
							if(fabs(dataset[i].calon.xx[k+1] - dataset[i].calon.xx[k] - mean[2]) > numSigma*sigma[2])
								{
								dataset[i].flagRFI[k] |= RFI_CALON_XX;
								dataset[i].flagRFI[k+1] |= RFI_CALON_XX;
								outlierFound = 1;
								}
							if(fabs(dataset[i].calon.yy[k+1] - dataset[i].calon.yy[k] - mean[3]) > numSigma*sigma[3])
								{
								dataset[i].flagRFI[k] |= RFI_CALON_YY;
								dataset[i].flagRFI[k+1] |= RFI_CALON_YY;
								outlierFound = 1;
								}
							if(fabs(dataset[i].caloff.xy[k+1] - dataset[i].caloff.xy[k] - mean[4]) > numSigma*sigma[4])
								{
								dataset[i].flagRFI[k] |= RFI_CALOFF_XY;
								dataset[i].flagRFI[k+1] |= RFI_CALOFF_XY;
								outlierFound = 1;									
								}
							if(fabs(dataset[i].caloff.yx[k+1] - dataset[i].caloff.yx[k] - mean[5]) > numSigma*sigma[5])
								{
								dataset[i].flagRFI[k] |= RFI_CALOFF_YX;
								dataset[i].flagRFI[k+1] |= RFI_CALOFF_YX;
								outlierFound = 1;									
								}
							if(fabs(dataset[i].calon.xy[k+1] - dataset[i].calon.xy[k] - mean[6]) > numSigma*sigma[6])
								{
								dataset[i].flagRFI[k] |= RFI_CALON_XY;
								dataset[i].flagRFI[k+1] |= RFI_CALON_XY;
								outlierFound = 1;									
								}
							if(fabs(dataset[i].calon.yx[k+1] - dataset[i].calon.yx[k] - mean[7]) > numSigma*sigma[7]) 
								{
								dataset[i].flagRFI[k] |= RFI_CALON_YX;
								dataset[i].flagRFI[k+1] |= RFI_CALON_YX;
								outlierFound = 1;
								}
							}
						}
					}
				}
			while(outlierFound);
			}
		}
}

void rfi_detection_time_domain1(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[])
{
// sigma
	int i, j, k, N;
	float mean[8], sigma[8], diff, delta, minf=hidrogenfreq-hidrogenband, maxf=hidrogenfreq+hidrogenband, nSigma = numSigma;
	char outlierFound;
	float signal[size][8];
	
	for(i = 0; i < size; i++)
		{
		if(!dataset[i].flagBAD)
			{
			N = 0; 
			for(j=0; j<8; j++){mean[j] = 0; sigma[j] = 0;}
			for(k = lowchan; k < highchan - 1; k++)
				{
				N++;
				diff = dataset[i].caloff.xx[k+1] - dataset[i].caloff.xx[k];
				delta = diff - mean[0];
				mean[0] += delta/N;
				sigma[0] += delta*(diff - mean[0]);
				diff = dataset[i].caloff.yy[k+1] - dataset[i].caloff.yy[k];
				delta = diff - mean[1];
				mean[1] += delta/N;
				sigma[1] += delta*(diff - mean[1]);
				diff = dataset[i].calon.xx[k+1] - dataset[i].calon.xx[k];
				delta = diff - mean[2];
				mean[2] += delta/N;
				sigma[2] += delta*(diff - mean[2]);
				diff = dataset[i].calon.yy[k+1] - dataset[i].calon.yy[k];
				delta = diff - mean[3];
				mean[3] += delta/N;
				sigma[3] += delta*(diff - mean[3]);
				diff = dataset[i].caloff.xy[k+1] - dataset[i].caloff.xy[k];
				delta = diff - mean[4];
				mean[4] += delta/N;
				sigma[4] += delta*(diff - mean[4]);
				diff = dataset[i].caloff.yx[k+1] - dataset[i].caloff.yx[k];
				delta = diff - mean[5];
				mean[5] += delta/N;
				sigma[5] += delta*(diff - mean[5]);
				diff = dataset[i].calon.xy[k+1] - dataset[i].calon.xy[k];
				delta = diff - mean[6];
				mean[6] += delta/N;
				sigma[6] += delta*(diff - mean[6]);
				diff = dataset[i].calon.yx[k+1] - dataset[i].calon.yx[k];
				delta = diff - mean[7];
				mean[7] += delta/N;
				sigma[7] += delta*(diff - mean[7]);
				}
			for(j=0; j<8; j++){if(N>1) signal[i][j] = sqrt(sigma[j]/(N-1)); else signal[i][j] = sqrt(sigma[j]);}
			}
		}
	do
		{
		N = 0; 
		for(j=0; j<8; j++){mean[j] = 0; sigma[j] = 0;} 
		for(i = 1; i < size - 1; i++)
			{
			if(dataset[i].flagBAD && dataset[i-1].flagBAD && dataset[i+1].flagBAD) continue;
			N++;
			for(j=0; j<8; j++)
				{
				diff = signal[i+1][j] - 2*signal[i][j] + signal[i-1][j];
				delta = diff - mean[j];
				mean[j] += delta/N;
				sigma[j] += delta*(diff - mean[j]);
				}
			}
		for(j=0; j<8; j++){if(N>1) sigma[j] = sqrt(sigma[j]/(N-1)); else sigma[j] = sqrt(sigma[j]);}
		outlierFound = 0;
		for(i = 1; i < size - 1; i++)
			{
			if(field[0] == 'N' && field[1] == '1' && ((dataset[i].RA>82.9333 && dataset[i].RA<84.0666) || (dataset[i].RA>68.95 && dataset[i].RA<69.6))) continue; // Just for N1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(dataset[i].flagBAD && dataset[i-1].flagBAD && dataset[i+1].flagBAD) continue;
			if( fabs(signal[i+1][0] - 2*signal[i][0] + signal[i-1][0] - mean[0]) > numSigma*sigma[0] ||
				fabs(signal[i+1][1] - 2*signal[i][1] + signal[i-1][1] - mean[1]) > numSigma*sigma[1] ||
				fabs(signal[i+1][2] - 2*signal[i][2] + signal[i-1][2] - mean[2]) > numSigma*sigma[2] ||
				fabs(signal[i+1][3] - 2*signal[i][3] + signal[i-1][3] - mean[3]) > numSigma*sigma[3] ||
				fabs(signal[i+1][4] - 2*signal[i][4] + signal[i-1][4] - mean[4]) > numSigma*sigma[4] ||
				fabs(signal[i+1][5] - 2*signal[i][5] + signal[i-1][5] - mean[5]) > numSigma*sigma[5] ||
				fabs(signal[i+1][6] - 2*signal[i][6] + signal[i-1][6] - mean[6]) > numSigma*sigma[6] ||
				fabs(signal[i+1][7] - 2*signal[i][7] + signal[i-1][7] - mean[7]) > numSigma*sigma[7])
				{
				dataset[i-1].flagBAD = 1;
				dataset[i].flagBAD = 1;
				dataset[i+1].flagBAD = 1;
				for(k = lowchan; k < highchan; k++)
					{
						{
						dataset[i-1].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i+1].flagRFI[k] |= RFI_OUTOFBAND;
						}
					}
				outlierFound = 1;
				}
			}
		}
	while(outlierFound);
}

void rfi_detection_time_domain2(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[])
{
// average
	int i, j, k, N;
	float mean[8], sigma[8], diff, delta, minf=hidrogenfreq-hidrogenband, maxf=hidrogenfreq+hidrogenband, nSigma = numSigma;
	char outlierFound;
	float signal[size][8];
	
	for(i = 0; i < size; i++)
		{
		if(!dataset[i].flagBAD)
			{
			for(j=0; j<8; j++) signal[i][j] = 0;
			N = 0; 
			for(k = lowchan; k < highchan; k++)
				{
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
			for(j=0; j<8; j++) signal[i][j] = signal[i][j]/N;
			}
		}
		
	do
		{
		N = 0; 
		for(j=0; j<8; j++){mean[j] = 0; sigma[j] = 0;}
		for(i = 1; i < size - 1; i++)
			{
			if(dataset[i].flagBAD && dataset[i-1].flagBAD && dataset[i+1].flagBAD) continue;
			N++;
			for(j=0; j<8; j++)
				{
				diff = signal[i+1][j] - 2*signal[i][j] + signal[i-1][j];
				delta = diff - mean[j];
				mean[j] += delta/N;
				sigma[j] += delta*(diff - mean[j]);
				}
			}
		for(j=0; j<8; j++){if(N>1) sigma[j] = sqrt(sigma[j]/(N-1)); else sigma[j] = sqrt(sigma[j]);}
		outlierFound = 0;
		for(i = 1; i < size - 1; i++)
			{
			if(field[0] == 'N' && field[1] == '1' && ((dataset[i].RA>82.9333 && dataset[i].RA<84.0666) || (dataset[i].RA>68.95 && dataset[i].RA<69.6))) continue; // Just for N1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(dataset[i].flagBAD && dataset[i-1].flagBAD && dataset[i+1].flagBAD) continue;
			if( fabs(signal[i+1][0] - 2*signal[i][0] + signal[i-1][0] - mean[0]) > numSigma*sigma[0] ||
				fabs(signal[i+1][1] - 2*signal[i][1] + signal[i-1][1] - mean[1]) > numSigma*sigma[1] ||
				fabs(signal[i+1][2] - 2*signal[i][2] + signal[i-1][2] - mean[2]) > numSigma*sigma[2] ||
				fabs(signal[i+1][3] - 2*signal[i][3] + signal[i-1][3] - mean[3]) > numSigma*sigma[3] ||
				fabs(signal[i+1][4] - 2*signal[i][4] + signal[i-1][4] - mean[4]) > numSigma*sigma[4] ||
				fabs(signal[i+1][5] - 2*signal[i][5] + signal[i-1][5] - mean[5]) > numSigma*sigma[5] ||
				fabs(signal[i+1][6] - 2*signal[i][6] + signal[i-1][6] - mean[6]) > numSigma*sigma[6] ||
				fabs(signal[i+1][7] - 2*signal[i][7] + signal[i-1][7] - mean[7]) > numSigma*sigma[7])
				{
				dataset[i-1].flagBAD = 1;
				dataset[i].flagBAD = 1;
				dataset[i+1].flagBAD = 1;
				for(k = lowchan; k < highchan; k++)
					{
						{
						dataset[i-1].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i+1].flagRFI[k] |= RFI_OUTOFBAND;
						}
					}
				outlierFound = 1;
				}
			}
		}
	while(outlierFound);
}

void rfi_detection_time_domain3(const char *field, SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, float hidrogenfreq, float hidrogenband, float freq[])
{
// difference
	int i, j, k, N;
	float mean[10], sigma[10], diff, delta, minf=hidrogenfreq-hidrogenband, maxf=hidrogenfreq+hidrogenband, nSigma = numSigma;
	char outlierFound;
	float signal[size][10];
	
	for(i = 0; i < size; i++)
		{
		if(!dataset[i].flagBAD)
			{
			for(j=0; j<8; j++) signal[i][j] = 0;
			N = 0; 
			for(k = lowchan; k < highchan; k++)
				{
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
			for(j=0; j<8; j++) signal[i][j] = signal[i][j]/N;
			signal[i][8] = signal[i][0] - signal[i][1];
			signal[i][9] = signal[i][2] - signal[i][3];
			}
		}
		
	do
		{
		N = 0; 
		for(j=0; j<10; j++){mean[j] = 0; sigma[j] = 0;}
		for(i = 1; i < size - 1; i++)
			{
			if(dataset[i].flagBAD && dataset[i-1].flagBAD && dataset[i+1].flagBAD) continue;
			N++;
			for(j=0; j<10; j++)
				{
				diff = signal[i+1][j] - 2*signal[i][j] + signal[i-1][j];
				delta = diff - mean[j];
				mean[j] += delta/N;
				sigma[j] += delta*(diff - mean[j]);
				}
			}
		for(j=0; j<10; j++){if(N>1) sigma[j] = sqrt(sigma[j]/(N-1)); else sigma[j] = sqrt(sigma[j]);}
		outlierFound = 0;
		for(i = 1; i < size - 1; i++)
			{
			if(field[0] == 'N' && field[1] == '1' && ((dataset[i].RA>82.9333 && dataset[i].RA<84.0666) || (dataset[i].RA>68.95 && dataset[i].RA<69.6))) continue; // Just for N1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(dataset[i].flagBAD && dataset[i-1].flagBAD && dataset[i+1].flagBAD) continue;
			if( fabs(signal[i+1][8] - 2*signal[i][8] + signal[i-1][8] - mean[8]) > numSigma*sigma[8] ||
				fabs(signal[i+1][9] - 2*signal[i][9] + signal[i-1][9] - mean[9]) > numSigma*sigma[9] ||
				fabs(signal[i+1][4] - 2*signal[i][4] + signal[i-1][4] - mean[4]) > numSigma*sigma[4] ||
				fabs(signal[i+1][5] - 2*signal[i][5] + signal[i-1][5] - mean[5]) > numSigma*sigma[5] ||
				fabs(signal[i+1][6] - 2*signal[i][6] + signal[i-1][6] - mean[6]) > numSigma*sigma[6] ||
				fabs(signal[i+1][7] - 2*signal[i][7] + signal[i-1][7] - mean[7]) > numSigma*sigma[7])
				{
				dataset[i-1].flagBAD = 1;
				dataset[i].flagBAD = 1;
				dataset[i+1].flagBAD = 1;
				for(k = lowchan; k < highchan; k++)
					{
						{
						dataset[i-1].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i].flagRFI[k] |= RFI_OUTOFBAND;
						dataset[i+1].flagRFI[k] |= RFI_OUTOFBAND;
						}
					}
				outlierFound = 1;
				}
			}
		}
	while(outlierFound);
}

void outofbandrfi_ann(SpecRecord dataset[], int size, int lowchan, int highchan)
{
	int i;

	FILE * file;
	file = fopen("outofbandrfi.ann", "w");
	fprintf(file, "COLOUR RED\n");

	for (i = 0; i < size; i++)
		{
		if(dataset[i].flagRFI[lowchan] && RFI_OUTOFBAND) fprintf(file, "LINE W %i %f %i %f\n", i, dataset[i].RA, 30, dataset[i].AST);
		}	
	fclose(file);
/*
	file = fopen("rfi_rd.ann", "w");
	for (i = 0; i < size; i++)
		{
		if(dataset[i].flagRFI[lowchan] && RFI_OUTOFBAND) fprintf(file, "%f %f %f\n", dataset[i].RA, dataset[i].DEC, dataset[i].AST);
		}	
	fclose(file);
*/
}

void rfi_ann(SpecRecord dataset[], int size, int lowchan, int highchan, float freq[])
{
	int i, k;

	FILE * file;
	file = fopen("rfi.ann", "w");
	
	for(i = 0; i < size; i++)
		{
		for(k = lowchan; k < highchan; k++)
			{
			if(dataset[i].flagRFI[k] != RFI_NONE) fprintf(file, "%f %f\n", freq[k], dataset[i].AST); 
			}
		}
	
	fclose(file);
}

#include "rfi.h"
#include <stdlib.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define SQR(X) ((X)*(X))

void compute_stats_on_stats(const PolStatistics * stats, int size, PolStatistics * pStatOut)
{
	int i;
	int count = 0;

	//working sum values on the stack
	double sumOffXX, sumOffYY;
	double sumOnXX, sumOnYY;

	//cached mean values 
	double meanOffXX, meanOffYY;
	double meanOnXX, meanOnYY;

	/* Calculate the means */

	sumOffXX = sumOffYY = 0;
	sumOnXX = sumOnYY = 0;
	for (i=0; i<size; i++) 
	{
		const PolStatistics * pStat = &(stats[i]);

		sumOffXX += pStat->sigmaOffXX;
		sumOffYY += pStat->sigmaOffYY;
		sumOnXX += pStat->sigmaOnXX;
		sumOnYY += pStat->sigmaOnYY;

		count+=1;
	}

	meanOffXX = sumOffXX / count;
	meanOffYY = sumOffYY / count;
	meanOnXX = sumOnXX / count;
	meanOnYY = sumOnYY / count;

	pStatOut->meanOffXX = meanOffXX;
	pStatOut->meanOffYY = meanOffYY;
	pStatOut->meanOnXX = meanOnXX;
	pStatOut->meanOnYY = meanOnYY;

	/* Calculate the sigmas */

	sumOffXX = sumOffYY = 0;
	sumOnXX = sumOnYY = 0;
	for (i=0; i<size; i++) 
	{
		const PolStatistics * pStat = &(stats[i]);
		sumOffXX += SQR(pStat->sigmaOffXX - meanOffXX);
		sumOffYY += SQR(pStat->sigmaOffYY - meanOffYY);
		sumOnXX += SQR(pStat->sigmaOnXX - meanOnXX);
		sumOnYY += SQR(pStat->sigmaOnYY - meanOnYY);

		count+=1;
	}

	pStatOut->sigmaOffXX = sqrt(sumOffXX / (count-1));
	pStatOut->sigmaOffYY = sqrt(sumOffYY / (count-1));
	pStatOut->sigmaOnXX = sqrt(sumOnXX / (count-1));
	pStatOut->sigmaOnYY = sqrt(sumOnYY / (count-1));

}


//compute the means using the non excluded channels
//sets the mean values in pStat using diffs in pStat
//only channels marked RFI_NONE in pRec are considered
void compute_stats_on_diffs(const SpecRecord * pRec, const PolDifferences * pDiff, PolStatistics * pStat, int lowchan, int highchan)
{
	int i;
	int count;

	//working sum values on the stack
	double sumOffXX, sumOffYY;
	double sumOnXX, sumOnYY;

	//cached mean values 
	double meanOffXX, meanOffYY;
	double meanOnXX, meanOnYY;

	/* Calculate the means */

	count = 0;
	sumOffXX = sumOffYY = 0;
	sumOnXX = sumOnYY = 0;
	for (i=lowchan; i<highchan; i++) 
	{
		if (pRec->flagRFI[i] == RFI_NONE)
		{
			sumOffXX += pDiff->OffXX[i];
			sumOffYY += pDiff->OffYY[i];
			sumOnXX += pDiff->OnXX[i];
			sumOnYY += pDiff->OnYY[i];

			count+=1;
		}
	}

	meanOffXX = sumOffXX / count;
	meanOffYY = sumOffYY / count;
	meanOnXX = sumOnXX / count;
	meanOnYY = sumOnYY / count;

	pStat->meanOffXX = meanOffXX;
	pStat->meanOffYY = meanOffYY;
	pStat->meanOnXX = meanOnXX;
	pStat->meanOnYY = meanOnYY;

	/* Calculate the sigmas */

	count = 0;
	sumOffXX = sumOffYY = 0;
	sumOnXX = sumOnYY = 0;
	for (i=lowchan; i<highchan; i++) 
	{
		if (pRec->flagRFI[i] == RFI_NONE) 
		{
			sumOffXX += SQR(pDiff->OffXX[i] - meanOffXX);
			sumOffYY += SQR(pDiff->OffYY[i] - meanOffYY);
			sumOnXX += SQR(pDiff->OnXX[i] - meanOnXX);
			sumOnYY += SQR(pDiff->OnYY[i] - meanOnYY);

			count+=1;
		}
	}

	pStat->sigmaOffXX = sqrt(sumOffXX / (count-1));
	pStat->sigmaOffYY = sqrt(sumOffYY / (count-1));
	pStat->sigmaOnXX = sqrt(sumOnXX / (count-1));
	pStat->sigmaOnYY = sqrt(sumOnYY / (count-1));

}


//compute the diffs for all channels
//sets the diff values in pStat using diffs and means in pStat
//differences are computed as "next channel minus the current channel"
void compute_diffs(const SpecRecord * pRec, PolDifferences * pDiff, int lowchan, int highchan)
{
	int i;

	//pre-compute the diffs for all adjacent channels for efficiency
	for (i=lowchan; i<highchan; i++) 
	{
		pDiff->OffXX[i] = pRec->caloff.xx[i+1] - pRec->caloff.xx[i];
		pDiff->OffYY[i] = pRec->caloff.yy[i+1] - pRec->caloff.yy[i];
		pDiff->OnXX[i] = pRec->calon.xx[i+1] - pRec->calon.xx[i];
		pDiff->OnYY[i] = pRec->calon.yy[i+1] - pRec->calon.yy[i];
	}
}

void print_means(FILE * file, const SpecRecord * pRec, const PolStatistics * pStat)
{
	fprintf(file, "%8.2f %8.6f %8.6f %8.6f %8.6f \n", 
			pRec->AST, pStat->meanOffXX, pStat->meanOffYY,
			pStat->meanOnXX, pStat->meanOnYY);
}


void print_sigmas(FILE * file, const SpecRecord * pRec, const PolStatistics * pStat)
{
	fprintf(file, "%8.2f %8.6f %8.6f %8.6f %8.6f\n", 
		pRec->AST, pStat->sigmaOffXX, pStat->sigmaOffYY,
		pStat->sigmaOnXX, pStat->sigmaOnYY);
}

void print_diffs(FILE * file, const SpecRecord * pRec, const PolDifferences * pDiff, int lowchan, int highchan)
{
	int i;
	for (i=lowchan; i<highchan; i++) {
		fprintf(file, "%8.2f %8.6f %8.6f %8.6f %8.6f \n", 
			pRec->AST, pDiff->OffXX[i], pDiff->OffYY[i],
			pDiff->OnXX[i], pDiff->OnYY[i]);
	}
}

/* data on either side of an RFI detection may not be detected as RFI, but is suspect 
 * nonetheless.  This routine broadens the the detections across the specified span.
 */
void rfi_spanning(SpecRecord dataset[], int size, int lowchan, int highchan, int span)
{
	int n, i, chan;

	//iterate over each time step and span RFI across timesteps
	for (chan=lowchan; chan<highchan; chan++) 
	{
		for (n=0; n<size; n++) 
		{
			if (dataset[n].flagRFI[chan]) {
				for (i=n-span/2; i<n+span/2; i++) {
					if (i<0 || i>=size) continue;
					dataset[i].flagRFI[chan] |= RFI_SPAN;
				}
				n += span/2+1;
			}
		}
	}
}

//counts the number of datapoints for each channel from highchan to lowchan.  
//If more than rfiTolerance percentage of the points are RFI,
//RFI flagged datapoints, then the entire channel is flagged as RFI.
void rfi_blanking(SpecRecord dataset[], int numRecords, int lowchan, int highchan, int rfiTolerance)
{
	int i;
	int chan;
	int rfiCount;
	int rfiPercent;

	for (chan = lowchan; chan < highchan; chan++) 
	{
		rfiCount = 0;
		for (i=0; i<numRecords; i++) {
			if (dataset[i].flagRFI[chan] != RFI_NONE) {
				rfiCount++;
			}
		}

		rfiPercent = rfiCount*100/numRecords;
		if (rfiPercent > rfiTolerance) {
			printf("Blanking channel %i since %i%% of the points are RFI\n", chan, rfiPercent);
			for (i=0; i<numRecords; i++) {
				dataset[i].flagRFI[chan] = RFI_OVERLIMIT;
			}			
		}
	}
}


// freq - an array of frequencies that correspond to the channels
void rfi_write(SpecRecord dataset[], int size, int lowchan, int highchan, float * freq)
{
	int n, chan;
	FILE * rfifile;

	rfifile = fopen("rfi.dat","w");
	fprintf(rfifile, "# chan freq RA DEC AST OffXX OffYY OnXX OnYY code\n");

	for (n=0; n<size; n++) 
	{
		SpecRecord * pRec = &(dataset[n]);

		//now that the channels are marked as RFI, write them out
		for (chan=lowchan; chan<highchan; chan++) { 
		 	if (pRec->flagRFI[chan] != RFI_NONE) {
				fprintf(rfifile,
					"%3i %9.6f %9.6f %9.6f %8.2f " 
					"%8.6f %8.6f %8.6f %8.6f "
					"%0#6x\n", 
					chan, freq[chan], pRec->RA, pRec->DEC, pRec->AST,  
					pRec->caloff.xx[chan], pRec->caloff.yy[chan],
					pRec->calon.xx[chan], pRec->calon.yy[chan],
					pRec->flagRFI[chan]);
			}
		}
	}

	fclose (rfifile);
}

// dataset - the channel data for a particular time step
// lowchan - the lowest channel number to do computations for
// highchan - the highest channel number to do computations for
// numSigma - the number of sigma away from the mean the difference should be rejected
// ignore?_[low|high] - two ranges of channels to ignore from an RFI perspective
void rfi_detection(SpecRecord dataset[], int size, int lowchan, int highchan, float numSigma, 
	float numSigmaThresh, int ignoreA_low, int ignoreA_high, 
	int ignoreB_low, int ignoreB_high)
{
	FILE * finalmeanfile;
	FILE * finalsigmafile;
	FILE * sigmathreshfile;
	
	PolStatistics sigmaStat;
	PolStatistics * stats;
	PolDifferences * diffs;
	int i;

	float sigmaThreshOffXX, sigmaThreshOffYY;
	float sigmaThreshOnXX, sigmaThreshOnYY;

	finalmeanfile = fopen("finalmean.dat", "w");
	fprintf(finalmeanfile, "# date offXX offYY onXX onYY\n");

	finalsigmafile = fopen("finalsigma.dat", "w");
	fprintf(finalsigmafile, "# date offXX offYY onXX onYY\n");

	sigmathreshfile = fopen("sigmathresh.dat", "w");
	fprintf(sigmathreshfile, "# offXX offYY onXX onYY\n");

	printf("Requesting malloc for %u bytes of memory\n",sizeof(PolStatistics)*size);
	stats = (PolStatistics *)malloc(sizeof(PolStatistics) * size);
	if (stats == NULL) {
		printf("ERROR: malloc failed in rfi_detection() !\n");
	}
	printf("Requesting malloc for %u bytes of memory\n",sizeof(PolDifferences)*size);
	diffs = (PolDifferences *)malloc(sizeof(PolDifferences) * size);
	if (diffs == NULL) {
		printf("ERROR: malloc failed in rfi_detection() !\n");
	}

	//iterate over each time step and calculate a final sigma over the
	//frequency domain
	for (i=0; i<size; i++) 
	{
		int chan;
		char outlierFound;

		SpecRecord * pRec = &(dataset[i]);
		PolStatistics * pStat = &(stats[i]);
		PolDifferences * pDiff = &(diffs[i]);

		//skip any timesteps flagged as bad
		if (pRec->flagBAD) continue;

		compute_diffs(pRec, pDiff, lowchan, highchan);

		//clear the channel flags
		for (chan=lowchan; chan<highchan; chan++) { 
			pRec->flagRFI[chan] = RFI_NONE;
		}	


		//now iterate over the channels 
		//compute sigma for channels not flagged as excluded
		//apply sigma to exclude any channels
		//repeat while we still find outliers
		do {

			//calculate means and sigmas from not flagged channels
			compute_stats_on_diffs(pRec, pDiff, pStat, lowchan, highchan);

			//we now have sigma values for a particular timestep
			//now iterate over the frequencies and mark outliers
			outlierFound = FALSE;
			for (chan=lowchan; chan<highchan; chan++) 
			{
				if ((chan >= ignoreA_low && chan <= ignoreA_high) ||
					(chan >= ignoreB_low && chan <= ignoreB_high)) {
					continue;
			}
				if (pRec->flagRFI[chan] == RFI_NONE)
				{
					if ( (fabs(pStat->meanOffXX - pDiff->OffXX[chan]) > (numSigma * pStat->sigmaOffXX)) 
						|| (fabs(pStat->meanOffYY - pDiff->OffYY[chan]) > (numSigma * pStat->sigmaOffYY)) 
						|| (fabs(pStat->meanOnXX - pDiff->OnXX[chan]) > (numSigma * pStat->sigmaOnXX)) 
						|| (fabs(pStat->meanOnYY - pDiff->OnYY[chan]) > (numSigma * pStat->sigmaOnYY)) ) {
						pRec->flagRFI[chan] = 1;
						pRec->flagRFI[chan+1] = 1;
						outlierFound = TRUE;
					}
				}
			}
		} while (outlierFound);
	
	}	//we now have a final sigma

	//determine the distribution of the sigmas and set thresholds
	compute_stats_on_stats(stats, size, &sigmaStat);	
	sigmaThreshOffXX = sigmaStat.sigmaOffXX * numSigmaThresh + sigmaStat.meanOffXX;
	sigmaThreshOffYY = sigmaStat.sigmaOffYY * numSigmaThresh + sigmaStat.meanOffYY;
	sigmaThreshOnXX = sigmaStat.sigmaOnXX * numSigmaThresh + sigmaStat.meanOnXX;
	sigmaThreshOnYY = sigmaStat.sigmaOnYY * numSigmaThresh + sigmaStat.meanOnYY;
	fprintf(stdout, "Sigma of Sigmas: %8.6f %8.6f %8.6f %8.6f \n",
			sigmaStat.sigmaOffXX, sigmaStat.sigmaOffYY,
			sigmaStat.sigmaOnXX, sigmaStat.sigmaOnYY);
	fprintf(stdout, "Mean of Sigmas: %8.6f %8.6f %8.6f %8.6f \n",
			sigmaStat.meanOffXX, sigmaStat.meanOffYY,
			sigmaStat.meanOnXX, sigmaStat.meanOnYY);
	fprintf(sigmathreshfile, "%8.6f %8.6f %8.6f %8.6f \n",
			sigmaThreshOffXX, sigmaThreshOffYY,
			sigmaThreshOnXX, sigmaThreshOnYY);
	
	//iterate over each time step and mark RFI in the channels
	for (i=0; i<size; i++) 
	{
		int chan;

		SpecRecord * pRec = &(dataset[i]);
		PolStatistics * pStat = &(stats[i]);
		PolDifferences * pDiff = &(diffs[i]);

		//skip any timesteps flagged as bad
		if (pRec->flagBAD) continue;

		//clear the channel flags
		for (chan=lowchan; chan<highchan; chan++) { 
			pRec->flagRFI[chan] = RFI_NONE;
		}	


		//do one last pass over the channels and mark the flags using the 
		//final values of sigma and the means
		for (chan=lowchan; chan<highchan; chan++) 
		{
			//skip over channels that are to be excluded
			if ((chan >= ignoreA_low && chan <= ignoreA_high) ||
				(chan >= ignoreB_low && chan <= ignoreB_high)) {
				continue;
			}

			//if the sigma is very large, then the sigma clipping will not work
			//so check for the sigma tolerance, and mark all channels as RFI
			//if the threshold is exceeded.
			if ((pStat->sigmaOffXX > sigmaThreshOffXX)
				|| (pStat->sigmaOffYY > sigmaThreshOffYY)
				|| (pStat->sigmaOnXX > sigmaThreshOnXX)
				|| (pStat->sigmaOnYY > sigmaThreshOnYY)) {
				pRec->flagRFI[chan] = RFI_SIGMA_EXCEEDED;
				continue;
			}

			//mark the channel and the next one as RFI if the difference is an outlier
			//flagRFI bitmask is updates with all reasons why it was RFI
			if (fabs(pStat->meanOffXX - pDiff->OffXX[chan]) > (numSigma * pStat->sigmaOffXX)) {
				pRec->flagRFI[chan] |= RFI_CALOFF_XX;
				pRec->flagRFI[chan+1] |= RFI_CALOFF_XX;
			}
			if (fabs(pStat->meanOffYY - pDiff->OffYY[chan]) > (numSigma * pStat->sigmaOffYY)) {
				pRec->flagRFI[chan] |= RFI_CALOFF_YY;
				pRec->flagRFI[chan+1] |= RFI_CALOFF_YY;
			}
			if (fabs(pStat->meanOnYY - pDiff->OnYY[chan]) > (numSigma * pStat->sigmaOnYY)) {
				pRec->flagRFI[chan] |= RFI_CALON_YY;
				pRec->flagRFI[chan+1] |= RFI_CALON_YY;
			}
			if (fabs(pStat->meanOnXX - pDiff->OnXX[chan]) > (numSigma * pStat->sigmaOnXX)) {
				pRec->flagRFI[chan] |= RFI_CALON_XX;
				pRec->flagRFI[chan+1] |= RFI_CALON_XX;
			}
		}


		//print the final sigmas and means for this time step
		print_sigmas(finalsigmafile, pRec, pStat);
		print_means(finalmeanfile, pRec, pStat);

	} //endfor timestep

	free(stats);
	free(diffs);
	fclose(finalmeanfile);
	fclose(finalsigmafile);
	fclose(sigmathreshfile);
}

//call this on unsmoothed raw data
void aerostat_rfi_blanking(SpecRecord dataset[], int size, int lowchan, int highchan)
{
	int i, j;
	double *offxx;
	float *diff1, *diff2;
	int chan;

	FILE * file;
	file = fopen("outofbandrfi.ann", "w");
	fprintf(file, "COLOUR RED\n");

	//average up the channels to reduce noise as these signals are broadband
	//they peak on the middle channel (128) which could also be used 
	printf("Requesting malloc for %u bytes of memory\n",sizeof(double)*size);
	offxx = calloc(size, sizeof(double));
	if (offxx == NULL) {
		printf("ERROR: malloc failed in aerostat_rfi_blanking() !\n");
	}
	printf("Requesting malloc for %u bytes of memory\n",sizeof(float)*size);	
	diff1 = calloc(size, sizeof(float));
	if (diff1 == NULL) {
		printf("ERROR: malloc failed in aerostat_rfi_blanking() !\n");
	}
	printf("Requesting malloc for %u bytes of memory\n",sizeof(float)*size);	
	diff2 = calloc(size, sizeof(float));
	if (diff2 == NULL) {
		printf("ERROR: malloc failed in aerostat_rfi_blanking() !\n");
	}

	for (i=0; i<size; i++)
	{
		int count = 0;
		for (chan = lowchan; chan<highchan; chan++) {
			if (dataset[i].flagRFI[chan]) {
				continue;
			}
			offxx[i] +=  dataset[i].caloff.xx[chan] - dataset[i].caloff.yy[chan];
			count++;
		}
		offxx[i] /= count;
	}
	

	//compute first differences
	for (i=0; i<size-1; i++) {
		diff1[i] = offxx[i] - offxx[i+1];
	}
	//compute second differences
	for (i=0; i<size-2; i++) {
		diff2[i] = diff1[i] - diff1[i+1];
	}

	//print out the differences
	{
		FILE * difffile = fopen("diff.dat", "w");
		for (i=0; i<size-2; i++) {
			fprintf(difffile, "%f %f %f %f %f\n", dataset[i].AST, offxx[i], offxx[i+1], diff1[i], diff2[i]);
		}
		fclose(difffile);
	}

	for (i=0; i<size-2; i++)
	{
		if (dataset[i].flagBAD || dataset[i+1].flagBAD) {
			continue;
		}

//		if (diff2[i] > 0.07 || diff2[i] < -0.07) 
		if ((diff2[i] > 0.07 || diff2[i] < -0.07) && i>0) //ssg fixed i > 0 
		{

			printf("Outofband rfi with diff %f %f\n", diff1[i], diff2[i]);
			fprintf(file, "LINE W %i %f %i %f\n", 0, dataset[i+1].AST, 30, dataset[i+1].AST);
			for (j=i-1; j<i+12*5 && j<size; j++) 
			{
				for (chan=0; chan<MAX_CHANNELS; chan++) 
				{
					dataset[j].flagBAD = 1;
					dataset[j].flagRFI[chan] |= RFI_OUTOFBAND;
				}
			}
			i = j-1; //start a little early for the next round
		}
	}

	fclose(file); 
	free(offxx);
	free(diff1);
	free(diff2);
}




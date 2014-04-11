#ifndef _STATS_H
#define _STATS_H

double compute_mean(double data[], int start, int end);
double compute_sigma(double data[], int size, double mean);
int reject_outliers(double data[], int size, float nsigma);
double compute_clean_mean(double data[], int size, float nsigma);

#endif //_STATS_H


/*
 * tcal_on_fly.h
 *
 *  Created on: 2013-04-11
 *      Author: csmith
 */
#include "spec.h"

#ifndef TCAL_ON_FLY_H_
#define TCAL_ON_FLY_H_

void compute_tcal(SpecRecord dataset[], int size, int lowchan, int highchan, float hif, float hiband, float freq[], int *badchannels, float tcalxx[],float tcalyy[],int r,int window);
void norm_one_tcal(int lowchan,int highchan,int *badchannels,float* tcalxx,float* tcalyy);
void write_tcal(float* tcalxx,float* tcalyy, int r,int low, int high);
#endif /* TCAL_ON_FLY_H_ */

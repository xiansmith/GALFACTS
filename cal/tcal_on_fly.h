/*
 * tcal_on_fly.h
 *
 *  Created on: 2013-04-11
 *      Author: csmith
 */
#include "spec.h"

#ifndef TCAL_ON_FLY_H_
#define TCAL_ON_FLY_H_

void compute_tcal(SpecRecord dataset[], int size, int lowchan, int highchan, float hif, float hiband, float freq[], int *badchannels, float tcalxx[],float tcalyy[]);



#endif /* TCAL_ON_FLY_H_ */

#include "common.h"
#include "fluxdata.h"

#ifndef _DECDEPENDENCE_H
#define _DECDEPENDENCE_H

void calculate_dec_dependence(FluxWappData * wappdata, float decmin, float decmax, float decgrain, int chan);
void remove_dec_dependence(FluxWappData * wappdata, float decmin, float decmax, int chan);
void average_dec_dependence(FluxWappData * wappdata, int lowchan, int highchan);

#endif


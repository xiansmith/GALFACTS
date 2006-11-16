#ifndef _GRID_H
#define _GRID_H

#include "common.h"
#include "map.h"


#define GRID_UP_SCANS 0x1
#define GRID_DOWN_SCANS 0x2
#define GRID_ALL_SCANS GRID_UP_SCANS|GRID_DOWN_SCANS

void init_psf_lookup_table(int size, float maxVal);

void grid_data(const FluxWappData *wappdata, const MapMetaData *md, float * dataI, float * dataQ, float * dataU, float * dataV, float *weight);


#endif


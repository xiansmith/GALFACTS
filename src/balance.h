#ifndef _BALANCE_H
#define _BALANCE_H

#include "common.h"
#include "fluxdata.h"
#include "map.h"

void balance_data(FluxWappData * wappdata, MapMetaData *md, int day_order, int scan_order, float loop_gain, float loop_epsilon, int progress);


#endif 

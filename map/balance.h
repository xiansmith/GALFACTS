#ifndef _BALANCE_H
#define _BALANCE_H

#include "common.h"
#include "fluxdata.h"
#include "map.h"
enum ShowProgress {NONE,I,Q,U,V,ALL};
void balance_data(FluxWappData * wappdata, int day_order, int scan_order, float loop_gain, float loop_epsilon, int bw_order);
void find_intersections(FluxWappData *wappdata);
void write_balance_data(FluxWappData * wappdata, int day_order, int scan_order, float loop_gain, float loop_epsilon, int order);
void read_apply_balance_data( FluxWappData * wappdata, int day_order, int scan_order, float loop_gain, float loop_epsilon, int bw_order );
//float day_weaveI(ScanDayData *daydata, float loop_gain, int day, int order, BWcoeff* coeffs );
//float day_weaveQ(ScanDayData *daydata, float loop_gain, int day, int order, BWcoeff* coeffs );
//float day_weaveU(ScanDayData *daydata, float loop_gain, int day, int order, BWcoeff* coeffs );
//float day_weaveV(ScanDayData *daydata, float loop_gain, int day, int order, BWcoeff* coeffs );
#endif

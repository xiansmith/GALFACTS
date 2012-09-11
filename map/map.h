#ifndef _MAP_H
#define _MAP_H

enum MapType {MAPTYPE_AVERAGE, MAP_TYPE_CUBE};

typedef struct {
	enum MapType type;	
	float ramin;
	float ramax;
	float decmin;
	float decmax;
	float RAcen;
	float DECcen;
	float RArange;
	float DECrange;
	int lowchan;
	int highchan;
	float cellsize;             // cell size of map in degrees
	int n1, n2, n3;
	int patch;
	int beamwidths;
	int gridtype;
	const char * title;
	float fcen;
	float fstart;
	float fwhm;
	float df;
	int avg;
	int avg_lowchan;
	int avg_highchan;
	int cavg;
	int band;
} MapMetaData;

void write_fits_maps(const char * wapp, MapMetaData *md, float dataI[], float dataQ[], float dataU[], float dataV[]);

void start_fits_cubes(const char * wapp, MapMetaData *md);

void write_fits_planes(float dataI[], float dataQ[], float dataU[], float dataV[], float dataW[]);

void finish_fits_cubes(void);

#endif

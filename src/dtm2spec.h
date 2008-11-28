#include "galfactsLib.h"
#define MAIN_HEADER 2880
#define BINTABLE_HEADER 20*2880
#define RAW_CHANNELS 4096
#define DUMPS_PER_ROW 3
#define LINELEN 80
typedef struct
{
	int A[RAW_CHANNELS];
	int B[RAW_CHANNELS];
	unsigned int U[RAW_CHANNELS];
	unsigned int V[RAW_CHANNELS];
}GFLIB_DATA;

typedef struct
{
	GFLIB_DATA dataon[DUMPS_PER_ROW];
	GFLIB_DATA dataoff[DUMPS_PER_ROW];
	GFLIB_PDEVSTAT staton[DUMPS_PER_ROW];
	GFLIB_PDEVSTAT statoff[DUMPS_PER_ROW];
	char object[16];
	double cf; // center frequency
	double fdelt; //delta frequency
	double RA;
	double DEC;
	double pol; //?? not used
	double AST;
	double tdelt; //time delta ?? integration time per point ??
	char fill[768]; //just read in the rest of the row, values currently not used for anything
}GFLIB_ROW;

static inline void cnvrt_end_sint(short int *x);
static inline void cnvrt_end_int(int *x);
static inline void cnvrt_end_db(double *x);


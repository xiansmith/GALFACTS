#include "common.h"
#include "spec.h"

#define BADFILE_LINE_LEN 100

/*
 * Reads the provided badfilename and marks bad datapoints as listed.
 * Format of the bad file is as follows:
 * # at the start of line for comment
 * [type] [low] [high]
 * type - one of RA, DEC, AST
 * low high - are integer or floating point numbers that delimit a range.  INF
 * can be used for infinity.
 * blank lines are allowed
 * lines are limited to BADFILE_LINE_LEN characters.
 */
void mark_bad_datapoints(const char * badfilename, SpecRecord dataset[], int numRecords);

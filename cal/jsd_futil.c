#include "jsd_futil.h"
//----------------------------------------------------------------------------------------------------------
int jsd_line_count(FILE * file)
{
	long pos;
	int count;

	if (file == NULL) {
		return -1;
	}

	pos = ftell(file);
	count = 0;
	while (!feof(file)) {
		if (fgetc(file) == '\n') count ++;
	}

	fseek(file, pos, SEEK_SET);
	clearerr(file);

	return count;
}
//----------------------------------------------------------------------------------------------------------

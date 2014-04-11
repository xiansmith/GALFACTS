
#include "jsd_futil.h"

int main(int argc, char *argv[])
{
	int count1;
	FILE * file1;
	
	file1 = fopen("file1", "r");
	count1 = jsd_futil_line_count(file1);
	printf("count1: %i\n", count1);
	if (file1 != NULL) {
		fclose(file1);
	}

	return EXIT_SUCCESS;
}

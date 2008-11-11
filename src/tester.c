#include <stdlib.h>
#include <stdio.h>

int main(int argc, char ** argv)
{
	int i,j;
	int segments = 4;
	int num = 100;

	for (i=0; i<segments; i++) 
	{
		for (j=0; j<segments; j++) {
			int start = (num*j)/segments;
			int end = (num*(j+1))/segments - 1;
			printf("%i %i\n", start, end);
		}
	}
	return EXIT_SUCCESS;
}


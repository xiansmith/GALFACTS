#include <stdio.h>

void main(int argc, char *argv[])
{
	char * filename = argv[1];
	int order = atoi(argv[2]);
	FILE *in = fopen(filename,"r");

	int i;
	for(i = 0 ;i < order+1;i++)
	{
		float C;
		fscanf(in,"%f",&C);
		printf("%f*x**%d +",C,i);
	}
	printf("\n");
	fclose(in);
}

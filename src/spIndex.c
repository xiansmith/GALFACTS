#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include <math.h>
#include "programs/fitsio.h"
#include "jsd/jsd_futil.h"
#include "jsd/jsd_fit.h"

//#define IS_BLANK_PIXEL(X) ( (X>(BLANK_PIXEL-BLANK_THRESH)) && (X<(BLANK_PIXEL+BLANK_THRESH)) )

static void print_usage(char * argv[])
{
	printf("\n");
	printf("Usage: %s <in.fits> <in.cat>\n", argv[0]);
}

int main(int argc, char *argv[]) 
{

	int i, j, k, n;
	char *infilename, *outfilename, *incubename;
	FILE *incat,*incube,*outfile;
	header_param_list inhpar;
	float *plane_data;
	double **fitdata;
	int plane_size;
	int num_planes;
	int num_comp;

	if (argc != 3) {
		print_usage(argv);
		return EXIT_FAILURE;
	}
	incubename = argv[1];
	infilename = argv[2];

	incube = fopen(incubename,"r");
	readfits_header(incube,&inhpar);

	plane_size = inhpar.naxis[0] * inhpar.naxis[1];
	num_planes = inhpar.naxis[2];
	plane_data = malloc(sizeof(float)*plane_size);
	printf("num_planes: %i\n", num_planes);
	printf("plane_size: %i\n", plane_size);

	incat = fopen(infilename,"r");
	num_comp = jsd_line_count(incat);
	printf("num_comp: %i\n", num_comp);
	
	fitdata = (double **)malloc (sizeof(double*) * num_comp);
	
	for(i =0;i<num_comp;i++)
	{
		fitdata[i] = (double *)malloc(sizeof(double)*num_planes);
	}
	
	int *coords[2];
	coords[0] = malloc(sizeof(int)* num_comp);
	coords[1] = malloc(sizeof(int)* num_comp);
	double *logfreq = malloc(sizeof(double *)*num_planes);

	for(i =0;i<num_comp;i++)
	{
		float RA,DEC,temp;
		fscanf(incat,"%f %f %f",&RA,&DEC,&temp);
		coords[0][i] = (int)(inhpar.crpix[0]+(RA-inhpar.crval[0])/inhpar.cdelt[0]+0.5)-1;
		coords[1][i] = (int)(inhpar.crpix[1]+(DEC-inhpar.crval[1])/inhpar.cdelt[1]+0.5)-1;
		//printf("%f %f %d %d\n",RA,DEC,coords[0][i],coords[1][i]);
	}

	for(i=0;i<num_planes;i++)
	{
		readfits_plane(incube,plane_data,&inhpar);
		logfreq[i] = log(inhpar.crval[2]+(i+1-inhpar.crpix[2])*inhpar.cdelt[2]);
		//printf("%f ",logfreq[i]);
		for(j=0;j<num_comp;j++)
		{
			fitdata[j][i] = log(plane_data[coords[0][j]+coords[1][j]*inhpar.naxis[0]]);
		//	if(!j)	
		//		printf("%f\n",fitdata[j][i]);
		}
	}
	//exit(1);
	for(i =0;i<num_comp;i++)
	{
		double C[2],chisq;
		jsd_poly_fit(fitdata[i],logfreq,num_planes,2.0,C,1,&chisq);
		printf("%d index=%lf\n",i,C[1]);
	}
	
	return EXIT_SUCCESS;
}


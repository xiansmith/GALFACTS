#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsio.h"
#include "jsd/jsd_futil.h"

int main(int argc, char *argv[]) 
{

	int i, j, k, n;
	char *Ifilename, *Vfilename,*duchampfile;
	FILE *Ifile,*Vfile,*xyfile,*outfile;
	header_param_list Ihpar, Vhpar;
	float *Idata, *Vdata;
	int plane_size;

	outfile = fopen("VvsIdata.txt","w");
	fprintf(outfile,"#RA DEC I V absV\n");
	if (argc < 3) {
		return EXIT_FAILURE;
	}
	Ifilename = argv[1];
	Vfilename = argv[2];

	Ifile = fopen(Ifilename,"r");
	Vfile = fopen(Vfilename,"r");
	//xyfile = fopen("XY.txt","r");
	xyfile = fopen("N1cc.dat","r");
        readfits_header (Ifile, &Ihpar);
        readfits_header (Vfile, &Vhpar);

	//int count = jsd_line_count(xyfile)-1;
	int count = jsd_line_count(xyfile);
	printf("number of sources :%d\n",count);

	plane_size = Ihpar.naxis[0] * Ihpar.naxis[1];
	printf("plane_size: %i xaxis: %i yaxis: %i\n", plane_size,Ihpar.naxis[0],Ihpar.naxis[1]);
	Idata = malloc (sizeof(float) * plane_size);
	Vdata = malloc (sizeof(float) * plane_size);

	int nx,ny;
	float crpx,crpy,delx,dely,crvx,crvy;

        nx = Ihpar.naxis[0];
        ny = Ihpar.naxis[1];
        crpx = Ihpar.crpix[0];
        crpy = Ihpar.crpix[1];
        delx = Ihpar.cdelt[0];
        dely = Ihpar.cdelt[1];
        crvx = Ihpar.crval[0];
        crvy = Ihpar.crval[1];

	readfits_plane (Ifile, Idata, &Ihpar);
	readfits_plane (Vfile, Vdata, &Vhpar);
	
	float RA,DEC,stokesI;
	char string[20];
	char header[81];
	fgets(header,80,xyfile);
	for(i=0;i<count;i++)
	{
		int arrx,arry;
		fscanf(xyfile,"%f %f %f",&RA,&DEC,&stokesI);
		
		arrx = (int)(crpx+(RA-crvx)/delx)-1;
		arry = (int)(crpy+(DEC-crvy)/dely)-1;

/*		if(X-(int)(X)>=0.5)
			xx = (int)(X)+1;
		else
			xx = (int)X;
		if(Y-(int)(Y)>=0.5)
			yy = (int)(Y)+1;
		else
			yy = (int)Y;
*//*		if(xx==Ihpar.naxis[0])
			xx--;
		if(yy==Ihpar.naxis[1])
			yy--;
*/		
//		printf("%f %d %f %d\n",X,xx,Y,yy);
//		printf("%d\n",(int)(Y)*Ihpar.naxis[0]+(int)(X));
//		printf("X %f Y %f I %f V %f absV %f\n",X+609,Y+220,Idata[(yy+220)*Ihpar.naxis[0]+xx+609],Vdata[(yy+220)*Ihpar.naxis[0]+xx+609],fabs(Vdata[(yy+220)*Ihpar.naxis[0]+xx+609]));
		//fprintf(outfile,"%f %f %f %f %f\n",X+830+1,Y+0+1,Idata[(yy+1)*Ihpar.naxis[0]+xx+830+1],Vdata[(yy+1)*Ihpar.naxis[0]+xx+830+1],fabs(Vdata[(yy+1)*Ihpar.naxis[0]+xx+830+1]));
		fprintf(outfile,"%2.6f %2.6f %2.6f %2.6f %2.6f %2.6f %d %d\n",fabs(Vdata[arry*nx+arrx])/Idata[arry*nx+arrx],RA,DEC,stokesI/1000.0,Idata[arry*nx+arrx],fabs(Vdata[arry*nx+arrx]),arrx,arry);
	}
	fclose(xyfile);
	fclose(Ifile);
	fclose(Vfile);
	fclose(outfile);
	return EXIT_SUCCESS;
}


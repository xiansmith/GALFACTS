#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "programs/fitsio.h"
#include "common.h"
#include "jsd/jsd_futil.h"

int main(int argc, char *argv[]) 
{

	int f,N1,N2,n1,n2,i,i1,j,k,numfiles;
	char outcube_name[40];
	char * incube_name;
	char * listfilename;
	FILE * incube_file;
	FILE * outcube_file;
	header_param_list incube_hpar;
	//header_param_list outcube_hpar;
	int inplane_len;
	int outplane_len;
	int outdata_len,xcount,ycount;
	float * inplane_data;
	float * outplane_data;
	float * outcube_data;
	char ** incubes;
	float start;

	if (argc != 4) 
	{
		printf("Usage: %s <filelist> <x> <y>\n", argv[0]);
		return EXIT_FAILURE;
	}
	else
	{
		listfilename = argv[1];
		xcount=atoi(argv[2]);
		ycount=atoi(argv[3]);
	}

        FILE * listfile;
        listfile = fopen(listfilename,"r");
        if(listfile != NULL)
        {
                printf("Found file %s\n",listfilename);
                numfiles = jsd_line_count(listfile);
                incubes = (char **) malloc(sizeof(char*) * numfiles);
                char tempstr[51];
                for(i=0;i<numfiles;i++)
                {
                        fscanf(listfile,"%s",tempstr);
                        incubes[i] = (char *) malloc(sizeof(char) * strlen(tempstr));
                        strcpy(incubes[i],tempstr);
                }
                fclose(listfile);
        }
	else
		exit(1);



	for(f=0;f<numfiles;f++)
	{
	incube_file = fopen(incubes[f], "r");
	if(incube_file == NULL)
		printf("Failed to open file\n");
	readfits_header(incube_file, &incube_hpar);
	printf("Reading %s\n",incubes[f]);
	//allocate memory for cube
	//outdata_len = incube_hpar.naxis[0]/10 * incube_hpar.naxis[1]/2 * 900;
	outdata_len = incube_hpar.naxis[0]/10 * incube_hpar.naxis[1]/2 * 1700;
	if(f ==0)
	{
		outcube_data = (float*) calloc(outdata_len, sizeof (float));
		outplane_len = incube_hpar.naxis[0]/10 * incube_hpar.naxis[1]/2;
		start =  incube_hpar.crval[2];
	}
	inplane_len = incube_hpar.naxis[0] * incube_hpar.naxis[1];
	inplane_data = (float*) calloc(inplane_len, sizeof (float));
	//outplane_data = (float*) calloc(outplane_len, sizeof (float));
	//printf("outplane %d outcube %d inplane %d offset %d\n",outplane_len,outdata_len,inplane_len,f*outplane_len*incube_hpar.naxis[2]);
	N1 = incube_hpar.naxis[0];
	n1 = incube_hpar.naxis[0]/10;
	N2 = incube_hpar.naxis[1];
	n2 = incube_hpar.naxis[1]/2;
	//iterate over the cube
	//for (i=0,i1=0; i<incube_hpar.naxis[2]; i=i+2,i1++) 
	for (i=0; i<incube_hpar.naxis[2];i++) 
	{
		//read a plane
		//printf("Plane %d\n",i);
		readfits_plane(incube_file, inplane_data, &incube_hpar);
		//readfits_plane(incube_file, inplane_data, &incube_hpar);
		//include in average only if stated in frames
		for(j=0;j<n2;j++)
		{
			//printf("%i \n", i);
			for (k=0; k<n1; k++) 
			{
				//add to average and increment count if its not blank
					outcube_data[j*n1+k+i*outplane_len+f*outplane_len*incube_hpar.naxis[2]] = inplane_data[j*incube_hpar.naxis[0]+k+n1*xcount+ycount*incube_hpar.naxis[0]*n2];
					//outcube_data[j*n1+k+i1*outplane_len+f*outplane_len*incube_hpar.naxis[2]/2] = inplane_data[j*incube_hpar.naxis[0]+k+n1*xcount+ycount*incube_hpar.naxis[0]*n2];
			}
		}
	}
	free(inplane_data);
	fclose(incube_file);
	}
	
	float Xcen,Ycen;
	Xcen = incube_hpar.crval[0];
	Ycen = incube_hpar.crval[1];
	//write the average fits
        incube_hpar.crval[0] = Xcen - N1*(incube_hpar.cdelt[0])/2+(2*xcount+1)*(incube_hpar.cdelt[0])*n1/2;         /* degrees */
        incube_hpar.crval[1] = Ycen - N2*(incube_hpar.cdelt[1])/2+(2*ycount+1)*(incube_hpar.cdelt[1])*n2/2;         /* degrees */
        incube_hpar.crpix[0] = 0.5 + n1 / 2.0;        /* image center in pixels */
        incube_hpar.crpix[1] = 0.5 + n2 / 2.0;
	//incube_hpar.naxis[2] = 900;
	incube_hpar.naxis[2] = 1700;
	incube_hpar.naxis[0] /= 10;
	incube_hpar.naxis[1] /= 2;
	incube_hpar.crval[2] = start;
	//incube_hpar.cdelt[2] *= 2;
	//strcat(cube_hpar.object, " Average");
	//printf("Writing output cube fstart %f\n",start);
	sprintf(outcube_name,"D%d_Qcube.fits",xcount+ycount*10+1);
	printf("Writing %s\n",outcube_name);
	writefits_cube(outcube_name, outcube_data, &incube_hpar,0);
 	FILE * annfile;
	char annfilename[41];
	sprintf(annfilename,"D%d.ann",xcount+ycount*10+1);
	annfile = fopen(annfilename,"w");
	fprintf(annfile,"BOX %f %f %f %f",Xcen - (N1*(incube_hpar.cdelt[0])/2)+xcount*n1*incube_hpar.cdelt[0],Ycen - (N2*(incube_hpar.cdelt[1])/2)+ycount*n2*incube_hpar.cdelt[1],Xcen-(N1*(incube_hpar.cdelt[0])/2)+(xcount+1)*n1*incube_hpar.cdelt[0],Ycen-(N2*(incube_hpar.cdelt[1])/2)+(ycount+1)*n2*incube_hpar.cdelt[1]);
	fclose(annfile);
	return EXIT_SUCCESS;
}

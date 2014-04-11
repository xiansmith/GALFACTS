#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "programs/fitsio.h"

int main(int argc, char * argv[])
{
	char *suf;
	int beam,lowchan,highchan;

	FILE *ifile,*vfile;

	char filename[100];

	suf = argv[1];
        beam = atoi(argv[2]);
        lowchan = atoi(argv[3]);
        //lowchan = 0;
        highchan = atoi(argv[4]);
        //highchan = 1;

	sprintf(filename,"%s_BEAM%d_%04i_%04i_Icube.fits",suf,beam,lowchan,highchan-1);
	//sprintf(filename,"%s_BEAM%d_Icube.fits",suf,beam);
	ifile = fopen(filename,"r");
	sprintf(filename,"%s_BEAM%d_%04i_%04i_Vcube.fits",suf,beam,lowchan,highchan-1);
	//sprintf(filename,"%s_BEAM%d_Vcube.fits",suf,beam);
	vfile = fopen(filename,"r");

        header_param_list ihpar,vhpar;
	int dtlen,idxmx,idxV;
	float *plnI,*plnV;

	readfits_header(ifile, &ihpar);
	readfits_header(vfile, &vhpar);
	dtlen = ihpar.naxis[0] * ihpar.naxis[1];

	plnI = (float*) calloc(dtlen, sizeof (float));
        plnV = (float*) calloc(dtlen, sizeof (float));

	int i,k;
	float imxI = -100.0;
	float imxV = +100.0;
	readfits_plane(ifile, plnI, &ihpar);
	readfits_plane(vfile, plnV, &vhpar);


        for (k=0; k<dtlen; k++)
        {
        	if (!IS_BLANK_PIXEL(plnI[k]) && isfinite(plnI[k]))
                { 
        	        if(imxI < plnI[k])
                        {
                	        imxI = plnI[k];
                                idxmx = k;
                        }
                }
        }

        for (k=0; k<dtlen; k++)
        {
        	if (!IS_BLANK_PIXEL(plnV[k]) && isfinite(plnV[k]))
                { 
        	        if(imxV > plnV[k])
                        {
                	        imxV = plnV[k];
                                idxV = k;
                        }
                }
        }
	printf("chan %04i i %2.6f v %2.6f  v/i %2.6f idx %d idxV %d\n",lowchan,plnI[idxmx],plnV[idxV],plnV[idxV]/plnI[idxmx],idxmx,idxV);

	for(i = lowchan+1;i<highchan;i++)
	{
		readfits_plane(ifile, plnI, &ihpar);
		readfits_plane(vfile, plnV, &vhpar);
		printf("chan %04i i %2.6f v %2.6f  v/i %2.6f idx %d idxV %d\n",i,plnI[idxmx],plnV[idxV],plnV[idxV]/plnI[idxmx],idxmx,idxV);
	}

	fclose(ifile);
	fclose(vfile);
}	

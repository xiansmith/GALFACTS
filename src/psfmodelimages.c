#include "jsd/jsd_futil.h"
#include "generic_grid.h"
#include "map.h"
#include "programs/fitsio.h"
#include <math.h>

int main(int argc, char * argv[])
{
	char filename[64+1];
	FILE *psffile;
	char header[80+1];
	char title[64+1];
	dim2data *read;
	float *phase_ap, *err_im, *mod_ap;
	MapMetaData md;
	float xmin,xmax,ymin,ymax;
	int linecount,j,beamno;

	if (argc != 2) {
		printf("Incorrect number of arguments\n");
		return EXIT_FAILURE;
	} 
	else
	{
		beamno = atoi(argv[1]);
	}

	md.fcen = 0;
	md.lowchan = 0;
	md.highchan = 0;
	md.gridtype = 1;
	sprintf(title,"Beam %d modelled aperture illumination.",beamno);
	md.title = title;
	md.fstart = 0;
	
	sprintf(filename,"beam%d_psfmodel_raw.dat",beamno);
//	sprintf(filename,"zernike_3_1.dat");
	psffile = fopen(filename,"r");
        fgets(header, 80, psffile);
	linecount = jsd_line_count(psffile);
	read = (dim2data *)malloc(sizeof (dim2data) * linecount);
	phase_ap = (float *)malloc(sizeof (float) * linecount);
	err_im = (float *)malloc(sizeof (float) * linecount);
	mod_ap = (float *)malloc(sizeof (float) * linecount);

	xmin = INFINITY;
	ymin = INFINITY;
	xmax = 0;
	ymax = 0;

	for(j = 0;j < linecount;j++)
	{
//		fscanf(psffile,"%f %f %f",&read[j].x,&read[j].y,&mod_ap[j]);
		fscanf(psffile,"%f %f %f %f",&read[j].x,&read[j].y,&mod_ap[j],&phase_ap[j]);
		if(xmin > read[j].x)
			xmin = read[j].x;
		if(ymin > read[j].y)
			ymin = read[j].y;
		if(xmax < read[j].x)
			xmax = read[j].x;
		if(ymax < read[j].y)
			ymax = read[j].y;
		read[j].z = mod_ap[j];
	}
	fclose(psffile);

	printf("xmax:%f xmin:%f ymax:%f ymin:%f \n",xmax,xmin,ymax,ymin);

	md.ramin = xmin; 
	md.ramax = xmax;
	md.decmin = ymin;
	md.decmax = ymax;
	md.RAcen = (md.ramax + md.ramin)/2.0;
	md.DECcen = (md.decmax + md.decmin)/2.0;
	md.RArange =  md.ramax - md.ramin;
	md.DECrange = md.decmax - md.decmin;
	md.cellsize = md.RArange/100;
	md.n1 = (int)(md.RArange/md.cellsize) + 1;
	md.n2 = (int)(md.DECrange/md.cellsize) + 1;
	md.n3 = md.highchan - md.lowchan;
	init_psf_lookup_table(65537, xmax, md.RArange/11);
	md.patch =  10;//(int)(md.RArange/md.cellsize*10);

	sprintf(filename,"beam%d_ap_ilum_raw.fits",beamno);
//	sprintf(filename,"zernike_3_1.fits");
	generic_grid_data(md,read,linecount,filename);

	for(j = 0;j < linecount;j++)
	{
		read[j].z = phase_ap[j];		
	}

	sprintf(title,"Beam %d modelled aperture phase.",beamno);
	md.title = title;
	sprintf(filename,"beam%d_ap_phase_raw.fits",beamno);
	generic_grid_data(md,read,linecount,filename);

	sprintf(filename,"beam%d_psf_err.dat_raw",beamno);
	psffile = fopen(filename,"r");
        fgets(header, 80, psffile);
	linecount = jsd_line_count(psffile);

	xmin = INFINITY;
	ymin = INFINITY;
	xmax = 0;
	ymax = 0;

	free_psf_lookup_table();

	for(j = 0;j < linecount;j++)
	{
		fscanf(psffile,"%f %f %f %f %f",&read[j].x,&read[j].y,&mod_ap[j],&err_im[j],&phase_ap[j]);
		read[j].x = read[j].x*57.3;
		read[j].y = read[j].y*57.3;
		if(xmin > read[j].x)
			xmin = read[j].x;
		if(ymin > read[j].y)
			ymin = read[j].y;
		if(xmax < read[j].x)
			xmax = read[j].x;
		if(ymax < read[j].y)
			ymax = read[j].y;
		read[j].z = err_im[j];
	}
	fclose(psffile);
	md.ramin = xmin; 
	md.ramax = xmax;
	md.decmin = ymin;
	md.decmax = ymax;
	md.RAcen = (md.ramax + md.ramin)/2.0;
	md.DECcen = (md.decmax + md.decmin)/2.0;
	md.RArange =  md.ramax - md.ramin;
	md.DECrange = md.decmax - md.decmin;
	md.cellsize = md.RArange/100;
	md.n1 = (int)(md.RArange/md.cellsize) + 1;
	md.n2 = (int)(md.DECrange/md.cellsize) + 1;
	md.n3 = md.highchan - md.lowchan;
	md.patch = 10;
	init_psf_lookup_table(65537, xmax, md.RArange/11);

	sprintf(title,"Beam %d modelled psf error.",beamno);
	md.title = title;
	sprintf(filename,"beam%d_psfmodel_err_raw.fits",beamno);
	printf("xmax:%f xmin:%f ymax:%f ymin:%f \n",xmax,xmin,ymax,ymin);
	generic_grid_data(md,read,linecount,filename);

	for(j = 0;j < linecount;j++)
	{
		read[j].z = mod_ap[j];
	}
	sprintf(title,"Beam %d modelled psf.",beamno);
	md.title = title;
	sprintf(filename,"beam%d_psf_model_raw.fits",beamno);
	printf("xmax:%f xmin:%f ymax:%f ymin:%f \n",xmax,xmin,ymax,ymin);
	generic_grid_data(md,read,linecount,filename);

	for(j = 0;j < linecount;j++)
	{
		read[j].z = phase_ap[j];
	}
	sprintf(title,"Beam %d measured psf.",beamno);
	md.title = title;
	sprintf(filename,"beam%d_psf_orig_raw.fits",beamno);
	printf("xmax:%f xmin:%f ymax:%f ymin:%f \n",xmax,xmin,ymax,ymin);
	generic_grid_data(md,read,linecount,filename);

	printf("Done.\n");
	return EXIT_SUCCESS;
}

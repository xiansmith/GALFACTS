#include "fluxdata.h"
#include "jsd/jsd_futil.h"
#include "psfmodel.h"
#include <gsl/gsl_multifit.h>
//#include <math.h>
#undef I
#define MODULUS 1
#define PHASE 0

int main(int argc, char * argv[])
{
	int linecount,j,k,counter,apcount,order,num_iter,beamno,flag,apcountsave;
	FluxRecord * rec;
	char filename[64+1],header[80+1];
	FILE * beammodelfile,*errfile,*psffile;
	dftdata *imgdata, *apdata;
	double *mod_img,*phase_img,*mod_ap,*phase_ap,*phasen_ap,*modn_ap, *pct_err, *psf_orig;
	double mse,mse_prev,maxerr,chisq,global_min;
	float *x,*y,*znk,*weight,apradius,apfield,epsilon,change,begin_r,end_r,min_r,step_r;
	gsl_matrix *x_gsl, *cov;
	gsl_vector *y_gsl, *c;
	
	if (argc != 6) {
		printf("Incorrect number of arguments\n");
		return EXIT_FAILURE;
	} 
	else
	{
		beamno = atoi(argv[1]);
		order = atoi(argv[2]);
		begin_r = atof(argv[3]);
		end_r = atof(argv[4]);
		step_r = atoi(argv[5]);
//		sigma = atoi(argv[6]);
	}
	
//	sprintf(filename, "beam%d_model/beam%d_model050.dat",beamno,beamno); ///for the time being
	sprintf(filename, "psftestcase_%d.dat",beamno); ///for the time being
	beammodelfile = fopen(filename, "r");
        fgets(header, 80, beammodelfile);
	linecount = jsd_line_count(beammodelfile);
	printf("Beginning PSF modelling.\n");
	printf("%d datapoints in the beam model.\n",linecount);


	rec = (FluxRecord *)malloc(sizeof (FluxRecord) * linecount);
	imgdata = (dftdata *)malloc(sizeof (dftdata) * linecount);
	apdata = (dftdata *)malloc(sizeof (dftdata) * linecount);
	mod_img = (double *)malloc(sizeof (double) * linecount);
	phase_img = (double *)malloc(sizeof (double) * linecount);
	mod_ap = (double *)malloc(sizeof (double) * linecount);
	phase_ap = (double *)malloc(sizeof (double) * linecount);
	phasen_ap = (double *)malloc(sizeof (double) * linecount);
	modn_ap = (double *)malloc(sizeof (double) * linecount);
	pct_err = (double *)malloc(sizeof (double) * linecount);
	psf_orig = (double *)malloc(sizeof (double) * linecount);

	apcount = 0;

	printf("Reading beam model.\n");
	for(j = 0;j < linecount;j++)
	{
		fscanf(beammodelfile,"%f %f %lf %lf",&rec[j].RA,&rec[j].DEC,&rec[j].stokes.I,&rec[j].stokes.Q);
//		fscanf(beammodelfile,"%f %f %f %lf %lf %lf %lf",&rec[j].RA,&rec[j].DEC,&rec[j].AST,&rec[j].stokes.I,&rec[j].stokes.Q,\
		&rec[j].stokes.U, &rec[j].stokes.V);

		if(rec[j].stokes.I < 0) //removing negatives before taking square root not sure how to handle this in a better way
		{
			mod_img[j] = 0;
			psf_orig[j] = 0.0;			
		}
		else
		{
			mod_img[j] = sqrt(rec[j].stokes.I);
//			mod_img[j] = rec[j].stokes.I;  
			psf_orig[j] = rec[j].stokes.I;
		}

		//setting up img coord in theta radians
		imgdata[j].x = rec[j].RA/57.3; //coordinates in radians
		imgdata[j].y = rec[j].DEC/57.3;	
//		imgdata[j].z = mod_img[j];
		//setting up aperture coord in lambdas	
/*		apdata[j].x = (rec[j].RA/0.15)*apfield; 
		apdata[j].y = (rec[j].DEC/0.15)*apfield;

		//setting up modulus of aperture (1 within the aperture 0 outside), setting up random phase for gerchberg saxton algorithm
		if((apdata[j].x*apdata[j].x+apdata[j].y*apdata[j].y) <= apradius*apradius)
		{
			mod_ap[j] = 1.0;
			modn_ap[j] = 1.0;
			phase_ap[j] = 0.0;
			phasen_ap[j] = 0.0;
			apcount++;
		}
		else
		{
			mod_ap[j] = 0.0;
			modn_ap[j] = 0.0;
			phase_ap[j] = 0.0;
			phasen_ap[j] = 0.0;
		}

		//calculating modulus in aperture plane from modulus and phase generated above
		apdata[j].z = mod_ap[j]*cexp(_Complex_I*phase_ap[j]);
*/	}

	fclose(beammodelfile);
	free(rec);

	global_min = INFINITY;
	for(apradius = begin_r;apradius <= end_r;apradius = apradius+step_r) //start for loop, searching thru varius apradius values
	{
		apfield = apradius;
		apcount = 0;
//		printf("--------------------------------------------\n");
		printf("Modelling using aperture radius value of %f\n",apradius);
//		printf("--------------------------------------------\n");

		for(j = 0;j < linecount;j++)
		{

			//setting up aperture coord in lambdas	
			apdata[j].x = (imgdata[j].x*57.3/0.15)*apfield; 
			apdata[j].y = (imgdata[j].y*57.3/0.15)*apfield;

			//setting up modulus of aperture (1 within the aperture 0 outside), setting up random phase for gerchberg saxton algorithm
			if((apdata[j].x*apdata[j].x+apdata[j].y*apdata[j].y) <= apradius*apradius)
			{
//				mod_ap[j] = 0.1;
				mod_ap[j] = (RAND_MAX-rand())/(RAND_MAX)+0.01;
//				modn_ap[j] = 1.0;
//				phase_ap[j] = 0.0;
//				phase_ap[j] = ((rand()-RAND_MAX/2)*3.14)/(RAND_MAX/2);
				phasen_ap[j] = 0.0;
				apcount++;
			}
			else
			{
				mod_ap[j] = 0.0;
//				modn_ap[j] = 0.0;
				phase_ap[j] = 0.0;
//				phasen_ap[j] = 0.0;
			}

			//calculating modulus in aperture plane from modulus and phase generated above
			apdata[j].z = mod_ap[j]*cexp(_Complex_I*phase_ap[j]);
		}
//	idft2d(imgdata,linecount,apdata);	
//	exit(1);

//	errfile = fopen("err.dat","w");	
//	fprintf(errfile,"#Img Ap\n");

	//Gerchberg Saxton main loop
//	printf("Starting modified Gerchberg-Saxton algorithm.\n");

	epsilon = 0.0001;
	counter = 0;
//	flag = PHASE;
	flag = MODULUS;	

	do
	{
		mse_prev = 0.0;
		do
		{
//			printf("Fourier transforming the aperture plane data.\n");
			dft2d(apdata,linecount,imgdata);

			mse = 0.0;
			maxerr = 0.0;
			
			//retreiving phase from the fourier transform and setting up for inverse transform
			for(j = 0;j < linecount;j++)
			{
				phase_img[j] = carg(imgdata[j].z);
//				phase_img[j] = atan2(creal(imgdata[j].z),cimag(imgdata[j].z));
				mse += fabs(psf_orig[j]-cabs(imgdata[j].z)*cabs(imgdata[j].z));//*(psf_orig[j]-cabs(imgdata[j].z)*cabs(imgdata[j].z));
				if(maxerr < fabs(psf_orig[j] - cabs(imgdata[j].z)*cabs(imgdata[j].z)))
					maxerr = fabs(psf_orig[j] - cabs(imgdata[j].z)*cabs(imgdata[j].z));
				imgdata[j].z = mod_img[j]*cexp(_Complex_I*phase_img[j]);
			}
			printf("Iteration %d, MSE:%1.9f MaxErr:%1.9f\n",counter+1,mse/linecount,maxerr);
			change = mse - mse_prev;
			mse_prev = mse;
			if(maxerr < global_min)
			{
				global_min = maxerr;
				min_r = apradius;
				apcountsave = apcount;
			}

//			fprintf(errfile,"%lf ",mse/denom1);

//			printf("Inverse Fourier transforming the intensity data.\n");
//			if(global_min > 0.119) //hack
			idft2d(imgdata,linecount,apdata);

			//retreiving phase from the fourier transform
			for(j = 0;j < linecount;j++)
			{
				if((apdata[j].x*apdata[j].x+apdata[j].y*apdata[j].y) <= apradius*apradius) //todo remove hard coding
				{
					if(flag == PHASE)
						phasen_ap[j] = carg(apdata[j].z);
//						phasen_ap[j] = atan2(creal(apdata[j].z),cimag(apdata[j].z));
					else
						modn_ap[j] = cabs(apdata[j].z);
				}
				else	
				{
					if(flag == PHASE)
						phasen_ap[j] = 0;
					else	
						modn_ap[j] = 0;
				}
			}	

			//setup new aperture plane
			for(j = 0;j < linecount;j++)
			{
				if(flag == PHASE)
					apdata[j].z = mod_ap[j]*cexp(_Complex_I*phasen_ap[j]);
				else
					apdata[j].z = modn_ap[j]*cexp(_Complex_I*phase_ap[j]);
			}
			counter++;
//		}while(counter < 20);
//		}while(global_min > 0.119 && fabs(change) > epsilon);
		}while(fabs(change) > epsilon);

		printf("Iteration %d, MSE:%1.9f MaxErr:%1.9f\n",counter+1,mse/linecount,maxerr);

		for(j = 0;j < linecount;j++)
		{
			if(flag == PHASE)
				phase_ap[j] = phasen_ap[j];
			else
				mod_ap[j] = modn_ap[j];
		}
		counter = 0;	
		epsilon/=10;
		if (flag==PHASE)
		{
			printf("Switching to optimizing modulus.\n");
			flag = MODULUS;
//			flag = PHASE;
		}
		else
		{
			printf("Switching to optimizing phase.\n");
//			flag = PHASE;
			flag = MODULUS;
		}
		
//	}while(fabs(change) > 0.000000001);
//	}while(epsilon > 0.1);
//	}while(global_min > 0.119);
	}while(epsilon > 0.0000000001); //number of zeros control whether the loop will end optimizing phase/modulus. It should always end with ??? phase.

//	printf("MSE:%1.9f MaxErr:%1.9f\n",mse/linecount,maxerr);
	
	}//for loop ends here

	printf("Global minimum of %lf at achieved at aperture radius value %f\n",global_min,min_r);

//	sprintf(filename,"beam%d_psfmodel.dat",beamno);
	sprintf(filename, "psftestcase_%d_model.dat",beamno); ///for the time being
	psffile = fopen(filename,"w");
	fprintf(psffile,"#x y mod phase\n");
	for(j = 0;j < linecount;j++)
	{
		fprintf(psffile,"%f %f %lf %lf\n",apdata[j].x,apdata[j].y,mod_ap[j],phase_ap[j]);
	}
	fclose(psffile);

//	sprintf(filename,"beam%d_psf.dat",beamno);
	sprintf(filename,"psfretrieved_%d.dat",beamno);
	psffile = fopen(filename,"w");
//	fprintf(psffile,"#Parameters for the run");
	fprintf(psffile,"#Aperture radius %f ",apradius);
	fprintf(psffile,"#Max error %f\n ",global_min);
//	fprintf(psffile,"#Aperture field %f\n",apfield);
	fprintf(psffile,"#x y rootI I err\n");
	dft2d(apdata,linecount,imgdata);
	for(j = 0;j < linecount;j++)
	{
	 fprintf(psffile,"%f %f %lf %lf %lf\n",imgdata[j].x,imgdata[j].y,cabs(imgdata[j].z),pow(cabs(imgdata[j].z),2),pow(cabs(imgdata[j].z),2)-psf_orig[j]);
	}
	fclose(psffile);
		
/*	free(mod_img);
	free(phase_img);
//	free(mod_ap);
//	free(phase_ap);
	free(modn_ap);
//	free(imgdata);
	free(pct_err);
	
*/	//setting up for least squares fitting
//	apcount = apcountsave;

	apradius = begin_r;
	apcount = 0;
//	sprintf(filename,"beam%d_psfmodel.dat",beamno);
	sprintf(filename, "psftestcase_%d_model.dat",beamno); ///for the time being
	psffile = fopen(filename,"r");
        fgets(header, 80, psffile);
	for(j = 0;j < linecount;j++)
	{
		fscanf(psffile,"%f %f %lf %lf",&apdata[j].x,&apdata[j].y,&mod_ap[j],&phase_ap[j]);
		if((apdata[j].x*apdata[j].x+apdata[j].y*apdata[j].y) <= apradius*apradius)
			apcount++;
		apdata[j].z = mod_ap[j]*cexp(_Complex_I*phase_ap[j]);
	}
	fclose(psffile);

	dft2d(apdata,linecount,imgdata);

	for(j = 0;j < linecount;j++)
	{
		if(maxerr < fabs(psf_orig[j] - cabs(imgdata[j].z)*cabs(imgdata[j].z)))
			maxerr = fabs(psf_orig[j] - cabs(imgdata[j].z)*cabs(imgdata[j].z));
	}
	printf("MaxErr in orig model is:%1.9f\n",maxerr);


	znk = (float*)malloc(sizeof(float)*apcount);
	x = (float*)malloc(sizeof(float)*apcount);
	y = (float*)malloc(sizeof(float)*apcount);
	weight = (float*)malloc(sizeof(float)*order);
	x_gsl = gsl_matrix_alloc(apcount,order);
	y_gsl = gsl_vector_alloc(apcount);
	c = gsl_vector_alloc(order);
	cov = gsl_matrix_alloc(order,order);
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (apcount,order);

//	printf("Setting up least square fitting matrices.\n");
	counter = 0;
//	apfield = apradius = min_r;
	
	for(j = 0;j < linecount;j++)
	{
		if((apdata[j].x*apdata[j].x+apdata[j].y*apdata[j].y) <= apradius*apradius)
		{
//			printf("\ncounter:%d",counter);
			x[counter]=apdata[j].x;
			y[counter]=apdata[j].y;
			gsl_vector_set(y_gsl,counter,phase_ap[j]);
			counter++;
		}
	}

	int n,m;
	n=0;
	m=0;
	for(j = 0;j < order;j++)
	{
//		printf("Calculating Zernike polynomial n:%d m:%d\n",n,m);
		zernike(n,m,x,y,apcount,znk,apradius);
		for(k = 0;k < apcount;k++)
		{
			gsl_matrix_set(x_gsl,k,j,znk[k]);
		}
		if(m==n)
		{
			n++;
			m = -n;
		}
		else
			m+=2;
	}
	
	//perform the fit
	printf("Performing the fit.\n");
	gsl_multifit_linear(x_gsl,y_gsl,c,cov,&chisq,work);

	//store the results in the output array
	n=0;
	m=0;
	for (j=0; j < order; j++)
	{
		weight[j]=gsl_vector_get(c,j);
		printf("n:%d m:%d calculated weight %f\n",n,m,weight[j]);
		if(m==n)
		{
			n++;
			m = -n;
		}
		else
			m+=2;
	}
///////////////////////////////////////////////

//	sprintf(filename,"beam%d_psfmodel_params.dat",beamno);
/*	sprintf(filename, "psftestcase_%d_model_params.dat",beamno); ///for the time being
	psffile = fopen(filename,"w");
	fprintf(psffile,"# weight\n");
	for(j = 0;j < order;j++)
	{
		fprintf(psffile,"%d %f\n",j,weight[j]);
	}
	fclose(psffile);

	n=0;
	m=0;

	for(j = 0;j < linecount;j++)
	{
		phase_ap[j] = 0;
	}

	for(j = 0;j < order;j++)
	{
		for(k=0;k<linecount;k++)
		{
			if((apdata[k].x*apdata[k].x+apdata[k].y*apdata[k].y) <= apradius*apradius)
			{
				x[0] = apdata[k].x;
				y[0] = apdata[k].y;
				zernike(n,m,x,y,1,znk,apradius);
				phase_ap[k] += weight[j]*znk[0];
			}
			else
				phase_ap[k] = 0;
		}

		if(m==n)
		{
			n++;
			m = -n;
		}
		else
			m+=2;
	}
	
	for(j = 0;j < linecount;j++)
	{
		apdata[j].z = mod_ap[j]*cexp(_Complex_I*phase_ap[j]);
	}

	dft2d(apdata,linecount,imgdata);

	maxerr = 0.0;
	//retreiving phase from the fourier transform and setting up for inverse transform

//	sprintf(filename,"beam%d_psf_err.dat",beamno);
	sprintf(filename, "psftestcase_%d_model_err.dat",beamno); ///for the time being
	psffile = fopen(filename,"w");
	fprintf(psffile,"#x y I(model) err orig\n");
	mse = 0.0;
	for(j = 0;j < linecount;j++)
	{
		fprintf(psffile,"%f %f %lf %lf %lf\n",imgdata[j].x,imgdata[j].y,cabs(imgdata[j].z)*cabs(imgdata[j].z),(psf_orig[j] - cabs(imgdata[j].z)*cabs(imgdata[j].z)),psf_orig[j]);
		mse += fabs(psf_orig[j]-cabs(imgdata[j].z)*cabs(imgdata[j].z));
		if(maxerr < fabs(psf_orig[j] - cabs(imgdata[j].z)*cabs(imgdata[j].z)))
			maxerr = fabs(psf_orig[j] - cabs(imgdata[j].z)*cabs(imgdata[j].z));
	}
	printf("MaxErr in model with order %d is:%1.9f, avg error is %1.9f\n",order,maxerr,mse/linecount);
	fclose(psffile);
*/////
	free(apdata);
//	gsl_matrix_free(x_gsl);
//	gsl_vector_free(y_gsl);
//	gsl_vector_free(c);
//	gsl_matrix_free(cov);
//	gsl_multifit_linear_free (work);
	free(x);
	free(y);
	free(znk);
	
	return EXIT_SUCCESS;
}

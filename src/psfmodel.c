#include <stdlib.h>
#include <stdio.h>
#include "psfmodel.h"
#include <math.h>

void dft2d(dftdata input[], int numpoints, dftdata output[])
{
	FILE *infile,*outfile; //for the time being
	int i,j;
	float max = 0.0; 

	infile = fopen("DFT_indata.dat","w");
	outfile = fopen("DFT_outdata.dat","w");
	fprintf(infile,"#X Y Z\n");
	fprintf(outfile,"#X Y Z\n");
	for(i = 0;i < numpoints;i++)
	{
		output[i].z = 0;
		for(j = 0;j < numpoints;j++)
		{
				output[i].z += input[j].z*cexp(-2*M_PI*_Complex_I*(input[j].x*output[i].x+input[j].y*output[i].y));
		}

//		output[i].z/=(sqrt(2*PI)*500);
		
		//for normalizing
		if(cabs(output[i].z) > max)
			max = cabs(output[i].z);

	}
	
	for(i = 0;i < numpoints;i++)
	{
		output[i].z = output[i].z/max;
		fprintf(infile,"%2.8f %2.8f %2.8lf+i%2.8lf %2.8lf %2.8lf\n",input[i].x,input[i].y,creal(input[i].z),cimag(input[i].z),cabs(input[i].z),carg(input[i].z));
//		fprintf(infile,"%2.8f %2.8f %2.8lf+i%2.8lf %2.8lf %2.8lf\n",input[i].x,input[i].y,creal(input[i].z),cimag(input[i].z),cabs(input[i].z),\
		atan2(creal(input[i].z),cimag(input[i].z)));
	       fprintf(outfile,"%2.8f %2.8f %2.8lf+i%2.8lf %2.8lf %2.8lf\n",output[i].x,output[i].y,creal(output[i].z),cimag(output[i].z),cabs(output[i].z),carg(output[i].z));
//	       fprintf(outfile,"%2.8f %2.8f %2.8lf+i%2.8lf %2.8lf %2.8lf\n",output[i].x,output[i].y,creal(output[i].z),cimag(output[i].z),cabs(output[i].z),\
		atan2(creal(output[i].z),cimag(output[i].z)));
	}

	fclose(infile);
	fclose(outfile);
}

void idft2d(dftdata input[], int numpoints, dftdata output[])
{
	FILE *infile,*outfile; //for the time being
	int i,j;
	float max = 0.0; 

	infile = fopen("iDFT_indata.dat","w");
	outfile = fopen("iDFT_outdata.dat","w");
	fprintf(infile,"#X Y Z\n");
	fprintf(outfile,"#X Y Z\n");
	for(i = 0;i < numpoints;i++)
	{
		output[i].z = 0;
		for(j = 0;j < numpoints;j++)
		{
				output[i].z += input[j].z*cexp(2*M_PI*_Complex_I*(input[j].x*output[i].x+input[j].y*output[i].y));
		}

//		output[i].z/=(sqrt(2*PI)*500);

		//for normalizing
		if(cabs(output[i].z) > max)
			max = cabs(output[i].z);

	}

	for(i = 0;i < numpoints;i++)
	{
		output[i].z = output[i].z/max;
		fprintf(infile,"%2.8f %2.8f %2.8lf+i%2.8lf %2.8lf %2.8lf\n",input[i].x,input[i].y,creal(input[i].z),cimag(input[i].z),cabs(input[i].z),carg(input[i].z));
//		fprintf(infile,"%2.8f %2.8f %2.8lf+i%2.8lf %2.8lf %2.8lf\n",input[i].x,input[i].y,creal(input[i].z),cimag(input[i].z),cabs(input[i].z),\
		atan2(creal(input[i].z),cimag(input[i].z)));
	       fprintf(outfile,"%2.8f %2.8f %2.8lf+i%2.8lf %2.8lf %2.8lf\n",output[i].x,output[i].y,creal(output[i].z),cimag(output[i].z),cabs(output[i].z),carg(output[i].z));
//	       fprintf(outfile,"%2.8f %2.8f %2.8lf+i%2.8lf %2.8lf %2.8lf\n",output[i].x,output[i].y,creal(output[i].z),cimag(output[i].z),cabs(output[i].z),\
		atan2(creal(output[i].z),cimag(output[i].z)));
	}

	fclose(infile);
	fclose(outfile);
}

float eval_dft(dftdata input[], int numpoints, float x, float y)
{
	complex out;
	int j;
	out = 0;
	for(j = 0;j < numpoints;j++)
	{
		out += input[j].z*cexp(-2*M_PI*_Complex_I*(input[j].x*x+input[j].y*y));
	}
	return(cabsf(out));
}

static int factorial(int n)
{
	int i,fact;
	if(n==0)
		return 1;
	fact = 1;
	for(i = 1;i<=n;i++)
		fact*=i;
	return(fact);
}

void zernike_pt(int n, int m, float x, float y, float *z,float radius)
{
	int j,jmax;
	float r, theta, norm,max;
	jmax = 0.5*(n-abs(m));
	norm=sqrt(2*(n+1)/(1+(m==0)));
	
	max = 0.0;
	*z = 0;
	if(n == 0)
		*z = 1;
	else if((n-m)%2)
		*z = 0;
	else
	{
		r = sqrt(x*x+y*y);			
		if((x>=0 && y>=0) || (x>=0 && y<0))
		        theta = atan(y/x);
		else
			theta = M_PI + atan(y/x);
		for(j = 0;j <= jmax;j++)
		{
			*z+=pow(-1,j)*factorial(n-j)*pow(r/radius,n-2*j)/(factorial(j)*factorial(0.5*(n+abs(m))-j)*factorial(0.5*(n-abs(m))-j));
		}
	        *z=norm*(*z)*((m>=0)*cos(m*theta)-(m<0)*sin(m*theta));
	}
}

void zernike(int n, int m, float x[], float y[], int numpoints, float z[],float radius)
{
	int i,j,jmax;
	float r, theta, norm,max;
//	FILE * output;
//	char filename[64];
//	sprintf(filename,"zernike_%d_%d.dat",n,m);
//	output = fopen(filename,"w");
//	fprintf(output,"#x y znk\n");
	jmax = 0.5*(n-abs(m));
	norm=sqrt(2*(n+1)/(1+(m==0)));
	
	max = 0.0;
	for (i = 0;i < numpoints;i++)
	{
		z[i] = 0;
		if(n == 0)
			z[i] = 1;
		else if((n-m)%2)
			z[i] = 0;
		else
		{
			r = sqrt(x[i]*x[i]+y[i]*y[i]);			
			if((x[i]>=0 && y[i]>=0) || (x[i]>=0 && y[i]<0))
			        theta = atan(y[i]/x[i]);
			else
				theta = M_PI + atan(y[i]/x[i]);

			if(x[i]==0 && y[i]==0)
				theta = M_PI/2;

			for(j = 0;j <= jmax;j++)
			{
				z[i]+=pow(-1,j)*factorial(n-j)*pow(r/radius,n-2*j)/(factorial(j)*factorial(0.5*(n+abs(m))-j)*factorial(0.5*(n-abs(m))-j));
			}
		        z[i]=norm*z[i]*((m>=0)*cos(m*theta)-(m<0)*sin(m*theta));
		}
//		if(z[i]>max)
//			max = z[i];
//		fprintf(output,"%lf %lf %lf\n",x[i],y[i],z[i]);
	}
	//normalizing
//	for (i = 0;i < numpoints;i++)
//	{
//		z[i]/=max;
//		fprintf(output,"%lf %lf %lf\n",x[i],y[i],z[i]);
//	}
//	fclose(output);
}


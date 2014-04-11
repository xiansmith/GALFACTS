#include<stdio.h>
#include<string.h>
#include<stdlib.h>

void main()
{
	FILE *f,*e,*p,*o;
	f = fopen("epsilonphi.list","r");
	if(NULL == f)
	{
		printf("epsilonphi.list not found !\n");
		exit(1);
	}

	int i;
	for(i=0;i<7;i++)
	{
		char a[20],b[20];
		fscanf(f,"%s %s",&a,&b);
		if(!strncmp(a,b))
		{
			printf("Beam %d same source for both, copy manualy.\n",i);
			continue;
		}
		printf("%s %s\n",a,b);
		char ename[30],pname[30],oname[30];
		sprintf(ename,"%s/epsilon_phi%d.dat",a,i);
		sprintf(pname,"%s/epsilon_phi%d.dat",b,i);
		e = fopen(ename,"r");
		if(NULL == e)
		{
			printf("%s not found !\n",ename);
			exit(1);
		}
		p = fopen(pname,"r");
		if(NULL == p)
		{
			printf("%s not found !\n",pname);
			exit(1);
		}
		int j;
		sprintf(ename,"epsilon_phi%d.dat",i);
		o = fopen(oname,"w");
		if(NULL == o)
		{
			printf("%s not opened !\n",oname);
			exit(1);
		}
		for(j=0;j<4096;j++)
		{
			float eps,phi,tmp;
			fscanf(e,"%f %f",&eps,&tmp);
			fscanf(p,"%f %f",&tmp,&phi);
			fprintf(o,"%f %f\n",eps,phi);
			printf("%f %f\n",eps,phi);
		}
		fclose(e);
		fclose(p);
		fclose(o);
	}
}

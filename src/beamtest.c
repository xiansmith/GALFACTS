#include<stdio.h>
#include<math.h>
#include"jsd/jsd_futil.h"

void main(int argc,char * argv[])
{
	int beamno;
	FILE *infile,*outfile;
        if(argc != 2)
        {
                printf("Usage: beamtest <beamno>\n");
                //printf("Usage: beammodels <radius> <startchannel> <endchannel>  <RA> <DEC> \n");
                exit(1);
        }
        else
        {
                beamno = atoi(argv[1]);
        }
	char infilename[64],outfilename[64];
	int i,lines;
	double x,y,data,model,err,power;
	sprintf(infilename,"beam%d_model/xxbm%d.dat",beamno,beamno);
	sprintf(outfilename,"beam%d_model/xxbm%d_err.dat",beamno,beamno);
	printf("%s %s\n",infilename,outfilename);
	infile = fopen(infilename,"r");
	if(infile == NULL)
	{
		printf("cannot open %s\n",infilename);
		exit(1);
	}
	else
		printf("opened %s\n",infilename);
	outfile = fopen(outfilename,"w");
	if(outfile == NULL)
	{
		printf("cannot open %s\n",outfilename);
		exit(1);
	}
	else
		printf("opened %s\n",outfilename);
	lines = jsd_line_count(infile);
	printf("no of points %d\n",lines);

	for(i = 0;i < lines;i++)
	{
		fscanf(infile,"%lf %lf %lf",&x,&y,&data);
		power = -786.673*(x+0.00238159)*(x+0.00238159)-817.957*(y+0.0037449)*(y+0.0037449);
		printf("%2.8f\n",power);
		model = 2.06007*exp(power);
		err = data-model;
		fprintf(outfile,"%2.8f %2.8f %2.8f %2.8f %2.8f\n",x,y,err,model,data);
	}
	fclose(infile);
	fclose(outfile);
	
	sprintf(infilename,"beam%d_model/yybm%d.dat",beamno,beamno);
	sprintf(outfilename,"beam%d_model/yybm%d_err.dat",beamno,beamno);
	infile = fopen(infilename,"r");
	if(infile == NULL)
	{
		printf("cannot open %s\n",infilename);
		exit(1);
	}
	outfile = fopen(outfilename,"w");
	if(outfile == NULL)
	{
		printf("cannot open %s\n",outfilename);
		exit(1);
	}
	lines = jsd_line_count(infile);
	printf("no of points %d\n",lines);

	for(i = 0;i < lines;i++)
	{
		fscanf(infile,"%lf %lf %lf",&x,&y,&data);
		power = -795.694*(x+0.00347361)*(x+0.00347361)-818.206*(y+0.00385754)*(y+0.00385754);
		model = 2.0785*exp(power);
		err = data-model;
		fprintf(outfile,"%lf %lf %lf %lf\n",x,y,err,model);
	}
	fclose(infile);
	fclose(outfile);

}

#include "fluxdata.h"
#include "jsd/jsd_futil.h"
#include "common.h"

int main(int argc, char *argv[])
{
	char infilename[64],outfilename[64],beamstring[5+1],header[80+1];
	int i,j,channel,beam,numDays,count,region;
	float ramax,ramin,decmax,decmin;
	int startchannel,endchannel;
	FILE * infile,* outfile;
	FluxRecord pRec;
	char **files;

	printf("Starting program, constructing subregion files....\n");
	if(argc != 8)
	{
		printf("Usage: get_subregion <startchannel> <endchannel> <ramin> <ramax> <decmin> <decmax> <#region>\n");
		return EXIT_FAILURE;
	}
	else
	{
		startchannel = atoi(argv[1]);
		endchannel = atoi(argv[2]);
		//ramin = atof(argv[3])*15;
		ramin = atof(argv[3]) -0.1;
		//ramax = atof(argv[4])*15;
		ramax = atof(argv[4]) +0.1;
		decmin = atof(argv[5])-0.1;
		decmax = atof(argv[6])+0.1;
		region = atoi(argv[7]);
	}
	numDays = get_date_dirs("./", &files);
//	numDays = get_date_dirs("/n/ras1/people/jeff/run8/", &files);

	for(channel = startchannel;channel < endchannel;channel++)
	{		
		for(j = 0;j < numDays;j++)
		{
			//for(beam = 0;beam < 7;beam++)
			for(beam = 0;beam < 6;beam++)
			{
	
				sprintf(beamstring,"beam%d",beam);
				//sprintf(infilename,"/n/ras1/people/jeff/run8/%s/%s/balance%03i.dat",files[j],beamstring,channel);
				sprintf(infilename,"%s/%s/balance%04i.dat",files[j],beamstring,channel);
				printf("Processing %s.\n",infilename);
				infile = fopen(infilename, "r");
			        if (infile == NULL) {
        				printf("ERROR: can't open input file %s\n", infilename);
			        	continue;
		        	}
				sprintf(outfilename,"%s/%s/subregion%04i_%02i.dat",files[j],beamstring,channel,region);
	
				outfile = fopen(outfilename, "w");
				if (outfile == NULL) {
        				printf("ERROR: can't open output file %s\n", outfilename);
			        	continue;
				}
				else
					fprintf(outfile,"#RA DEC AST I Q U V\n");

				count = jsd_line_count(infile);
		        	fgets(header, 80, infile);
	
				for(i = 0;i < count;i++)
				{
					fscanf(infile,"%f %f %f %lf %lf %lf %lf",&pRec.RA,&pRec.DEC,&pRec.AST,&pRec.stokes.I,&pRec.stokes.Q,&pRec.stokes.U,\
					&pRec.stokes.V);
					if(pRec.RA < ramax && pRec.RA >= ramin && pRec.DEC < decmax && pRec.DEC >= decmin)
						fprintf(outfile,"%2.8f %2.8f %f %lf %lf %lf %lf\n",pRec.RA,pRec.DEC,pRec.AST,pRec.stokes.I,pRec.stokes.Q,\
						pRec.stokes.U,pRec.stokes.V);
				}
			        fclose(infile);
				fclose(outfile);			
				printf("Wrote %s.\n",outfilename);
			}//beam loop
		}

	}
	printf("Done.\n");
	return EXIT_SUCCESS;
}

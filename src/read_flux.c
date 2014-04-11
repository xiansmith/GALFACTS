#include "read_flux.h"
#include "jsd/jsd_futil.h"
#include <string.h>
#include <values.h>
#include <math.h>

static int fluxrecord_read(FluxRecord * pRec, FILE * file)
{
    	return fscanf(file,"%f %f %f %lf %lf %lf %lf",&pRec->RA, &pRec->DEC, &pRec->AST,&pRec->stokes.I, &pRec->stokes.Q, &pRec->stokes.U, &pRec->stokes.V);

}


static int fluxdata_read(FluxData *dataset, FILE *infile, int beam)
{
    	int k;
    	int numRecords;
    	char header[80+1];
    	float RAmax = FLT_MIN;
    	float RAmin = FLT_MAX;
    	
	// Each line is one record, count the number of lines to 
	// determine the number of Records, allocate the dataser record array
    	numRecords = jsd_line_count(infile)-1;
	printf("There are %d records\n",numRecords);
	dataset->records = (FluxRecord*) malloc(numRecords * sizeof(FluxRecord));

	if(dataset->records == NULL)
	{
		printf("ERROR: malloc failed !\n");
		exit(0);
	}

	//First line of each file is just the column names
    	// read out the header on the fluxtime files
    	fgets(header, 80, infile);

    	k = 0;
	//loop through the whole file ewading a line at a time
    	while (!feof(infile) && k<numRecords)
    	{
		FluxRecord *pRec = &dataset->records[k];
        	int num = fluxrecord_read(pRec, infile);
		//A correct line in the file should have 7 items RA,DEC,AST,I,Q,U,V
        	if (num == 7) 
		{
			//Determine the RA Max and Min values
	            	float RA = dataset->records[k].RA;
        	    	if (RA > RAmax) RAmax = RA;
            		if (RA < RAmin) RAmin = RA;
			
			//The IGALFA WAPP data has pointing errors due to incorrect rotation angles
			//in first stage processing. These error should ideally be caculated during the cal
			//program, but currently they are propagated through and are corrected for here.
			//The MOCK data for first run also has pointing problems. Below are the correction
			//factors for each case. Please uncomment the appropriate correction factors
			//before running for each beam.
			switch(beam)
			{
				case 0:
		        		break;
		    		case 1:
					//MOCK
	        			//dataset->records[k].RA -= 0.05086899;
		        		//dataset->records[k].DEC -= 0.09247589;
		        		//A Field
	        			//dataset->records[k].RA += 0.04854810;
		        		//dataset->records[k].DEC += 0.01256719;
		        		//B Field 
	        			//dataset->records[k].RA += 0.00432813;
		        		//dataset->records[k].DEC += 0.05169168;
				        //C Field
				        //dataset->records[k].RA -= 0.00390088;
				        //dataset->records[k].DEC -= 0.05169645;
				        //D Field
			        	dataset->records[k].RA -= 0.05422437;
				        dataset->records[k].DEC -= 0.01237836;
				        //ZW Field
			        	//dataset->records[k].RA -= 0.10672534;
				        //dataset->records[k].DEC += 0.08725516;
			        	break;
				case 2:
			        	//dataset->records[k].RA -= 0.10275722;
				        //dataset->records[k].DEC -= 0.00045204;
			        	//dataset->records[k].RA += 0.00427891;
				        //dataset->records[k].DEC += 0.05073190;
			        	//dataset->records[k].RA -= 0.08580898;
				        //dataset->records[k].DEC += 0.10434259;
			        	//dataset->records[k].RA += 0.08623622;
				        //dataset->records[k].DEC -= 0.10433782;
	        			dataset->records[k].RA -= 0.01685216;
				        dataset->records[k].DEC -= 0.05022825;
			        	//dataset->records[k].RA -= 0.12859462;
				        //dataset->records[k].DEC -= 0.03564440;
	        			break;
				case 3:
			        	//dataset->records[k].RA -= 0.05189181;
				        //dataset->records[k].DEC += 0.09227562;
			        	//dataset->records[k].RA -= 0.04443035;
				        //dataset->records[k].DEC += 0.03798867;
			        	//dataset->records[k].RA -= 0.09032879;
				        //dataset->records[k].DEC += 0.05237175;
			        	//dataset->records[k].RA += 0.09032879;
				        //dataset->records[k].DEC -= 0.05235840;
			        	dataset->records[k].RA += 0.03753338;
				        dataset->records[k].DEC -= 0.03779388;
	        			//dataset->records[k].RA -= 0.02193891;
				        //dataset->records[k].DEC -= 0.12381958;
	        			break;
				case 4:
			        	//dataset->records[k].RA += 0.05200982 ;
				        //dataset->records[k].DEC += 0.09222222;
	        			//dataset->records[k].RA -= 0.04753914;
				        //dataset->records[k].DEC -= 0.01311164;
			        	//dataset->records[k].RA -= 0.00295296;
				        //dataset->records[k].DEC -= 0.05213325;
			        	//dataset->records[k].RA += 0.00249520;
				        //dataset->records[k].DEC += 0.05213802;
			        	dataset->records[k].RA += 0.05388680;
				        dataset->records[k].DEC += 0.01303928;
			        	//dataset->records[k].RA += 0.10806999;
				        //dataset->records[k].DEC -= 0.08763950;
			        	break;
		    		case 5:
			        	//dataset->records[k].RA += 0.10275006;
				        //dataset->records[k].DEC -= 0.00055885;
			        	//dataset->records[k].RA -= 0.00363614;
				        //dataset->records[k].DEC -= 0.05003117;
			        	//dataset->records[k].RA += 0.08645175;
				        //dataset->records[k].DEC -= 0.10322927;
			        	//dataset->records[k].RA -= 0.08690952;
				        //dataset->records[k].DEC += 0.10322641;
			        	dataset->records[k].RA += 0.01623990;
				        dataset->records[k].DEC += 0.05000566;
			        	//dataset->records[k].RA += 0.12853549;
				        //dataset->records[k].DEC += 0.03735663;
			        	break;
		    		case 6:
			        	//dataset->records[k].RA += 0.05074024;
				        //dataset->records[k].DEC -= 0.09252929;
			        	//dataset->records[k].RA += 0.04377896;
				        //dataset->records[k].DEC -= 0.03709938;
			        	//dataset->records[k].RA += 0.08925015;
				        //dataset->records[k].DEC -= 0.05137899;
	       			 	//dataset->records[k].RA -= 0.08925015;
				        //dataset->records[k].DEC += 0.05137136;
			        	dataset->records[k].RA -= 0.03743130;
				        dataset->records[k].DEC += 0.03702619;
	        			//dataset->records[k].RA += 0.02046353;
				        //dataset->records[k].DEC += 0.12407229;
	        			break;
				default:
	        			printf("WARN: invalid beam (%i) specified.\n", beam);
			}
		}
        	else if (num < 0)
		{
			printf("ERROR %d while reading record\n",num);
			break;
		}
        	else{
			printf("ERROR: flux file record only read %i fields\n", num);
			exit(1);
		}
		k++;	
    	}

    	dataset->numRecords = k;
    	dataset->RAmin = RAmin;
    	dataset->RAmax = RAmax;
    	return k;
}

static void print_fluxdata(FluxData dataset,int start,int end)
{
	int i;
	printf("RA	   DEC	     AST      I         Q         U        V\n");
	for(i=start;i<end;i++)
	{
		FluxRecord *pRec = &dataset.records[i];
		printf("%2.6f %2.6f %6.2f %2.6f %2.6f %2.6f %2.6f\n",pRec->RA,pRec->DEC,pRec->AST,\
		pRec->stokes.I,pRec->stokes.Q,pRec->stokes.U,pRec->stokes.V);
	}
}

void main(int argc, char *argv[])
{
        char fluxfilename[32+1],*day;
        int numRecords,beam,chan;
        FluxData dataset;
        FILE *fluxfile;
        if (argc != 4) 
	{
                printf("Usage: %s <day> <beam> <channel>\n", argv[0]);
                exit(1);
        } 
        else
        { 
                day = argv[1];
                beam = atoi(argv[2]);
                chan = atoi(argv[3]);
        }
	sprintf(fluxfilename,"%s/beam%d/fluxtime%04i.dat",day,beam,chan);
        printf("Reading data file %s. \n", fluxfilename);
        fluxfile = fopen(fluxfilename,"r");
        numRecords = fluxdata_read(&dataset, fluxfile, beam);
        fclose(fluxfile);
        printf("Read %i records\n", numRecords);
        print_fluxdata(dataset, 0, 22);
}


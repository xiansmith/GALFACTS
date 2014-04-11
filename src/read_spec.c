#include "read_spec.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

/*
Helper Function Declarations
*/
static int fread_record(SpecRecord * rec, FILE * datafile, int beam);



/*
Read the .cfg metadata from file into the ConfigData structure.
@param pFile - open file handle in binary mode for reading the data
@param pCfg - pointer to the pre-allocated object to store the config data
@param beam - the beam number to read for [0..6]
@return true if successful
*/
int read_cfgfile(FILE * pFile, ConfigData * pCfg)
{
	int num = 0;
	num += fscanf(pFile, "%f\n", &(pCfg->integrationTime)); 
	num += fscanf(pFile, "%i\n", &(pCfg->stokesSetSize)); 
	num += fscanf(pFile, "%f\n", &(pCfg->centerMHz)); 
	num += fscanf(pFile, "%f\n", &(pCfg->bandwitdhkHz)); 
	num += fscanf(pFile, "%i%i%i%i%i\n", &(pCfg->startChanNum), &(pCfg->samplesPerStokesSet), 
										&(pCfg->stokesProducts), &(pCfg->numStokesSets), &(pCfg->crap)); 
	num += fscanf(pFile, "%s\n", pCfg->observingTag); 
	num += fscanf(pFile, "%i%i%i\n", &(pCfg->day), &(pCfg->month), &(pCfg->year)); 
	num += fscanf(pFile, "%i%i%f\n", &(pCfg->hour), &(pCfg->minute), &(pCfg->second)); 
	num += fscanf(pFile, "%s\n", pCfg->observatoryCode); 
	num += fscanf(pFile, "%i\n", &(pCfg->bandFlip)); 			


	return (num == 18);
}

/*
pDataset is malloced and loaded from file
number of records is returned
pFile must be open for reading
caller is responsible for freeing the dataset
@param pFile - open file handle in binary mode for reading the data
@param pDataset - pointer to the array to store the spec records in.  malloc is used
to allocate these data sturcutres and must be externally freed.
@param beam - the beam number to read for [0..6]
@return the number of records read
*/
int read_datafile(FILE * pFile, SpecRecord ** pDataset, int beam)
{
	int numRead;
	int numExpected;
	long fileSize;

	assert(sizeof(SpecPointingBlock) == SIZEOF_SPECPOINTINGBLOCK);

//  printf("Using pointing from beam %i\n", beam);

	// obtain file size.			for (i=0; i<numRecords; i++) {
	fseek (pFile , 0 , SEEK_END);
	fileSize = ftell (pFile);
	rewind (pFile);

	// malloc the memory for the expected umber of records
	numExpected = fileSize / (512 + ((MAX_CHANNELS*4*4) * 2));
//	printf("Requesting malloc for %ld bytes of memory\n",sizeof(SpecRecord)*numExpected);
	//*pDataset = (SpecRecord *)malloc(sizeof(SpecRecord) * numExpected);
	*pDataset = (SpecRecord *)malloc(sizeof(SpecRecord) * 300);
	if (*pDataset == NULL) {
		printf("ERROR: malloc failed in read_datafile() !\n");
		return 0;
	}

	// read from file to memory
	numRead = 0;
	do {
		SpecRecord * rec = &((*pDataset)[numRead]);
		numRead += fread_record(rec, pFile, beam);
		rec->RA *= 15.0; //convert to degrees
	} while(!feof(pFile) && numRead < 300);
	//} while(!feof(pFile) && numRead < numExpected);

	//return results
	return numRead;
}

/*
Read a single SpecRecord from disk including the pointing block
and the raw spectral data.  
  
@param rec the structure to read the data into.
@param datafile the open file for reading
@param beam the beam number to read the pointing info from
@return 1 if record read successfully, 0 otherwise
*/
static int fread_record(SpecRecord * rec, FILE * datafile, int beam)
{
	int itemsRead = 0;
	SpecPointingBlock block;

    //read and handle the header block
    itemsRead += fread(&block, sizeof(SpecPointingBlock), 1, datafile);

    switch (beam) {
    case 0:
        rec->RA = block.centralBeam.raj_true_in_hours;
        rec->DEC = block.centralBeam.decj_true_in_degrees;
        break;
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
        rec->RA = block.outerBeams[beam-1].raj_true_in_hours;
        rec->DEC = block.outerBeams[beam-1].decj_true_in_degrees;
        break;
    default:
        printf("WARN: invalid beam (%i) specified.\n", beam);
    }
    
    rec->AST = block.centralBeam.atlantic_solar_time_now_in_sec;

//fix for bad beam6 days commit ???
/*    int i;
    if(beam == 6)
    {
	    for(i = 0; i < MAX_CHANNELS;i++)
	    {
			rec->calon.yy[i] = rec->calon.xx[i];
			rec->caloff.yy[i] = rec->caloff.xx[i];
	    }	
    }*/
//fix end

    //read the data
    //itemsRead += fread(&(rec->cal),  sizeof(PolSet), 1, datafile);
    itemsRead += fread(&(rec->calon), sizeof(PolSet), 1, datafile);
    itemsRead += fread(&(rec->caloff), sizeof(PolSet), 1, datafile);

    //initialze the rest of the record
	memset(rec->flagRFI, RFI_NONE, MAX_CHANNELS);
	rec->flagBAD = 0;

    //return result
	if (itemsRead == 3) {
		return 1; //successfully read 1 record
	} else if (itemsRead == 0) {
		return 0; //no more records to read
	} else { //partial read of a record
		printf("ERROR: partial record read %i items\n", itemsRead);
		return 0;
	}
}

/*
Write the spectral data to file.  For testing purposes only.  Output is
in plain text rows and columns.
@param array of SpecRecords to write
@param numRecords the number of records in the dataset
@param lowchan the low channel 
@param highchan the high channel
@param freq an array of frequencies which map to the channel numbers
*/
void print_data(SpecRecord dataset[], int numRecords, int lowchan, int highchan)
{
	int i,j;
	SpecRecord * pRec;
	FILE * file;

	file = fopen("specfile.dat", "w");
	fprintf(file, "#chan RA DEC AST xx yy xy yx\n");

	for (i=0; i<numRecords; i++)
	{ 
		pRec = &(dataset[i]);
		for (j=lowchan; j<highchan; j++) 
		{
			fprintf(file,"%4d %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f\n",
					j, pRec->RA, pRec->DEC,pRec->AST, pRec->calon.xx[j], pRec->calon.yy[j],
					pRec->calon.xy[j], pRec->calon.yx[j]);
		}
	}
	fclose(file);
}


void main(int argc, char *argv[])
{
	char *specfilename;
	int numRecords,beam;
	SpecRecord * dataset;
	FILE *specfile;
	if (argc != 3) {
		printf("Usage: %s <specfilename> <beam>\n", argv[0]);
		exit(1);
	}
	else
	{
		specfilename = argv[1];
		beam = atoi(argv[2]);
	}
	printf("Reading data file %s. \n", specfilename);
	specfile = fopen(specfilename,"r");
	numRecords = read_datafile(specfile, &dataset, beam);
	fclose(specfile);
	printf("Read %i records\n", numRecords);
	//print_data(dataset, numRecords, 0, 128);
	print_data(dataset, numRecords, 1000, 1001);
}

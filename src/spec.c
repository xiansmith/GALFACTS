#include "spec.h"
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

    printf("Using pointing from beam %i\n", beam);

	// obtain file size.			for (i=0; i<numRecords; i++) {
	fseek (pFile , 0 , SEEK_END);
	fileSize = ftell (pFile);
	rewind (pFile);

	// malloc the memory for the expected umber of records
	numExpected = fileSize / (512 + ((MAX_CHANNELS*4*4) * 2));
	*pDataset = (SpecRecord *)malloc(sizeof(SpecRecord) * numExpected);
	if (*pDataset == NULL) {
		printf("ERROR: malloc failed!\n");
		return 0;
	}

	// read from file to memory
	numRead = 0;
	do {
		SpecRecord * rec = &((*pDataset)[numRead]);
		numRead += fread_record(rec, pFile, beam);
		rec->RA *= 15.0; //convert to degrees
	} while(!feof(pFile) && numRead < numExpected);

	//return results
	return numRead;
}

/*
Read a single SpecRecord from disk including the pointing block
and the raw spectral data.  
TODO: Due to a bug in the Arecibo software, the pointing information in the
header is systematically wrong.  There is a hack in this code to correct for
this problem.  After it is fixed in Arecibo, the hack must be removed here.
  
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

    //read the data
    itemsRead += fread(&(rec->calon),  sizeof(PolSet), 1, datafile);
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
void print_data(SpecRecord dataset[], int numRecords, int lowchan, int highchan, float freq[])
{
	int i,j;
	SpecRecord * pRec;
	FILE * file;

	file = fopen("datafile.dat", "w");
	fprintf(file, "#chan freq AST offxx offyy onxx onyy\n");

	for (i=0; i<numRecords; i++)
	{ 
		pRec = &(dataset[i]);
		for (j=lowchan; j<highchan; j++) 
		{
			fprintf(file,"%4d %9.3f %8.2f %8.3f %8.3f %8.3f %8.3f\n", 
					j, freq[j], pRec->AST, pRec->caloff.xx[j], pRec->caloff.yy[j], 
					pRec->calon.xx[j], pRec->calon.yy[j]);
		}
	}
	fclose(file);
}



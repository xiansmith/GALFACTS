#include "fluxdata.h"
#include "jsd/jsd_futil.h"
#include <string.h>
#include <values.h>

int fluxrecord_read(FluxRecord * pRec, FILE * file)
{
    return fscanf(file,"%f %f %f %lf %lf %lf %lf",
            &pRec->RA, &pRec->DEC, &pRec->AST,
            &pRec->stokes.I, &pRec->stokes.Q, &pRec->stokes.U, &pRec->stokes.V);

}

int fluxrecord_write(FluxRecord * pRec, FILE * file)
{
    return fprintf(file,"%2.8f %2.8f %8.2f %4.6f %4.6f %4.6f %4.6f\n",
            pRec->RA, pRec->DEC, pRec->AST,
            pRec->stokes.I, pRec->stokes.Q, pRec->stokes.U, pRec->stokes.V);
}

FluxWappData * fluxwappdata_alloc(const char *wapp, char **days, int numDays)
{
    int i;
    FluxWappData * wappdata;

    wappdata = (FluxWappData*) malloc(sizeof(FluxWappData));
    strncpy(wappdata->wapp, wapp, WAPP_LEN);
    wappdata->numDays = numDays;
    wappdata->daydata = (FluxDayData*) malloc(sizeof(FluxDayData) * numDays);
    wappdata->scanDayData = (ScanDayData*) malloc(sizeof(ScanDayData) * numDays);
    for (i=0; i<numDays; i++) {
        strncpy(wappdata->daydata[i].mjd, days[i], MJD_LEN);
        wappdata->daydata[i].numRecords = 0;
        wappdata->daydata[i].records = NULL;
        wappdata->scanDayData[i].numScans = 0;
        wappdata->scanDayData[i].scans = NULL;
    }

    return wappdata;
}

void fluxwappdata_free(FluxWappData * wappdata)
{
    int i;

    if (wappdata == NULL) return;

    if (wappdata->daydata != NULL) {
        for (i=0; i<wappdata->numDays; i++) {
            if (wappdata->daydata[i].records != NULL) {
                free(wappdata->daydata[i].records);
                wappdata->daydata[i].records = NULL;
            }
        }
        free(wappdata->daydata);
        wappdata->daydata = NULL;
    }
    free(wappdata);
}

int fluxwappdata_writeavg(FluxWappData * wappdata)
{
    int m;
    int count;

    count = 0;
    for (m=0; m<wappdata->numDays; m++)
    {
        int k;
        FILE *file;
        int numRecords;
        char filename[64+1];

        FluxDayData * daydata = &wappdata->daydata[m];
        sprintf(filename, "%s/%s/balance.dat", daydata->mjd, wappdata->wapp);
        file = fopen(filename, "w");
        if (file == NULL) {
            printf("ERROR: can't open output file %s\n", filename);
            continue;
        }

        numRecords = daydata->numRecords;
        for (k=0; k<numRecords; k++) {
            fluxrecord_write(&daydata->records[k], file);
        }

        fclose(file);
        count++;
    }
    return count;
}


int fluxwappdata_writechan(FluxWappData * wappdata, int chan)
{
    int m;
    int count;

    count = 0;
    for (m=0; m<wappdata->numDays; m++)
    {
        int k;
        FILE *file;
        int numRecords;
        char filename[64+1];

        FluxDayData * daydata = &wappdata->daydata[m];
        sprintf(filename, "%s/%s/balance%03i.dat", daydata->mjd, wappdata->wapp, chan);
        file = fopen(filename, "w");
        if (file == NULL) {
            printf("ERROR: can't open output file %s\n", filename);
            continue;
        }

        numRecords = daydata->numRecords;
        for (k=0; k<numRecords; k++) {
            fluxrecord_write(&daydata->records[k], file);
        }

        fclose(file);
        count++;
    }
    return count;
}


static int fluxdaydata_read(FluxDayData *daydata, FILE *infile)
{
    int k;
    int numRecords;
    char header[80+1];
    float RAmax = FLT_MIN;
    float RAmin = FLT_MAX;

    // allocate this days record array
    numRecords = jsd_line_count(infile);
    if (daydata->records != NULL) {
        free(daydata->records);
    }
printf("allocating %i FluxRecords\n", numRecords);
    daydata->records = (FluxRecord*) malloc(numRecords * sizeof(FluxRecord));

    // read out the # header on the fluxtime files
    fgets(header, 80, infile);

    k = 0;
    while (!feof(infile) && k<numRecords)
    {
        int num = fluxrecord_read(&daydata->records[k], infile);
        if (num == 7) {
            float RA = daydata->records[k].RA;
            if (RA > RAmax) RAmax = RA;
            if (RA < RAmin) RAmin = RA;
            k++;
        }
        else if (num < 0) break;
        else printf("ERROR: flux file record only read %i fields\n", num);
    }
    daydata->numRecords = k;
    daydata->RAmin = RAmin;
    daydata->RAmax = RAmax;

    return numRecords;
}




//----------- read in data from input flux files
//dataset must be allocated of size ndays
//dataset[].records will be malloced new memory
int fluxwappdata_readchan(FluxWappData * wappdata, int chan)
{
    int  m;
    int count;
    FILE *infile;

    count = 0;
    for (m=0; m<wappdata->numDays; m++)
    {
        char filename[64+1];
        FluxDayData * daydata = &wappdata->daydata[m];

        sprintf(filename, "%s/%s/fluxtime%03i.dat", daydata->mjd, wappdata->wapp, chan);
        infile = fopen(filename, "r");
        if (infile == NULL) {
            printf("ERROR: can't open input file %s\n", filename);
            continue;
        }

        fluxdaydata_read(daydata, infile);
        fclose(infile);
        count++;

    }

    return count;
}





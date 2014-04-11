#include "fluxdata.h"
#include "jsd/jsd_futil.h"
#include <string.h>
#include <values.h>
#include <math.h>

extern int multibeam;//SSG
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
    int i,j;
    FluxWappData * wappdata;

    wappdata = (FluxWappData*) malloc(sizeof(FluxWappData));
    strncpy(wappdata->wapp, wapp, WAPP_LEN);
    wappdata->numDays = numDays;
//	printf("Requesting malloc for %u bytes\n",sizeof(FluxDayData) * numDays);
    wappdata->daydata = (FluxDayData*) malloc(sizeof(FluxDayData) * numDays);
//	printf("Requesting malloc for %u bytes\n",sizeof(ScanDayData) * numDays);
    wappdata->scanDayData = (ScanDayData*) malloc(sizeof(ScanDayData) * numDays);
    for (i=0; i<numDays; i++) {
		if(!strcmp(wapp,"multibeam"))//SSG
		{
			j = i/7;
	       	strncpy(wappdata->daydata[i].mjd, days[j], MJD_LEN);//SSG
		}
		else
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
		else //SSG
			fprintf(file, "#RA DEC AST I Q U V\n"); //SSG
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
	char tempstring[6];
        FluxDayData * daydata = &wappdata->daydata[m];
	//SSG
	if(multibeam)
	{
		sprintf(tempstring,"beam%d",m%7);
	        //sprintf(filename, "%s/%s/balance%03i.dat", daydata->mjd, tempstring, chan);
	        sprintf(filename, "%s/%s/clean%04i.dat", daydata->mjd, tempstring, chan);
	}
	else
	        //sprintf(filename, "%s/%s/balance%03i.dat", daydata->mjd, wappdata->wapp, chan);
	        sprintf(filename, "%s/%s/balanceraw%04i.dat", daydata->mjd, wappdata->wapp, chan);
	//SSG
        file = fopen(filename, "w");
        if (file == NULL) {
            printf("ERROR: can't open output file %s\n", filename);
            continue;
        }
	else //SSG
		fprintf(file, "#RA DEC AST I Q U V\n"); //SSG

        numRecords = daydata->numRecords;
        for (k=0; k<numRecords; k++) {
            fluxrecord_write(&daydata->records[k], file);
        }

        fclose(file);
        count++;
    }
    return count;
}


static int fluxdaydata_read(FluxDayData *daydata, FILE *infile,int beam)
{
    int k;
    int numRecords;
    char header[80+1];
    float RAmax = FLT_MIN;
    float RAmin = FLT_MAX;
    // allocate this days record array
    numRecords = jsd_line_count(infile);

//    printf("has %d records\n",numRecords);
    if (daydata->records != NULL) {
        free(daydata->records);
    }

//	printf("Requesting malloc for %u bytes\n",numRecords * sizeof(FluxRecord));
    daydata->records = (FluxRecord*) malloc(numRecords * sizeof(FluxRecord));

	if(daydata->records == NULL)
	{
		printf("ERROR: malloc failed !\n");
		exit(0);
	}

    // read out the # header on the fluxtime files
    fgets(header, 80, infile);

	//pointing fix once and for all dec > 18
	//static float RAp[200000],DECp[200000],ASTp[200000];


    k = 0;
   int flag = 0;
	int m = 0;
    while (!feof(infile) && k<numRecords)
    {
		FluxRecord *pRec = &daydata->records[k];

        int num = fluxrecord_read(pRec, infile);
	
	//fix for certain missing data in beam 0 (a real pain!)
/*	if(ASTp[m] != daydata->records[k].AST && beam != 0)
		continue;

	//pointing fix once and for all dec > 18
	if(beam == 0)
	{
		RAp[m] = pRec->RA;
		DECp[m] = pRec->DEC;
		ASTp[m] = pRec->AST;
		m++;
		if(m > 200000)
			printf("Pointing fix array overrun !\n");
	}
*/
        if (num == 7) {
            float RA = daydata->records[k].RA;
  //          float DECR = daydata->records[k].DEC*M_PI/180.0;
            if (RA > RAmax) RAmax = RA;
            if (RA < RAmin) RAmin = RA;
/*		switch(beam)
		{
		    case 0:
		        break;
		    case 1:
                        daydata->records[k].RA = RAp[m] + 2.7417/(60*cos(DECR));
			m++;
			//printf("%f %f\n",DECR,daydata->records[k].RA);
			//MOCK
	        	//daydata->records[k].RA += 0.05086899;
		        //daydata->records[k].DEC -= 0.09247589;
		        //A Field
	        	//daydata->records[k].RA += 0.04854810;
		        //daydata->records[k].DEC += 0.01256719;
		        //B Field ?
	        	//daydata->records[k].RA += 0.00432813;
		        //daydata->records[k].DEC += 0.05169168;
		        //C Field
		        //daydata->records[k].RA -= 0.00390088;
		        //daydata->records[k].DEC -= 0.05169645;
		        //D Field
	        	//daydata->records[k].RA -= 0.05422437;
		        //daydata->records[k].DEC -= 0.01237836;
		        //ZW Field
	        	//daydata->records[k].RA -= 0.10672534;
		        //daydata->records[k].DEC += 0.08725516;
	        	break;
		    case 2:
                        daydata->records[k].RA =RAp[m]+ 5.4833/(60*cos(DECR));
			m++;
	        	//daydata->records[k].RA += 0.10275722;
		        //daydata->records[k].DEC -= 0.00045204;
	        	//daydata->records[k].RA += 0.00427891;
		        //daydata->records[k].DEC += 0.05073190;
	        	//daydata->records[k].RA -= 0.08580898;
		        //daydata->records[k].DEC += 0.10434259;
	        	//daydata->records[k].RA += 0.08623622;
		        //daydata->records[k].DEC -= 0.10433782;
	        	//daydata->records[k].RA -= 0.01685216;
		        //daydata->records[k].DEC -= 0.05022825;
	        	//daydata->records[k].RA -= 0.12859462;
		        //daydata->records[k].DEC -= 0.03564440;
	        	break;
		    case 3:
                        daydata->records[k].RA =RAp[m]+ 2.7417/(60*cos(DECR));
			m++;
	        	//daydata->records[k].RA += 0.05189181;
		        //daydata->records[k].DEC += 0.09227562;
	        	//daydata->records[k].RA -= 0.04443035;
		        //daydata->records[k].DEC += 0.03798867;
	        	//daydata->records[k].RA -= 0.09032879;
		        //daydata->records[k].DEC += 0.05237175;
	        	//daydata->records[k].RA += 0.09032879;
		        //daydata->records[k].DEC -= 0.05235840;
	        	//daydata->records[k].RA += 0.03753338;
		        //daydata->records[k].DEC -= 0.03779388;
	        	//daydata->records[k].RA -= 0.02193891;
		        //daydata->records[k].DEC -= 0.12381958;
	        	break;
		    case 4:
                        daydata->records[k].RA =RAp[m] -2.7417/(60*cos(DECR));
			m++;
	        	//daydata->records[k].RA -= 0.05200982 ;
		        //daydata->records[k].DEC += 0.09222222;
	        	//daydata->records[k].RA -= 0.04753914;
		        //daydata->records[k].DEC -= 0.01311164;
	        	//daydata->records[k].RA -= 0.00295296;
		        //daydata->records[k].DEC -= 0.05213325;
	        	//daydata->records[k].RA += 0.00249520;
		        //daydata->records[k].DEC += 0.05213802;
	        	//daydata->records[k].RA += 0.05388680;
		        //daydata->records[k].DEC += 0.01303928;
	        	//daydata->records[k].RA += 0.10806999;
		        //daydata->records[k].DEC -= 0.08763950;
	        	break;
		    case 5:
                        daydata->records[k].RA =RAp[m] -5.4833/(60*cos(DECR));
			m++;
	        	//daydata->records[k].RA -= 0.10275006;
		        //daydata->records[k].DEC -= 0.00055885;
	        	//daydata->records[k].RA -= 0.00363614;
		        //daydata->records[k].DEC -= 0.05003117;
	        	//daydata->records[k].RA += 0.08645175;
		        //daydata->records[k].DEC -= 0.10322927;
	        	//daydata->records[k].RA -= 0.08690952;
		        //daydata->records[k].DEC += 0.10322641;
	        	//daydata->records[k].RA += 0.01623990;
		        //daydata->records[k].DEC += 0.05000566;
	        	//daydata->records[k].RA += 0.12853549;
		        //daydata->records[k].DEC += 0.03735663;
	        	break;
		    case 6:
                        daydata->records[k].RA =RAp[m] -2.7417/(60*cos(DECR));
			m++;
	        	//daydata->records[k].RA -= 0.05074024;
		        //daydata->records[k].DEC -= 0.09252929;
	        	//daydata->records[k].RA += 0.04377896;
		        //daydata->records[k].DEC -= 0.03709938;
	        	//daydata->records[k].RA += 0.08925015;
		        //daydata->records[k].DEC -= 0.05137899;
	        	//daydata->records[k].RA -= 0.08925015;
		        //daydata->records[k].DEC += 0.05137136;
	        	//daydata->records[k].RA -= 0.03743130;
		        //daydata->records[k].DEC += 0.03702619;
	        	//daydata->records[k].RA += 0.02046353;
		        //daydata->records[k].DEC += 0.12407229;
	        	break;
		    default:
	        	printf("WARN: invalid beam (%i) specified.\n", beam);
		}


			//if(k > 0 && daydata->records[k].DEC < 37.6 && daydata->records[k].DEC > 19.7  && (daydata->records[k].RA - daydata->records[k-1].RA) < 0)

			//	printf("%2.8f %2.8f %2.8f\n",daydata->records[k].RA,daydata->records[k].DEC,daydata->records[k].AST);		
	    int moonflag = 0;	
            if (isfinite(pRec->stokes.I))
            //if (pRec->stokes.I < 100.0)
	    {
		if(!strcmp(daydata->mjd,"54811") && (daydata->records[k].AST > 9505)  && (daydata->records[k].AST < 9580))
			moonflag = 1;
		if(!strcmp(daydata->mjd,"54812") && (daydata->records[k].AST > 13475)  && (daydata->records[k].AST < 13555))
			moonflag = 1;
	//	if(!strcmp(daydata->mjd,"54817") && (daydata->records[k].RA < 19.18)) //discontinuity not moon
	//		moonflag = 1;

//		Remove alien spacecraft !!
		if(k > 0)
		{
			if(daydata->records[k].DEC < 37.6001 && daydata->records[k].DEC > 19.6999  && fabs(daydata->records[k].DEC - daydata->records[k-1].DEC) > 0.0003 && (daydata->records[k].RA - daydata->records[k-1].RA) > 0 && !moonflag)		
	           	k++;
		}
		if(k==0)
			k++;
	    }
*/
     	    //printf("RA %2.8f  DEC %2.8f I % 2.8f\n",daydata->records[k].RA,daydata->records[k].DEC,daydata->records[k].stokes.I);
            if (isfinite(pRec->stokes.I))
	    k++;

//hack for testing
		//pRec->stokes.I = 0.0;
		//pRec->stokes.Q = 0.0;
		//pRec->stokes.U = 0.0;
		//pRec->stokes.V = 0.0;

//hack for testing


//          if((pRec->stokes.I+pRec->stokes.U) > 1000 && flag ==0)
//          {
//            	printf("Stokes value too large %f\n",pRec->stokes.I);
//           	flag = 1;
//	    }
//            if (isfinite(pRec->stokes.I+pRec->stokes.U))  k++;
//		if(k == 1)
//			printf("RA %2.8f  DEC %2.8f \n",daydata->records[0].RA,daydata->records[0].DEC);
        }
        else if (num < 0) break;
        else{
		 printf("ERROR: flux file record only read %i fields\n", num);
		exit(1);
	}
    }

    daydata->numRecords = k;
    daydata->RAmin = RAmin;
    daydata->RAmax = RAmax;
    return k;
}

//----------- read in data from input flux files
//dataset must be allocated of size ndays
//dataset[].records will be malloced new memory
int fluxwappdata_readchan(FluxWappData * wappdata, int chan, int id)
{
    int  m,j;
    int count;
    FILE *infile;
    char beamno[6];
    count = 0;
    for (m=0; m<wappdata->numDays; m++)
    {
        char filename[64+1];
        FluxDayData * daydata = &wappdata->daydata[m];

	if(!strcmp(wappdata->wapp,"multibeam"))
	{
		j = m%7;
		if(j == 6)
		{
			daydata->numRecords = 0;
			//break;
		}
		sprintf(beamno,"beam%d",j);
		if(id == CLEAN)
	        	//sprintf(filename, "%s/%s/balance%04i.dat", daydata->mjd, beamno, chan);
	        	sprintf(filename, "/n/fox/processed/FIELD1/run8/band0/%s/%s/balance%04i.dat", daydata->mjd, beamno, chan);
			//sprintf(filename, "%s/%s/subregion%04i_01.dat", daydata->mjd, beamno, chan);
    	 	if(id == BASKETWEAVE)
       			sprintf(filename, "%s/%s/fluxtime%04i.dat", daydata->mjd, beamno, chan);
    	 	if(id == CORR)
       			sprintf(filename, "%s/%s/corr%04i.dat", daydata->mjd, beamno, chan);
 	//       	sprintf(filename, "%s/%s/clean%03i.dat", daydata->mjd, beamno, chan); //SSG hack for cleanmain
	}
	else
	{
		//printf("%s\n",wappdata->wapp);
		j = atoi(&wappdata->wapp[4]);
		if(id == CLEAN)
			sprintf(filename, "%s/%s/balance%04i.dat", daydata->mjd, wappdata->wapp, chan);
	        if(id == BASKETWEAVE)
		    	sprintf(filename, "%s/%s/fluxtime%04i.dat", daydata->mjd, wappdata->wapp, chan);
	        if(id == CORR)
		    	sprintf(filename, "%s/%s/corr%04i.dat", daydata->mjd, wappdata->wapp, chan);
	}


    	infile = fopen(filename, "r");
    	if (infile == NULL) {
//    		printf("ERROR: can't open input file %s\n", filename);
	        continue;
    	}
	//else//SSG
	//	printf("DIAG: Opened file %s ",filename);//SSG
//    	printf(".");
	int numread;
	if (j!= 6)
        numread = fluxdaydata_read(daydata, infile,j);
	else
	numread = 0;
	//printf("read %d records.\n",numread);
        //fflush(stdout);
        fclose(infile);
        count++;
    }
//   printf("\n");
    return count;
}

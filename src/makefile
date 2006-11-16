CC = gcc
CFLAGS = \
-Wall -W \
-Wmissing-prototypes \
-Wshadow \
-Wcast-qual -Wcast-align \
-Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common -g \
-D_GNU_SOURCE \


BIN=../bin

all: cal map

scan: scanmain.c
	$(CC) $(CFLAGS) -lm -o scan common.o fluxdata.o jsd/jsd_futil.o scanmain.c

cube2avg: cube2avg.c
	$(CC) $(CFLAGS) -lm -o cube2avg cube2avg.c programs/fitsio.a
	mv cube2avg $(BIN)

pointing: pointingmain.c common.o
	$(CC) $(CFLAGS) -lm -o pointing common.o pointingmain.c
	mv pointing $(BIN)

calanalyze: calanalyze.c
	$(CC) $(CFLAGS) -lm -o calanalyze calanalyze.c programs/fitsio.a
	cp calanalyze $(BIN)

fits_smooth: fits_smooth.c
	$(CC) $(CFLAGS) -lm -o fits_smooth fits_smooth.c programs/fitsio.a
	cp fits_smooth $(BIN)

spec2fits: spec2fits.c common.o smooth.o spec.o rfi.o markdata.o
	$(CC) $(CFLAGS) -lm -o spec2fits common.o smooth.o spec.o rfi.o markdata.o spec2fits.c programs/fitsio.a
	cp spec2fits $(BIN)

spectra: spectramain.o common.o spec.o
	$(CC) $(CFLAGS) -lm -o spectra common.o spec.o spectramain.o
	mv spectra $(BIN)

calhack: calhack.c rfi.o stokes.o calibrate.o common.o smooth.o scan.o spec.o
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o calhack rfi.o calibrate.o common.o smooth.o calhack.c spec.o
	mv calhack $(BIN)

cal: calmain.o rfi.o stokes.o calibrate.o common.o smooth.o spec.o markdata.o
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o cal rfi.o stokes.o calibrate.o common.o smooth.o calmain.o spec.o markdata.o programs/fitsio.a jsd/jsd_fit.o jsd/nrutil.o jsd/nrfit.o
	mv cal $(BIN)

spechack : 
	$(CC) -lm -o spechack spechack.c common.o

map: mapmain.o common.o jsd/jsd_futil.o jsd/jsd_util.o balance.o grid.o map.o scan.o fluxdata.o jsd/jsd_fit.o
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o map common.o jsd/jsd_futil.o jsd/jsd_util.o jsd/jsd_fit.o jsd/nrutil.o jsd/nrfit.o balance.o grid.o map.o scan.o fluxdata.o mapmain.o programs/fitsio.a
	mv map $(BIN)

fitscat: fitscat.c programs/fitsio.a
	$(CC) $(CFLAFGS) -lm -o fitscat fitscat.c programs/fitsio.a
	mv fitscat $(BIN)


spec.o : spec.h spec.c
scan.o : scan.h scan.c
smooth.o : smooth.h smooth.c
fluxdata.o : fluxdata.h fluxdata.c
map.o : map.h map.c
grid.o : grid.c grid.h
balance.o : balance.c balance.h
common.o : common.h common.c
rfi.o : common.h rfi.h rfi.c
stokes.o : common.h stokes.h stokes.c
calibrate.o : common.h calibrate.h calibrate.c
calmain.o : common.h rfi.h cal.h calmain.c
mapmain.o : common.h cal.h mapmain.c
spectramain.o : common.h spectramain.c
markdata.o : markdata.c markdata.h
poly_harm_fit.o : poly_harm_fit.f
	g77 -c poly_harm_fit.f
fetch_from_best_fit_model.o: fetch_from_best_fit_model.f
	g77 -c fetch_from_best_fit_model.f



polyeval : polyeval.c
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o polyeval polyeval.c

.PHONY : clean
clean : empirical_clean
	rm *.o


empirical: empirical_determine empirical_apply empirical_intensity

empirical_determine: empirical_determine.c
	gcc -g -lm -o empirical_determine empirical_determine.c programs/fitsio.a

empirical_apply: empirical_apply.c
	gcc -g -lm -o empirical_apply empirical_apply.c

empirical_intensity: empirical_intensity.c
	gcc -g -lm -o empirical_intensity empirical_intensity.c

empirical_clean:
	rm -f empirical_determine empirical_apply empirical_intensity

tcal: tcal.c
	gcc -g -lm -lgsl -lgslcblas -o tcal tcal.c jsd/jsd_fit.o

CC = gcc
CFLAGS = \
-Wall -W \
-Wmissing-prototypes \
-Wshadow \
-Wcast-qual -Wcast-align \
-Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common -g \
-D_GNU_SOURCE \

DIRS = jsd programs leakage

BIN=../bin/

all:	cal map cimafits2spec new_spec2fits spec2fits

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

new_spec2fits: new_spec2fits.c common.o smooth.o new_spec.o
	$(CC) $(CFLAGS) -lm -o new_spec2fits common.o new_spec.o new_spec2fits.c programs/fitsio.a
	cp new_spec2fits $(BIN)
new_spec.o : new_spec.h new_spec.c
spectra: spectramain.o common.o spec.o
	$(CC) $(CFLAGS) -lm -o spectra common.o spec.o spectramain.o
	mv spectra $(BIN)

calhack: calhack.c rfi.o stokes.o calibrate.o common.o smooth.o scan.o spec.o
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o calhack rfi.o calibrate.o common.o smooth.o calhack.c spec.o
	mv calhack $(BIN)

cal: calmain.o rfi.o stokes.o calibrate.o common.o smooth.o spec.o markdata.o stats.o smoothing.o
	$(CC) $(CFLAGS) -lm -o cal rfi.o stokes.o calibrate.o common.o smooth.o smoothing.o calmain.o spec.o stats.o markdata.o programs/fitsio.a jsd/jsd_fit.o jsd/nrutil.o jsd/nrfit.o
	mv cal $(BIN)

spechack : 
	$(CC) -lm -o spechack spechack.c common.o

map: mapmain.o common.o jsd/jsd_futil.o jsd/jsd_util.o balance3.o grid.o map.o scan.o fluxdata.o jsd/jsd_fit.o decdependence.o fieldflat.o stats.o
	$(CC) $(CFLAGS) -lm -o map common.o jsd/jsd_futil.o jsd/jsd_util.o jsd/jsd_fit.o jsd/nrutil.o jsd/nrfit.o balance3.o decdependence.o fieldflat.o grid.o map.o scan.o stats.o fluxdata.o mapmain.o programs/fitsio.a
	#$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o map common.o jsd/jsd_futil.o jsd/jsd_util.o jsd/jsd_fit.o jsd/nrutil.o jsd/nrfit.o balance3.o decdependence.o fieldflat.o grid.o map.o scan.o stats.o fluxdata.o mapmain.o programs/fitsio.a
	mv map $(BIN)

fitscat: fitscat.c programs/fitsio.a
	$(CC) $(CFLAFGS) -lm -o fitscat fitscat.c programs/fitsio.a
	mv fitscat $(BIN)
#SSG
cimafits2spec: cimafits2spec.c cimafits2.h
	$(CC) $(CFLAGS) -lm -o cimafits2spec cimafits2spec.c
	cp cimafits2spec $(BIN)
beamgainratios: beamgainratios.c
	$(CC) $(CFLAGS) -lm -o beamgainratios beamgainratios.c programs/fitsio.a
	mv beamgainratios $(BIN)
beamgains: beamgains.c
	$(CC) $(CFLAGS) -lm -o beamgains beamgains.c programs/fitsio.a
	mv beamgains $(BIN)

beammodels: beammodels.c beammodels.o common.o
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o beammodels common.o beammodels.c beammodelsmain.c programs/fitsio.a jsd/jsd_futil.o
	mv beammodels $(BIN)

cleanmap: cleanmain.c 
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o cleanmap common.o beammodels.o fluxdata.o map.o grid.o cleanmain.c programs/fitsio.a jsd/jsd_futil.o
	mv cleanmap $(BIN)

beammodels.o : beammodels.h beammodels.c
beammodelsmain.o : beammodelsmain.c
#SSG
spec.o : spec.h spec.c
scan.o : scan.h scan.c
smooth.o : smooth.h smooth.c
fluxdata.o : fluxdata.h fluxdata.c
map.o : map.h map.c
grid.o : grid.c grid.h
balance.o : balance.c balance.h
balance2.o : balance2.c balance.h
balance3.o : balance3.c balance.h
common.o : common.h common.c
rfi.o : common.h rfi.h rfi.c
stokes.o : common.h stokes.h stokes.c
calibrate.o : common.h calibrate.h calibrate.c
calmain.o : common.h rfi.h cal.h calmain.c
mapmain.o : common.h cal.h mapmain.c
spectramain.o : common.h spectramain.c
markdata.o : markdata.c markdata.h
decdependence.o : decdependence.c decdependence.h
fieldflat.o : fieldflat.c fieldflat.h
smoothing.o : smoothing.c smoothing.h
stats.o : stats.c stats.h

.PHONY : clean
clean : empirical_clean
	rm *.o 
	for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


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

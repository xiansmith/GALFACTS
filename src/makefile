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

BIN=bin/

all:	cal map cimafits2spec new_spec2fits spec2fits

scan: scanmain.c
	$(CC) $(CFLAGS) -lm -o scan common.o fluxdata.o jsd/jsd_futil.o scanmain.c

cube2avg: cube2avg.c
	$(CC) $(CFLAGS) -lm -o cube2avg cube2avg.c programs/fitsLib.a
	mv cube2avg $(BIN)

quickcube: quickcube.c common.o jsd/jsd_futil.o
	$(CC) $(CFLAGS) -lm -o quickcube quickcube.c common.o jsd/jsd_futil.o programs/fitsLib.a
	mv quickcube $(BIN)

pointing: pointingmain.c common.o
	$(CC) $(CFLAGS) -lm -o pointing common.o pointingmain.c
	mv pointing $(BIN)

vvsi: vvsi.c
	$(CC) $(CFLAGS) -lm -o vvsi vvsi.c programs/fitsLib.a
	mv vvsi $(BIN)

poster: poster.c
	$(CC) $(CFLAGS) -lm -o poster poster.c programs/fitsLib.a jsd/jsd_futil.o
	mv poster $(BIN)

posterbins: posterbins.c jsd/jsd_futil.o jsd/jsd_fit.o stats.o
	$(CC) $(CFLAGS) -lm -o posterbins posterbins.c programs/fitsLib.a jsd/jsd_fit.o jsd/jsd_futil.o stats.o jsd/nrfit.o jsd/nrutil.o
	mv posterbins $(BIN)

calanalyze: calanalyze.c
	$(CC) $(CFLAGS) -lm -o calanalyze calanalyze.c programs/fitsLib.a
	cp calanalyze $(BIN)

fits_smooth: fits_smooth.c
	$(CC) $(CFLAGS) -lm -o fits_smooth fits_smooth.c programs/fitsLib.a
	cp fits_smooth $(BIN)

spIndex: spIndex.c
	$(CC) $(CFLAGS) -lm -o spIndex spIndex.c programs/fitsLib.a jsd/jsd_futil.o jsd/jsd_fit.o jsd/nrfit.o jsd/nrutil.o
	cp spIndex $(BIN)

epsilonPhi: epsilonPhi.c
	$(CC) $(CFLAGS) -lm -o epsilonPhi epsilonPhi.c  denoising.c programs/fitsLib.a
	cp epsilonPhi $(BIN)

spec2fits: spec2fits.c common.o smooth.o spec.o rfi.o markdata.o
	$(CC) $(CFLAGS) -lm -o spec2fits common.o smooth.o spec.o rfi.o markdata.o spec2fits.c programs/fitsLib.a
	cp spec2fits $(BIN)

spectral_analysis: spectral_analysis.c stats.o fluxdata.o common.o jsd/jsd_futil.o
	$(CC) $(CFLAGS) -lm -o spectral_analysis stats.o fluxdata.o common.o jsd/jsd_futil.o spectral_analysis.c
	cp spectral_analysis $(BIN)

dec_polyfits: dec_polyfits.c jsd/jsd_fit.o
		$(CC) $(CFLAGS) -lm -o dec_polyfits jsd/jsd_fit.o jsd/nrfit.o jsd/nrutil.o   dec_polyfits.c
		cp dec_polyfits $(BIN)

spectra: spectramain.o common.o spec.o
	$(CC) $(CFLAGS) -lm -o spectra common.o spec.o spectramain.o
	mv spectra $(BIN)

calhack: calhack.c rfi.o stokes.o calibrate.o common.o smooth.o scan.o spec.o
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o calhack rfi.o calibrate.o common.o smooth.o calhack.c spec.o
	mv calhack $(BIN)

cal: calmain.o rfi.o stokes.o calibrate.o common.o smooth.o spec.o markdata.o stats.o smoothing.o
	$(CC) $(CFLAGS) -lm -o cal rfi.o stokes.o calibrate.o common.o smooth.o smoothing.o calmain.o spec.o stats.o markdata.o jsd/jsd_futil.o jsd/jsd_util.o jsd/jsd_fit.o  jsd/nrfit.o jsd/nrutil.o programs/fitsLib.a
	mv cal $(BIN)

cal_sdfits:
	$(CC) $(CFLAGS) calmain_sdfits.c -o cal_sdfits /u/sguram/cfitsio/libcfitsio.a -lm rfi.o stokes.o calibrate.o common.o smooth.o smoothing.o spec.o stats.o markdata.o jsd/jsd_fit.o  jsd/nrfit.o jsd/nrutil.o jsd/jsd_futil.o  # programs/fitsLib.a
	mv cal_sdfits $(BIN)
spechack : 
	$(CC) -lm -o spechack spechack.c common.o

map: mapmain.o common.o jsd/jsd_futil.o jsd/jsd_util.o balance.o grid.o map.o scan.o fluxdata.o jsd/jsd_fit.o decdependence.o fieldflat.o stats.o leakagefit.o
	$(CC) $(CFLAGS) -lm -o map common.o jsd/jsd_futil.o jsd/jsd_util.o jsd/jsd_fit.o jsd/nrutil.o jsd/nrfit.o balance.o decdependence.o fieldflat.o grid.o map.o scan.o stats.o fluxdata.o mapmain.o programs/fitsLib.a leakagefit.o
	mv map $(BIN)

fitscat: fitscat.c programs/fitsio.a
	$(CC) $(CFLAFGS) -lm -o fitscat fitscat.c programs/fitsLib.a
	mv fitscat $(BIN)

dtm2spec: dtm2spec.c dtm2spec.h jsd/jsd_futil.o
	$(CC) $(CFLAGS) -I. -lm -o dtm2spec dtm2spec.c jsd/jsd_futil.o
	cp dtm2spec $(BIN)

caldeflection: caldeflection.c dtm2spec.h jsd/jsd_futil.o
	$(CC) $(CFLAGS) -I. -lm -o caldeflection caldeflection.c jsd/jsd_futil.o
	mv caldeflection $(BIN)

calamptest: calamptest.c dtm2spec.h jsd/jsd_futil.o
	$(CC) $(CFLAGS) -I. -lm -o calamptest calamptest.c jsd/jsd_futil.o
	cp calamptest $(BIN)

dfrq2spec: dfrq2spec.c dfrq2spec.h jsd/jsd_futil.o
	$(CC) $(CFLAGS) -I. -lm -o dfrq2spec dfrq2spec.c jsd/jsd_futil.o
	cp dfrq2spec $(BIN)

read_spec: read_spec.c read_spec.h
	$(CC) $(CFLAGS) -I. -lm -o read_spec read_spec.c
	cp read_spec $(BIN)
	
cimafits2spec: cimafits2spec.c cimafits2.h
	$(CC) $(CFLAGS) -lm -o cimafits2spec cimafits2spec.c
	cp cimafits2spec $(BIN)

beamgainratios: beamgainratios.c
	$(CC) $(CFLAGS) -lm -o beamgainratios beamgainratios.c programs/fitsLib.a
	mv beamgainratios $(BIN)

beamgains: beamgains.c
	$(CC) $(CFLAGS) -lm -o beamgains beamgains.c programs/fitsLib.a
	mv beamgains $(BIN)

beammodels: beammodels.c beammodels.o common.o
	$(CC) $(CFLAGS) -lm  -o beammodels common.o beammodels.c beammodelsmain.c programs/fitsLib.a jsd/jsd_futil.o
	mv beammodels $(BIN)

beamtest: beamtest.c
	$(CC) $(CFLAGS) -lm  -o beamtest beamtest.c jsd/jsd_futil.o
	mv beamtest $(BIN)

rmtest: rmtest.c 
	$(CC) $(CFLAGS) -lm  -o rmtest rmtest.c programs/fitsLib.a 
	mv rmtest $(BIN)

cleanmap: cleanmain.c 
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas -o cleanmap common.o beammodels.o fluxdata.o map.o grid.o cleanmain.c programs/fitsLib.a jsd/jsd_futil.o
	mv cleanmap $(BIN)

quleakage: QUleakage.c 
	$(CC) $(CFLAGS) QUleakage.c -o quleakage common.o jsd/jsd_futil.o
	mv quleakage $(BIN)

leakage_pattern: leakage_pattern.c 
	$(CC) $(CFLAGS) -lm leakage_pattern.c -o leakage_pattern common.o jsd/jsd_futil.o
	mv leakage_pattern $(BIN)

corrmain: corrmain.c common.o beammodels.o fluxdata.o map.o grid.o corrmain.c programs/fitsLib.a jsd/jsd_futil.o cfit/chebyshev.c cfit/pfit.c
	$(CC) $(CFLAGS) -lm -O2 -ffast-math -o corrmain common.o beammodels.o fluxdata.o map.o grid.o corrmain.c programs/fitsLib.a jsd/jsd_futil.o cfit/chebyshev.c cfit/pfit.c
	mv corrmain $(BIN)

corrtest: corrtest.c common.o beammodels.o fluxdata.o map.o grid.o corrtest.c programs/fitsLib.a jsd/jsd_futil.o cfit/chebyshev.c cfit/pfit.c leakagefit.c jsd/jsd_fit.o jsd/nrfit.o jsd/nrutil.o
	$(CC) $(CFLAGS) -lm -O2 -ffast-math -o corrtest common.o beammodels.o fluxdata.o map.o grid.o corrtest.c programs/fitsLib.a jsd/jsd_futil.o cfit/chebyshev.c cfit/pfit.c leakagefit.c jsd/jsd_fit.o jsd/nrfit.o jsd/nrutil.o
	mv corrtest $(BIN)

beammodels.o : beammodels.h beammodels.c
beammodelsmain.o : beammodelsmain.c
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
leakagefit.o: leakagefit.c
.PHONY : clean
clean : empirical_clean
	rm *.o 
	for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


empirical: empirical_determine empirical_apply empirical_intensity

empirical_determine: empirical_determine.c
	gcc -g -lm -o empirical_determine empirical_determine.c programs/fitsLib.a

empirical_apply: empirical_apply.c
	gcc -g -lm -o empirical_apply empirical_apply.c

empirical_intensity: empirical_intensity.c
	gcc -g -lm -o empirical_intensity empirical_intensity.c

empirical_clean:
	rm -f empirical_determine empirical_apply empirical_intensity

tcal: tcal.c
	gcc -g -lm  -o tcal tcal.c jsd/jsd_fit.o programs/fitsLib.a jsd/nrfit.o jsd/nrutil.o
	mv tcal $(BIN)
#	gcc -g -lm -lgsl -lgslcblas -o tcal tcal.c jsd/jsd_fit.o programs/fitsLib.a jsd/nrfit.o jsd/nrutil.o
get_subregion:
	gcc -o get_subregion get_subregion.c jsd/jsd_futil.c common.c
	mv get_subregion bin/

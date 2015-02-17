all: map
#I have only changed this Makefile and not mk
# -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
programs: 
	gcc -O3 -c programs/fitsLib.c -o programs/fitsLib.o
	ar rcs programs/fitsLib.a programs/fitsLib.o
map: mapmain.c balance.c common.c scan.c grid.c map.c fluxdata.c fix_pointing.c decdep.c chebyshev.c abscal.c correctUV.c corrections.c QUVcorr.c jsd/jsd_futil.c jsd/jsd_util.c programs/io.c programs/fitsio.c programs/fitsLib.a
	gcc  -std=c99 -fshort-enums -fno-common -g -D_GNU_SOURCE \
	-O3 -ffast-math \
	-o map \
	mapmain.c balance.c common.c scan.c grid.c map.c fluxdata.c decdep.c chebyshev.c \
	abscal.c correctUV.c corrections.c QUVcorr.c jsd/jsd_futil.c jsd/jsd_util.c \
	fix_pointing.c programs/io.c programs/fitsio.c programs/fitsLib.a \
	-lm
	chmod 755 map
	mv map /opt/exp_soft/galfacts
programs/fitsLib.a:
	gcc -O3 -c programs/fitsLib.c -o programs/fitsLib.o
	ar rcs programs/fitsLib.a programs/fitsLib.o
clean:
	rm -f map programs/*.o programs/*.a
	echo Clean complete

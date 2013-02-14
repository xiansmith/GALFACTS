all: map
map:
	gcc -O3 -c programs/fitsLib.c -o programs/fitsLib.o
	ar rcs programs/fitsLib.a programs/fitsLib.o
	gcc -fshort-enums -fno-common -g -D_GNU_SOURCE \
	-O3 -ffast-math \
	-o map \
	mapmain.c balance.c common.c scan.c grid.c map.c fluxdata.c decdep.c chebyshev.c \
	jsd/jsd_futil.c jsd/jsd_util.c \
	programs/io.c programs/fitsio.c programs/fitsLib.a \
	-lm

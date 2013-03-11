all: map
#I have only changed this Makefile and not mk
map:
	gcc -O3 -c programs/fitsLib.c -o programs/fitsLib.o
	ar rcs programs/fitsLib.a programs/fitsLib.o
	gcc -fshort-enums -fno-common -g -D_GNU_SOURCE \
	-O3 -ffast-math \
	-o map \
	mapmain.c balance.c common.c scan.c grid.c map.c fluxdata.c decdep.c chebyshev.c \
	correctUV.c jsd/jsd_futil.c jsd/jsd_util.c \
	programs/io.c programs/fitsio.c programs/fitsLib.a \
	-lm
clean:
	rm -f map programs/*.o programs/*.a
	echo Clean complete

CC = gcc
CFLAGS = \
-Wall -W \
-Wmissing-prototypes \
-Wshadow \
-Wcast-qual -Wcast-align \
-Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common -g \
-D_GNU_SOURCE \


objects = io.o misc_math.o fitsio.o


all : $(objects) fitsio.a

fitsio.a :
	ar r fitsio.a $(objects)
	rm -rf $(objects)

io.o : io.h io.c
misc_math.o : misc_math.h misc_math.c mathdefs.h
fitsio.o : fitsio.h fitsio.c io.h


.PHONY clean :
clean :
	rm -rf $(objects) fitsio.a


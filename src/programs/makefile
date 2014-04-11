CC = gcc
CFLAGS = \
-Wall -W \
-Wmissing-prototypes \
-Wshadow \
-Wcast-qual -Wcast-align \
-Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common -g \
-D_GNU_SOURCE \


objects = io.o misc_math.o fitsLib.o


all : $(objects) fitsLib.a

fitsLib.a :
	ar r fitsLib.a $(objects)
	rm -rf $(objects)

io.o : io.h io.c
misc_math.o : misc_math.h misc_math.c mathdefs.h
fitLib.o : fitLib.h fitLib.c io.h


.PHONY clean :
clean :
	rm -rf $(objects) fitsio.a


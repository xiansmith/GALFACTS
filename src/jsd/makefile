HOME=/home/jeff
CC = gcc
CFLAGS = \
-Wall -W \
-Wmissing-prototypes \
-Wcast-qual -Wcast-align \
-Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common -g \
-D_GNU_SOURCE \


all : jsd_futil.o jsd_util.o jsd_fit.o

jsd_futil.o : jsd_futil.c jsd_futil.h
jsd_util.o : jsd_util.c jsd_util.h
jsd_fit.o : jsd_fit.c jsd_fit.h
nrfit.o : nrfit.h nrfit.c
nrutil.o : nrutil.h nrutil.c
jsd_fit_test : jsd_fit_test.c jsd_fit.o jsd_futil.o nrfit.o nrutil.o
	$(CC) $(CFLAGS) -o jsd_fit_test jsd_fit_test.c jsd_fit.o jsd_futil.o -lm -lgsl -lgslcblas nrutil.o nrfit.o

clean :
	rm *.o 

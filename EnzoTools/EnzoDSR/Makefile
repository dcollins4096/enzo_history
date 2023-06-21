SITE=THUNDERHEAD
SYSTEM=$(shell uname -s)


ifeq ($(SITE),THUNDERHEAD)
HDF5_INSTALL = /sw
HDF4_INSTALL = 
endif


ifeq ($(SYSTEM),Darwin)
CC         = g++
CFLAGS     = -O2
LIB        = -L/usr/local/lib -lz
LIBSHDF = $(HDF5_INSTALL)/lib/libhdf5.a
endif


INCLUDE   = -I$(HDF5_INSTALL)/include

enzodsr:  EnzoDSR.C
	$(CC) $(CFLAGS) -o $@ EnzoDSR.C $(INCLUDE) $(LIBSHDF) $(LIB)

all:
	make enzodsr

clean: 
	rm -f *.o *~ enzodsr

.SUFFIXES:.o.c

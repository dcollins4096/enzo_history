SITE=SDSC_SP4
SYSTEM=$(shell uname -s)

 # yes or no (this is for cluster extraction - make two projections)
CLU_SPL := no

ifeq ($(SITE),SDSC_SP4)
#SDSC
#HDF4_INSTALL = /users02/ucsd/ux453739/hdf4-64
#HDF5_INSTALL = /users02/ucsd/ux453739/hdf5-64
endif

ifeq ($(SITE),FNORD)
HDF5_INSTALL = /sw
HDF4_INSTALL = 
endif

ifeq ($(SITE),SDSC_SP3)
#SDSC
HDF4_INSTALL = /usr/local/apps/apps64/hdf4
HDF5_INSTALL = /usr/local/apps/apps64/hdf5
endif

ifeq ($(SYSTEM),AIX)
CC          = xlC_r
CFLAGS      = -q64 -qtune=pwr4 -qarch=pwr4 -O3 -qstrict
LIB         = -lm
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a $(ZLIB_INSTALL)/lib/libz.a $(SZIP_INSTALL)/lib/libsz.a
endif

ifeq ($(SYSTEM),Darwin)
CC         = g++
CFLAGS	   = -O2
LIB	   = -L/usr/local/lib -lz -L/sw/lib -lsz
LIBSHDF = $(HDF5_INSTALL)/lib/libhdf5.a
endif

ifeq ($(CLU_SPL),yes)
  CLUST_SP := USE_CLUSTER_SPLIT
else
  CLUST_SP := NO_USE_CLUSTER_SPLIT
endif

INCLUDE   = -I$(HDF5_INSTALL)/include

collage:  main.o UserDefines.o Utilities.o ProjectionRoutines.o 
	$(CC) $(CFLAGS) -o $@ \
	main.o UserDefines.o Utilities.o ProjectionRoutines.o \
	$(INCLUDE) $(LIBSHDF) $(LIB)


main.o:
	$(CC) $(CFLAGS) $(INCLUDE) -D $(CLUST_SP) -c main.C   

Utilities.o:
	$(CC) $(CFLAGS) $(INCLUDE) -D $(CLUST_SP) -c Utilities.C   

UserDefines.o:
	$(CC) $(CFLAGS) $(INCLUDE) -D $(CLUST_SP) -c UserDefines.C   

ProjectionRoutines.o:
	$(CC) $(CFLAGS) $(INCLUDE) -D $(CLUST_SP) -c ProjectionRoutines.C   

clean: 
	touch temp.o
	rm -f *.o *~ collage

.SUFFIXES:.o.c

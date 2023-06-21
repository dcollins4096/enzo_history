#
#  Makefile for example enzo dataset reader.  The user will
#  have to download and install the HDF 5 library from 
#  http://hdf.ncsa.uiuc.edu/HDF5.  Make sure to download
#  the 1.4.5-post2 version of the library, NOT the most 
#  current version.
#
#  The code can be compiled by typing "make" or "gmake".
#
#  Email questions to Brian O'Shea (bwoshea@cosmos.ucsd.edu)
#


# User must modify this to reflect where their HDF 5 libary 
# has been installed, the desired C++ compiler and the
# location of libz.a (the szip library, required for HDF 5).
HDF5_INSTALL = /users02/ucsd/ux453739/hdf5-64
HDF4_INSTALL = /users02/ucsd/ux453739/hdf4-64
CC         = xlC_r
CFLAGS	   = -q64 -qtune=pwr4 -qarch=pwr4 -O3 -qstrict
#LIB	   = -L/usr/lib -lz


LIBSHDF = $(HDF5_INSTALL)/lib/libhdf5.a $(HDF4_INSTALL)/lib/libz.a
INCLUDE   = -I$(HDF5_INSTALL)/include

hdfreader: enzo_dataset_reader.o
	$(CC) $(CFLAGS) -o $@ \
	enzo_dataset_reader.o \
	$(INCLUDE) $(LIBSHDF) $(LIB)

enzo_dataset_reader.o :
	$(CC) $(CFLAGS) -c enzo_dataset_reader.C $(INCLUDE) $(LIBSHDF) $(LIB)

clean: 
	touch temp.o
	rm -f *.o *~ hdfreader

.SUFFIXES:.o.c

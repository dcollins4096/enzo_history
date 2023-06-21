SITE=SDSC_SP4
SYSTEM=$(shell uname -s)

ifeq ($(SITE),FNORD)
HDF5_INSTALL = /usr/local
HDF4_INSTALL = 
endif

ifeq ($(SITE),SDSC_SP4)
#SDSC
HDF4_INSTALL = /users02/ucsd/ux453739/hdf4-64
HDF5_INSTALL = /users02/ucsd/ux453739/hdf5-64
endif


ifeq ($(SITE),NCSA_O2K)
#NCSA
HDF4_INSTALL = /afs/ncsa/packages/hdf/4.1r5-irix64-64
HDF5_INSTALL = /afs/ncsa/packages/hdf5/5-1.4.5-irix64
endif

ifeq ($(SITE),SDSC_SP3)
#SDSC
HDF4_INSTALL = /usr/local/apps/apps64/hdf4
HDF5_INSTALL = /usr/local/apps/apps64/hdf5
#HDF5_INSTALL = /projects/paci/SciComp/harkness/hdf5/hdf5-1.4.2-alpha
#HDF5_INSTALL = /projects/paci/SciComp/harkness/hdf5/hdf5-1.4.2-alpha-pp
endif

ifeq ($(SITE),SDSC_SUN)
#SDSC
HDF5_INSTALL = /projects/paci/SciComp/harkness/hdf5/5-1.4.3-sol-64
HDF4_INSTALL = /projects/paci/SciComp/harkness/hdf4/4.1r5-solaris8
endif

ifeq ($(SYSTEM),IRIX64)
CC          = CC -mips4
CFLAGS      = -prelink -DMPI_NO_CPPBIND -O2
#LIB         = -lm -lC -lmpi
LIB         = -lm -lC 
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a $(HDF4_INSTALL)/lib/libz.a
#CFLAGS      = -prelink -TARG:platform=ip27 -OPT:Olimit=0,alias=typed -r10000 -DMPI_NO_CPPBIND
#LIB         = -lm -lC -lmalloc_ss -lmpi
endif

ifeq ($(SYSTEM),AIX)
CC          = xlC_r
CFLAGS      = -q64 -qarch=pwr4 -qtune=pwr4 # -bnoquiet
LIB         = -lm
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a $(HDF4_INSTALL)/lib/libz.a
endif

ifeq ($(SYSTEM),SunOS)
CC          = CC -xarch=v9
CFLAGS      = -O2
LIB         = -L/usr/local/lib -lsunmath -lnsl -lm -lsocket -lz
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a
endif

ifeq ($(SYSTEM),Darwin)
CC         = g++
CFLAGS     = -O2
LIB        = -L/usr/local/lib -lz
LIBSHDF = $(HDF5_INSTALL)/lib/libhdf5.a
endif


INCLUDE   = -I$(HDF5_INSTALL)/include

denpeak:  DenPeak.C
	$(CC) $(CFLAGS) -o $@ DenPeak.C $(INCLUDE) $(LIBSHDF) $(LIB)

clean: 
	rm -f *.o *~ denpeak 

.SUFFIXES:.o.c

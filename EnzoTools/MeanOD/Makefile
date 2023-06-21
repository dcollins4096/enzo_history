SITE= SDSC_SP4
SYSTEM=$(shell uname -s)

ifeq ($(SITE),FNORD)
HDF5_INSTALL = /usr/local
HDF4_INSTALL = 
endif

ifeq ($(SITE),SDSC_SP4)
#SDSC
#HDF4_INSTALL = /users02/ucsd/ux453739/hdf4-64
#HDF5_INSTALL = /users02/ucsd/ux453739/hdf5-64

#HDF4_INSTALL = /users02/ucsd/ux453739/hdf4-64
#HDF5_INSTALL = /users02/ucsd/ux453739/HDF5/5-1.6.3-aix5.2-64
endif


ifeq ($(SITE),NCSA_SP4)
#SDSC
HDF4_INSTALL = /u/ac/bwoshea/hdf4/
HDF5_INSTALL = /u/ac/bwoshea/hdf5/
endif


ifeq ($(SITE),NCSA_O2K)
#NCSA
HDF4_INSTALL = /afs/ncsa/packages/hdf/4.1r5-irix64-64
HDF5_INSTALL = /afs/ncsa/packages/hdf5/IRIX64_MIPS4_6.5/
endif

ifeq ($(SITE),SDSC_SP3)
#SDSC
HDF4_INSTALL = /usr/local/apps/apps64/hdf4
HDF5_INSTALL = /usr/local/apps/apps64/hdf5
endif

ifeq ($(SITE),SDSC_SUN)
#SDSC
HDF5_INSTALL = /projects/paci/SciComp/harkness/hdf5/5-1.4.3-sol-64
HDF4_INSTALL = /projects/paci/SciComp/harkness/hdf4/4.1r5-solaris8
endif



ifeq ($(SYSTEM),IRIX64)
CC          = CC -mips4
#CFLAGS      = -prelink -DMPI_NO_CPPBIND -O2 -g
CFLAGS      = -prelink -DMPI_NO_CPPBIND -g
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
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a $(ZLIB_INSTALL)/lib/libz.a $(SZIP_INSTALL)/lib/libsz.a
endif

ifeq ($(SYSTEM),SunOS)
CC          = CC -xarch=v9
CFLAGS      = -O2
LIB         = -L/usr/local/lib -lsunmath -lnsl -lm -lsocket -lz
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a
endif

ifeq ($(SYSTEM),Linux)
  ifeq ($(UNAMEM),ia64)
    ENZO_ARCH = tg64_sdsc_v7
    ENZO_MESS = mpich
  endif
  ifeq ($(ENZO_INT),LARGE_INTS)
    ifeq ($(ENZO_XPR),r8)
      FORTPR = -r8 -i8
    else
      FORTPR = -r4 -i8
    endif
  else
    ifeq ($(ENZO_XPR),r8)
      FORTPR = -r8 -i4
    else
      FORTPR = -r4 -i4
    endif
  endif
endif

ifeq ($(SYSTEM),Linux)
CC      = ecc
CFLAGS = -O2
LIB = -L$(HDF5_INSTALL)/lib -lhdf5 -L/usr/lib -lz -ljpeg
endif


ifeq ($(SYSTEM),Darwin)
CC         = g++
CFLAGS	   = -O2
LIB	   = -L/usr/local/lib -lz
LIBSHDF = $(HDF5_INSTALL)/lib/libhdf5.a
endif

INCLUDE   = -I$(HDF5_INSTALL)/include

meanod: MeanOD.o 
	$(CC) $(CFLAGS) -o $@ \
	MeanOD.o  \
	$(INCLUDE) $(LIBSHDF) $(LIB)

MeanOD.o:
	$(CC) $(CFLAGS) -c MeanOD.C $(INCLUDE) $(LIBSHDF) $(LIB)

clean: 
	touch temp.o
	rm -f *.o *~ meanod

.SUFFIXES:.o.c

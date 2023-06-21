SITE=SDSC_SP4    # this corresponds to the appropriate site (so the code knows where
                 # the HDF and associated libraries are).
PRECISION = R4   # this can be R4 or R8 depending on whether files have 32 or 64 bit data, respectively

AMRPACK = PACK_AMR  # this can be PACK_AMR or NO_PACK_AMR  

MEKAL = USE_MEKAL  # this can be USE_MEKAL or NO_USE_MEKAL

SYSTEM=$(shell uname -s)

ifeq ($(SITE),SDSC_SP4)
#SDSC
HDF4_INSTALL = /users02/ucsd/ux453739/hdf4-64
HDF5_INSTALL = /users02/ucsd/ux453739/hdf5-64
endif

ifeq ($(SITE),FNORD)
HDF5_INSTALL = /sw
HDF4_INSTALL = /sw
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
CFLAGS      = -prelink -DMPI_NO_CPPBIND -g -D $(PRECISION) -D $(AMRPACK) -D $(MEKAL)
LIB         = -lm -lC 
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a $(HDF4_INSTALL)/lib/libz.a
endif

ifeq ($(SYSTEM),AIX)
CC          = xlC_r
CFLAGS      = -q64 -qtune=pwr4 -qarch=pwr4 -O3 -qstrict -D $(PRECISION) -D $(AMRPACK) -D $(MEKAL)
LIB         = -lm
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a $(HDF4_INSTALL)/lib/libz.a
endif

ifeq ($(SYSTEM),SunOS)
CC          = CC -xarch=v9
CFLAGS      = -O2 -D $(PRECISION) -D $(AMRPACK) -D $(MEKAL)
LIB         = -L/usr/local/lib -lsunmath -lnsl -lm -lsocket -lz
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a
endif

ifeq ($(SYSTEM),Darwin)
CC         = g++
CFLAGS	   = -O2 -D $(PRECISION) -D $(AMRPACK) -D $(MEKAL)
LIB	   = -L/usr/local/lib -lz
LIBSHDF = $(HDF5_INSTALL)/lib/libhdf5.a
endif

INCLUDE   = -I$(HDF5_INSTALL)/include

enzoproj:  EP_Main.o  EP_ParseInputs.o EP_GridStuff.o EP_Misc.o EP_MakeProj.o EP_PartProj.o EP_MEKAL.o
	$(CC) $(CFLAGS) -o $@ \
	EP_Main.o EP_ParseInputs.o EP_GridStuff.o  \
	EP_Misc.o EP_MakeProj.o EP_PartProj.o EP_MEKAL.o \
	$(INCLUDE) $(LIBSHDF) $(LIB)

EP_Main.o:
	$(CC) $(CFLAGS) -c EP_Main.C $(INCLUDE) $(LIBSHDF) $(LIB)

EP_ParseInputs.o:
	$(CC) $(CFLAGS) -c EP_ParseInputs.C $(INCLUDE) $(LIBSHDF) $(LIB)

EP_GridStuff.o:
	$(CC) $(CFLAGS) -c EP_GridStuff.C $(INCLUDE) $(LIBSHDF) $(LIB)

EP_Misc.o:
	$(CC) $(CFLAGS) -c EP_Misc.C $(INCLUDE) $(LIBSHDF) $(LIB)

EP_MEKAL.o:
	$(CC) $(CFLAGS) -c EP_MEKAL.C $(INCLUDE) $(LIBSHDF) $(LIB)

EP_MakeProj.o:
	$(CC) $(CFLAGS) -c EP_MakeProj.C $(INCLUDE) $(LIBSHDF) $(LIB)

EP_PartProj.o:
	$(CC) $(CFLAGS) -c EP_PartProj.C $(INCLUDE) $(LIBSHDF) $(LIB)

clean: 
	touch temp.o
	rm -f *.o *~ enzoproj

.SUFFIXES:.o.c

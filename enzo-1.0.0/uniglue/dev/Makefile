SITE=NCSA_O2K
SYSTEM=$(shell uname -s)

HDF4_INSTALL = ../../libs/HDF4
HDF5_INSTALL = ../../libs/HDF5

ifeq ($(SYSTEM),IRIX64)
CC          = CC -mips4
CFLAGS      = -prelink -DMPI_NO_CPPBIND -O2 -DHDF5_BE -DIO_LOG
LIB         = -lm -lC -lmpi
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a $(HDF4_INSTALL)/lib/libz.a
#CFLAGS      = -prelink -TARG:platform=ip27 -OPT:Olimit=0,alias=typed -r10000 -DMPI_NO_CPPBIND
#LIB         = -lm -lC -lmalloc_ss -lmpi
endif

ifeq ($(SYSTEM),AIX)
CC          = mpCC_r
CFLAGS      = -bmaxdata:0x80000000 -bmaxstack:0x70000000
LIB         = -lm
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a $(HDF4_INSTALL)/lib/libz.a
endif

ifeq ($(SYSTEM),SunOS)
CC          = CC -xarch=v9
CFLAGS      = -O2
LIB         = -L/usr/local/lib -lsunmath -lnsl -lm -lsocket -lz
LIBSHDF     = $(HDF5_INSTALL)/lib/libhdf5.a
endif

INCLUDE   = -I$(HDF5_INSTALL)/include

glop: glop.C
	$(CC) $(CFLAGS) -c glop.C $(INCLUDE) $(LIBSHDF) $(LIB)
	$(CC) $(CFLAGS) -c sub.C $(INCLUDE) $(LIBSHDF) $(LIB)
	$(CC) $(CFLAGS) -o $@ glop.o sub.o $(INCLUDE) $(LIBSHDF) $(LIB)

select: select.C
	$(CC) $(CFLAGS) -o $@ select.C $(INCLUDE) $(LIBSHDF) $(LIB)

dumper: dumper.C
	$(CC) $(CFLAGS) -o $@ dumper.C $(INCLUDE) $(LIBSHDF) $(LIB)

logger: logger.C
	$(CC) $(CFLAGS) -o $@ logger.C $(INCLUDE) $(LIBSHDF) $(LIB)

differ: differ.C
	$(CC) $(CFLAGS) -o $@ differ.C $(INCLUDE) $(LIBSHDF) $(LIB)

clean: 
	rm -f *.o glop select dumper logger differ

.SUFFIXES:.o.c

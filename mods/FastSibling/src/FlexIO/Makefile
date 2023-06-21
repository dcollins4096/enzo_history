# Makefile for the IEEE library
# You *MUST* use GNU gmake for this to work!!!!!
# Go get it from ftp://prep.ai.mit.edu/gnu/make
# The 9000 names of the T3E...
TMPUN   := $(shell uname)
ifeq ($(TMPUN), CYGWIN32/95)
UNAME   = CYGWIN32/95
else
# treat NT the same as 95 (historically, 95 came first...)
ifeq ($(TMPUN), CYGWIN32/NT)
UNAME   = CYGWIN32/95
else
UNAME            := $(shell uname | perl -pe 's/(sn\d\d\d\d|jsimpson)/UNICOS\/mk/')
endif
endif

# F77 = f77 # default fortran compiler

# Various architecture defaults
#Solaris/SunOS
ifeq ($(UNAME), SunOS)
#CC = gcc 
#C++ = g++ 
CC = cc 
C++ = CC -pto
F77 = f77
INCLUDES = -I.
HDFROOT = /afs/ncsa.uiuc.edu/packages/hdf/SunOS_5.4
HDFINCLUDES = -I$(HDFROOT)/include
HDFLIBPATH = $(HDFROOT)/lib
HDFLIBS = -L$(HDFROOT)/lib -lmfhdf -ldf -lz -ljpeg -lrpcsvc -lnsl
LIBS = -L. -lieeeio
CFLAGS = -O $(INCLUDES) -DANSI -DSUN -DSOLARIS 
CCFLAGS = $(CFLAGS) 
FFLAGS = -O $(INCLUDES) -DANSI -DSUN -DSOLARIS
CCFLAGS = $(CFLAGS)
ARCHIVE = ar crs
OBJ_SUFF = .o
C2OBJ = -c
endif

# 64 bit IRIX
ifeq ($(UNAME), IRIX64)
CC = cc
C++ = CC
INCLUDES = -I.
HDFROOT = /afs/ncsa.uiuc.edu/packages/hdf/IRIX_6.4
#HDFROOT = /usr/local
HDFINCLUDES = -I$(HDFROOT)/include
HDFLIBS = -L$(HDFLIBPATH) -lmfhdf -ldf -lz -ljpeg
LIBS = -L. -lieeeio 
HDFLIBPATH = $(HDFROOT)/lib
CFLAGS = -64 -mips4 -LANG:ansi-for-init-scope=ON -IPA -LNO -O3 -ptused $(INCLUDES) -DANSI -DSGI
ifeq ($(SGI_BITS), n32)
HDFROOT = /afs/ncsa.uiuc.edu/packages/ieeeio/IRIXn32_6.4
HDFLIBPATH = $(HDFROOT)/lib
CFLAGS = -mips4 -n32 -LANG:ansi-for-init-scope=ON -O -ptused $(INCLUDES) -DANSI -DSGI
endif
ifeq ($(SGI_BITS), 32)
HDFROOT = /afs/ncsa.uiuc.edu/packages/ieeeio/IRIX_5.3
HDFLIBPATH = $(HDFROOT)/lib
CFLAGS = -32 -O -ptused $(INCLUDES) -DANSI -DSGI
endif
CCFLAGS= $(CFLAGS)
FFLAGS = $(CFLAGS) 
ARCHIVE = /bin/ar crs 
OBJ_SUFF = .o
C2OBJ = -c
endif

# 32/n32 bit Irix
ifeq ($(UNAME), IRIX)
CC = cc
C++ = CC
INCLUDES = -I.
HDFROOT = /afs/ncsa.uiuc.edu/packages/hdf/IRIX_5.3
HDFINCLUDES = -I$(HDFROOT)/include
LIBS = -L. -lieeeio
ifeq ($(SGI_BITS), n32)
HDFROOT = /afs/ncsa.uiuc.edu/packages/ieeeio/IRIXn32_6.2
HDFLIBPATH = $(HDFROOT)/lib
CFLAGS = -n32 -mips3 -LANG:ansi-for-init-scope=ON -ptused -O $(INCLUDES) -DANSI
FFLAGS = $(CFLAGS) 
CCFLAGS = $(CFLAGS)
else
HDFLIBPATH = $(HDFROOT)/lib
CFLAGS = -32 -ptused -O $(INCLUDES) -DANSI
FFLAGS = $(CFLAGS) 
CCFLAGS = $(CFLAGS)
endif
HDFLIBS = -L$(HDFLIBPATH) -lmfhdf -ldf -ljpeg -lz
ARCHIVE = /bin/ar crs 
OBJ_SUFF = .o
C2OBJ = -c
endif

# HP-UX / Convex SPP-2000
ifeq ($(UNAME), HP-UX)
CC = cc
C++ = CC
F77 = f77
INCLUDES = -I.
HDFROOT = /afs/ncsa.uiuc.edu/packages/hdf/SPP2000
HDFINCLUDES = -I$(HDFROOT)/include
HDFLIBPATH = $(HDFROOT)/lib
HDFLIBS = -L$(HDFROOT)/lib -lmfhdf -ldf -ljpeg -lz
LIBS = -L. -lieeeio
# Kludge for HP S/X-class
CFLAGS = -O $(INCLUDES) -DANSI -DHP +a1
CCFLAGS = $(CFLAGS)
FFLAGS = -O $(INCLUDES) -DANSI -DHP
CCFLAGS = $(CFLAGS)
ARCHIVE = ar crs
OBJ_SUFF = .o
C2OBJ = -c
endif

# DEC Alpha
ifeq ($(UNAME), OSF1)
CC = gcc
C++ = g++
INCLUDES = -I.
HDFROOT = /afs/.aei-potsdam.mpg.de/alpha_osf32/local/hdf/lib
HDFINCLUDES = -I$(HDFROOT)/include
HDFLIBPATH = $(HDFROOT)/lib
HDFLIBS = -L$(HDFROOT)/lib -lmfhdf -ldf -ljpeg -lz
LIBS = -L. -lieeeio
CFLAGS = $(INCLUDES) -DANSI -O 
FFLAGS = $(CFLAGS) 
CCFLAGS = $(CFLAGS)
ARCHIVE = ar crs 
OBJ_SUFF = .o
C2OBJ = -c
endif

# Linux PowerMac/Intel
ifeq ($(UNAME), Linux)
CC = gcc
C++ = g++
INCLUDES = -I.
HDFROOT = /usr/local
HDFINCLUDES = -I$(HDFROOT)/include
HDFLIBPATH = $(HDFROOT)/lib
HDFLIBS = -L$(HDFROOT)/lib -lmfhdf -ldf -ljpeg -lz
LIBS = -L. -lieeeio
CFLAGS = $(INCLUDES) -DANSI -O2 -g
FFLAGS = $(CFLAGS)
CCFLAGS = $(CFLAGS)
ARCHIVE = ar crs
OBJ_SUFF = .o
C2OBJ = -c
endif

# Macintosh / PowerMach-MachTen
ifeq ($(UNAME), machten)
CC = gcc
C++ = g++
INCLUDES = -I.
HDFROOT = 
HDFINCLUDES = -I$(HDFROOT)/include
HDFLIBPATH = $(HDFROOT)/lib
HDFLIBS = -L$(HDFROOT)/lib -lmfhdf -ldf -ljpeg -lz
LIBS = -L. -lieeeio
CFLAGS = $(INCLUDES) -DANSI -O 
FFLAGS = $(CFLAGS) 
CCFLAGS = $(CFLAGS)
ARCHIVE = ar crs 
OBJ_SUFF = .o
C2OBJ = -c
endif

# Cygwin / Win32
ifeq ($(UNAME), CYGWIN32/95)
CC = cl
C++ = cl
INCLUDES = 
HDFROOT = /afs/.aei-potsdam.mpg.de/sgi_62/local/hdf
HDFINCLUDES = -I$(HDFROOT)/include
HDFLIBPATH = $(HDFROOT)/lib
HDFLIBS = -L$(HDFROOT)/lib -lmfhdf -ldf -ljpeg -lz
LIBS = -L. -lieeeio
CFLAGS = $(INCLUDES) -DANSI -DWIN32
FFLAGS = $(CFLAGS) 
CCFLAGS = $(CFLAGS)
ARCHIVE = ar crs 
OBJ_SUFF = .obj
C2OBJ = -c /Tp
endif

# T3E
ifeq ($(UNAME), UNICOS/mk)
CC = cc
C++ = CC
INCLUDES = -I.
LIBS = -L. -lieeeio
CFLAGS = $(INCLUDES) -DANSI -DT3E -hinstantiate=used
ARCHIVE = ar crs
OBJ_SUFF = .o
C2OBJ = -c
endif

OBJECTS = IO$(OBJ_SUFF) IEEEIO$(OBJ_SUFF) 
# HLOBJECTS = Reader$(OBJ_SUFF) Writer$(OBJ_SUFF)
HLOBJECTS = Writer$(OBJ_SUFF) 
HDFOBJECTS = HDFIO$(OBJ_SUFF)
AMROBJECTS = AmrGridReader$(OBJ_SUFF) AmrUcdFileReader$(OBJ_SUFF) \
	AmrFileReader$(OBJ_SUFF) Bounds$(OBJ_SUFF) \
	AmrUcdGridHierarchy$(OBJ_SUFF) AMRwriter$(OBJ_SUFF)
TCPINCLUDES = -ITCP++ -LTCP++ 
SOCKOBJECTS = SockIOreader$(OBJ_SUFF) SockIOwriter$(OBJ_SUFF)
MPOBJECTS = MPIO$(OBJ_SUFF) MPIutils$(OBJ_SUFF)

all: archive 

objects: $(OBJECTS)

tests:  data3d.raw

archive: libieeeio.a libhlio.a libAMR.a libhdfio.a 

# always hits libTCP makefile for updates
libTCP++.a:
	(cd TCP++ ; $(MAKE) libTCP++.a)

libAMR.a: IO.hh IO.h $(AMROBJECTS)
	$(ARCHIVE) libAMR.a $(AMROBJECTS)

libhdfio.a: IO.hh IO.h HDFIO.hh $(HDFOBJECTS)
	$(ARCHIVE) libhdfio.a $(HDFOBJECTS) 

#		$(HDFLIBPATH)/libmfhdf.a \
#		$(HDFLIBPATH)/libdf.a \
#		$(HDFLIBPATH)/libjpeg.a \
#		$(HDFLIBPATH)/libz.a

libieeeio.a: $(OBJECTS) IO.h IO.hh IEEEIO.hh
	$(ARCHIVE) libieeeio.a $(OBJECTS)

libhlio.a: $(HLOBJECTS) IO.h IO.hh Writer.hh AMRwriter.hh
	$(ARCHIVE) libhlio.a $(HLOBJECTS)

libsockio.a: $(SOCKOBJECTS) SockIOwriter.hh SockIOreader.hh libTCP++.a
	$(ARCHIVE) libsockio.a $(SOCKOBJECTS)

AVSreadIEEE: libieeeio.a
	(gmake -f Makefile.AVS AVSreadIEEE)

AppendTest:	AppendTest.c libieeeio.a
	$(CC) -c $(CFLAGS) AppendTest.c
	$(C++) -o AppendTest $(CFLAGS) AppendTest.o -L. -lieeeio
	rm AppendTest.o

libmpio.a: MPIO.o MPIutils.o
	$(ARCHIVE) libmpio.a $(MPOBJECTS)

MPIutils.o: MPIutils.hh MPIutils.cc
	$(C++) -c  MPIutils.cc $(CFLAGS) 

MPIO.o: IO.hh MPIutils.hh MPIO.hh MPIO.cc Timer.hh
	$(C++) -c MPIO.cc $(CFLAGS) 

MPIOspeed: MPIOspeed.cc IEEEIO.hh MPIutils.o MPIO.o Timer.o vatoi.o libieeeio.a 	
	$(C++) $(CFLAGS) -o MPIOspeed MPIOspeed.cc MPIO.o MPIutils.o Timer.o vatoi.o  -L. -lieeeio -lmpi
 
vatoi.o: vatoi.cc vatoi.hh
	$(C++) -c vatoi.cc $(CFLAGS) 

Timer.o: Timer.cc Timer.hh
	$(C++) -c Timer.cc $(CFLAGS) 

chunkSpeed: chunkSpeed.o vatoi.o libieeeio.a 	
	$(C++) $(CFLAGS) -o chunkSpeed chunkSpeed.o vatoi.o -L. -lieeeio

chunkSpeed.o: chunkSpeed.cc IEEEIO.hh vatoi.hh	
	$(C++) $(CFLAGS) -c chunkSpeed.cc

rawSpeed: rawSpeed.o vatoi.o libieeeio.a 	
	$(C++) $(CFLAGS) -o rawSpeed rawSpeed.o vatoi.o -L. -lieeeio

rawSpeed.o: rawSpeed.cc IEEEIO.hh vatoi.hh	
	$(C++) $(CFLAGS) -c rawSpeed.cc -L. -lieeeio

mpspeed.csh: MPIOspeed chunkSpeed rawSpeed
	touch mpspeed.csh
	chmod +rx mpspeed.csh

MPItestIsend: MPItestIsend.cc MPIutils.o 
	$(C++) $(CFLAGS) -o MPItestIsend MPItestIsend.cc MPIutils.o  -L. -lmpi

data3d.raw: testllWriter
	./testllWriter
	touch data3d.raw

IO$(OBJ_SUFF): IO.cc IO.hh IO.h
	$(C++) $(CFLAGS) $(C2OBJ) IO.cc

Writer$(OBJ_SUFF): Writer.cc Writer.hh IO.hh
	$(C++) $(CFLAGS) $(C2OBJ) Writer.cc   

WildWriter$(OBJ_SUFF): WildWriter.cc WildWriter.hh IO.hh
	$(C++) $(CFLAGS) $(C2OBJ) -I/nfs/origin/data/wild/3D/ WildWriter.cc 

AmrGridReader$(OBJ_SUFF): AmrGridReader.cc AmrGridReader.hh \
	IO.hh FlexArrayTmpl.H AmrGrid.h
	$(C++) $(CFLAGS) $(C2OBJ) AmrGridReader.cc 

AmrFileReader$(OBJ_SUFF): AmrFileReader.cc AmrFileReader.hh \
	IO.hh FlexArrayTmpl.H AmrGrid.h AmrGridReader.hh
	$(C++) $(CFLAGS) $(C2OBJ) AmrFileReader.cc

AmrUcdFileReader$(OBJ_SUFF): AmrUcdFileReader.cc AmrUcdFileReader.hh \
	AmrFileReader.hh AmrUcdGridHierarchy.hh \
	IO.hh FlexArrayTmpl.H AmrGrid.h AmrGridReader.hh
	$(C++) $(CFLAGS) $(C2OBJ) AmrUcdFileReader.cc 

AmrUcdGridHierarchy$(OBJ_SUFF): AmrUcdGridHierarchy.cc \
	AmrUcdGridHierarchy.hh AmrGrid.h \
	AmrNode.hh Bounds.hh IO.hh FlexArrayTmpl.H 
	$(C++) $(CFLAGS) $(C2OBJ) AmrUcdGridHierarchy.cc 

Bounds$(OBJ_SUFF): Bounds.hh Bounds.cc
	$(C++) $(CFLAGS) $(C2OBJ) Bounds.cc

AMRwriter$(OBJ_SUFF): AMRwriter.cc AMRwriter.hh AMRwriter.h \
	Writer.hh IO.hh FlexArrayTmpl.H
	$(C++) $(CFLAGS) $(C2OBJ) AMRwriter.cc

SockIOreader$(OBJ_SUFF): SockIOreader.cc SockIOreader.hh IO.hh IEEEIO.hh
	$(C++) $(CFLAGS) $(C2OBJ) $(TCPINCLUDES) SockIOreader.cc

SockIOwriter$(OBJ_SUFF): SockIOwriter.cc SockIOwriter.hh IO.hh
	$(C++) $(CFLAGS) $(C2OBJ) $(TCPINCLUDES) SockIOwriter.cc

testSockwrite: testSockwrite.cc libieeeio.a libsockio.a
	$(C++) $(CFLAGS) $(TCPINCLUDES) -o testSockwrite testSockwrite.cc \
	-lsockio $(LIBS) -lTCP++ 

testSockread: testSockread.cc libieeeio.a libsockio.a
	$(C++) $(CFLAGS) $(TCPINCLUDES) -o testSockread testSockread.cc \
	-lsockio $(LIBS) -lTCP++ 

testLeak: testLeak.cc
	$(C++) $(CFLAGS) -o testLeak testLeak.cc $(LIBS)

IEEEIO$(OBJ_SUFF): IEEEIO.cc IEEEIO.hh IO.hh FlexArrayTmpl.H
	$(C++) $(CFLAGS) $(C2OBJ) IEEEIO.cc

HDFIO$(OBJ_SUFF): HDFIO.hh HDFIO.cc IO.hh
	$(C++) $(CCFLAG) $(CFLAGS) $(HDFINCLUDES) $(C2OBJ) HDFIO.cc

Reader$(OBJ_SUFF): Reader.cc Reader.hh IO.hh FlexArrayTmpl.H
	$(C++) $(CFLAGS) $(C2OBJ) Reader.cc   

testWriter: testWriter.cc libieeeio.a libhlio.a libhlio.a
	$(C++) $(CFLAGS) -o testWriter testWriter.cc $(LIBS) -lhlio

testllWriter: testllWriter.cc libieeeio.a libhlio.a
	$(C++) $(CFLAGS) -o testllWriter testllWriter.cc $(LIBS) -lhlio

testllReader: testllReader.cc libieeeio.a
	$(C++) $(CFLAGS) -o testllReader testllReader.cc $(LIBS)

testChunkReader: testChunkReader.cc libieeeio.a
	$(C++) $(CFLAGS) -o testChunkReader testChunkReader.cc $(LIBS)

testChunkWriter: testChunkWriter.cc libieeeio.a
	$(C++) $(CFLAGS) -o testChunkWriter testChunkWriter.cc $(LIBS)

testhdfReader: testhdfReader.cc libieeeio.a libhdfio.a
	$(C++) $(CFLAGS) $(HDFINCLUDES) -o testhdfReader testhdfReader.cc -lhdfio $(LIBS) $(HDFLIBS)

convert2native: convert2native.cc libieeeio.a
	$(C++) $(CFLAGS) -o convert2native convert2native.cc $(LIBS)

ioconvert: ioconvert.cc libieeeio.a libhdfio.a
	$(C++) $(CFLAGS) $(HDFINCLUDES) -o ioconvert ioconvert.cc -lhdfio $(LIBS) $(HDFLIBS)

ioinfo: ioinfo.cc libieeeio.a libhdfio.a
	$(C++) $(CFLAGS) $(HDFINCLUDES) -o ioinfo ioinfo.cc -lhdfio $(LIBS) $(HDFLIBS)

writef77.o: writef77.f
	$(F77) $(FFLAGS) -c writef77.f

IOspeed: IOspeed.cc writef77.o libieeeio.a libhdfio.a
	$(C++) $(CCFLAGS) $(HDFINCLUDES) -o IOspeed IOspeed.cc writef77.o -lftn $(LIBS) -lhdfio $(HDFLIBS) 

utilities: ioinfo ioconvert convert2native

testtio: testio.cc libieeeio.a
	$(C++) $(CFLAGS) -o testio testio.cc $(LIBS)

AmrUcd: libAMR.a AmrUcd.c AmrUcdCompute.cc
	(gmake -f Makefile.AVS AmrUcd)

#AVSreadIEEE: libhlio.a libIEEEIO.a
#	(gmake -f Makefile.AVS AVSreadIEEE)

AMRdata: AMRwriter.o Writer.o archive
	$(C++) $(CFLAGS) GenerateFakeAMR.cc AMRwriter.o -lhlio $(LIBS) -o GenerateFakeAMR 

GenerateFakeAMR.o: GenerateFakeAMR.cc
	$(C++) $(CFLAGS) -c GenerateFakeAMR.cc

GenerateFakeAMR: GenerateFakeAMR.o archive
	$(C++) $(CFLAGS) GenerateFakeAMR.o -lhlio $(LIBS) -o GenerateFakeAMR

readf77.o: readf77.f
	f77 $(CFLAGS) -c readf77.f

readf77: readf77.o archive
	$(C++) $(CFLAGS) readf77.o $(LIBS) -lftn -o readf77

clean:
	touch t.o ii_files Templates.DB
	rm -rf *.o ii_files Templates.DB
	touch data3d.raw testio testWriter
	rm -rf data3d.raw testio testWriter

#	(cd TCP++ ; make distclean)

distclean: clean
	touch t.a t.so
	rm -rf *.a t.so
	touch ioinfo ioconvert convert2native IOspeed testllWriter
	rm -f ioinfo ioconvert convert2native IOspeed testllWriter


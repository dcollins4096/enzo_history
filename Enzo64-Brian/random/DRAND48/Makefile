#=======================================================================
#  MAKEFILE FOR SPRNG RANDOM FIELD FILE
#  Robert Harkness, December 2005
#=======================================================================
#
# ENZO_INT is SMALL_INTS for 32 bits or LARGE_INTS for 64 bits
# ENZO_XPR is r4 (32 bits) or r8 (64 bits) for float and real
#
# FORTPR   defines how Fortran real precision is specified to the
#          compiler e.g. SGI -r8; IBM -qrealsize=8
#
# The process of building Inits executables can be controlled by setting
# the following environment variables:
#
# name       default     allowed values
# ENZO_XPR   r8          r4, r8
# ENZO_INT   LARGE_INTS  SMALL_INTS, LARGE_INTS
# ENZO_HDF   hdf5        hdf5
# ENZO_END   HDF5_BE     HDF5_BE, HDF5_LE (not recommended)
#
# The makefile will attempt to determine the operating system and
# its default message-passing library (in most cases MPI), but by
# setting the following environment variables one can easily
# experiment with new systems. All that is required is consistency
# in the naming convention, and machine-specific include files
# for HDF5.  It may also be necessary to specify FORTPR.
#
# ENZO_ARCH  auto        anything you like (e.g. cray_x1)
# ENZO_MESS  auto        mpi, vmi, mpich etc.

null :=
hdf4 := hdf4
hdf5 := hdf5

ENZO_PPR =
ENZO_XPR = r8
#ENZO_INT = SMALL_INTS
ENZO_INT = LARGE_INTS

#ENZO_MESS    = nompi

SYSTEM := $(shell uname -s)
UNAMEM := $(shell uname -m)

ifeq ($(ENZO_ARCH),$(null))

ifeq ($(SYSTEM),UNICOS/mp)
  ENZO_ARCH = crayx1
  ENZO_MESS = nompi
  ifeq ($(ENZO_XPR),r8)
    FORTPR = -s real64
  endif
endif

ifeq ($(SYSTEM),IRIX64)
  ENZO_ARCH = irix64
  ENZO_MESS = nompi
  ifeq ($(ENZO_XPR),r8)
    FORTPR = -r8 -r8const -d8
  endif
endif

ifeq ($(SYSTEM),OSF1)
  ENZO_ARCH = osf64
  ENZO_MESS = nompi
  ifeq ($(ENZO_XPR),r8)
    FORTPR = -real_size 64 -double_size 64 -fpconstant
  endif
endif

ifeq ($(SYSTEM),AIX)
  ifeq ($(OBJECT_MODE),$(null))
    ENZO_ARCH = aix32
    ENZO_MESS = nompi
  endif
  ifeq ($(OBJECT_MODE),64)
    ENZO_ARCH = aix64
    ENZO_MESS = nompi
  endif
  ifeq ($(ENZO_INT),LARGE_INTS)
    ifeq ($(ENZO_XPR),r8)
      FORTPR = -qintsize=8 -qrealsize=8 -qdpc=e
    else
      FORTPR = -qintsize=8
    endif
  else
    ifeq ($(ENZO_XPR),r8)
      FORTPR = -qintsize=4 -qrealsize=8 -qdpc=e
    else
      FORTPR = -qintsize=4
    endif
  endif
endif

ifeq ($(SYSTEM),SunOS)
  ENZO_ARCH = sun64
  ENZO_MESS = nompi
  ifeq ($(ENZO_XPR),r8)
    FORTPR = -xtypemap=real:64,double:64
  endif
endif

ifeq ($(SYSTEM),Linux)
  ifeq ($(UNAMEM),ia64)
    ENZO_ARCH = altix_v8
    ENZO_MESS = omp
  endif
  ifeq ($(UNAMEM),i686)
    ENZO_ARCH = xeon_tacc
    ENZO_MESS = nompi
  endif
  ifeq ($(ENZO_XPR),r8)
    FORTPR = -r8
  endif
endif

endif


ifeq ($(ENZO_HDF),$(null))
  ENZO_HDF := hdf5
endif

ifeq ($(ENZO_END),$(null))
  ENZO_END := HDF5_BE
endif

ifeq ($(SYSTEM),Linux)
  CPP_FLAG := -traditional
endif
ifeq ($(SYSTEM),OSF1)
  CPP_FLAG := -traditional
endif

#-------------------  object files  ------------------------------------

OBJ =	Main.o				\
	make_random_field.o		\
	write_random_number_table.o	\
	enzo_seed.o			\
	enzo_ranf.o

#-------------------  macro definitions -------------------------------

.SUFFIXES:
.SUFFIXES: .src .src90 .c .f .o .C
EXENAME = mrtab
HEADER  =
EXEDIR  =

#------------------ include machine-dependent macros -----------------

INCLUDE_MAKEFILE = $(HOME)/MakeMacros/Single/Make.${ENZO_ARCH}.${ENZO_MESS}
include ${INCLUDE_MAKEFILE}

#--------------------  implicit rules ---------------------------------

.src90.o:
	${CPP} ${DEFINES} -D${ENZO_INT} -D${ENZO_XPR} -P ${CPP_FLAG} $< > $*.f90
	${FC} ${F90SRC} ${FC_OPT} ${FREE90} ${FORTPR} ${DEBUG} ${OPTS} $*.f90
.src.o:
	${CPP} ${DEFINES} -D${ENZO_INT} -D${ENZO_XPR} -P ${CPP_FLAG} $< > $*.f
	${FC} ${FC_OPT} ${FORTPR} ${DEBUG} ${OPTS} $*.f 
.C.o:
	${CCxx} ${CCxx_OPT} -D${ENZO_INT} -D${ENZO_XPR} -D${ENZO_END} ${DEBUG} ${OPTS} $*.C
.c.o:
	${CC} ${CC_OPT} -D${ENZO_INT} -D${ENZO_XPR} -D${ENZO_END} ${DEBUG} ${OPTS} -c $<

#-------------------------  main target  -------------------------------

${EXENAME}: ${OBJ}
	${LDR} ${LDR_OPT} ${DEBUG} -o ${EXEDIR}${EXENAME} ${OBJ} ${LIB}

#---------------------  special compilation rules ---------------------


#----------------------------- help ------------------------------------

help:
	@echo Type 'make' to generate executable
	@echo Type 'make clean' to remove object files, etc

#----------------------------- clean ------------------------------------

clean:
	rm -f *.o *.mod *.f ${EXEDIR}${EXENAME}

#-----------------------------------------------------------------------

${OBJ}: ${HEADER}

#======================================================================
#
# ENZO Makefile: inits/Makefile.in
#
#======================================================================
#
# Autoconf variable values:
#
#    top_srcdir    = ..
#
#---------------------------------------------------------------------
#
# Targets:
#
#    all:          make executables
#    install:      install executables in top-level bin
#    help:         help for make targets
#    clean:        clean directory
#    distclean:    clean directory and reset configuration
#
# Also, all configuration targets available in config/MAKE.CONFIG
#
#---------------------------------------------------------------------
#
# James Bordner (jbordner@cosmos.ucsd.edu)
#
# 2003-05-19  jb  created
#
#======================================================================

TOP_DIR = ..

FLAGS_MAIN = 
DEFINES_MAIN = 

ENZO_SRC_CC = FastFourierTransform.C c_message.C jbMem.C
ENZO_SRC_FC = fourn.src power_of_2.src prefort.src fft_rotate.src f_message.src fortfft.src 
ENZO_SRC_F90 = fft90.src90 fft_f90.src90
ENZO_SRC_H = macros_and_parameters.h jbMem.h message.h
ENZO_SRC_DEF = error.def

ENZO_SRC = $(ENZO_SRC_CC) \
           $(ENZO_SRC_FC) \
           $(ENZO_SRC_F90) \
           $(ENZO_SRC_DEF) \
           $(ENZO_SRC_H)

EXES = inits

.PHONY: all
# NOTE: apparently need to re-call make, otherwise gmake doesn't recognize
#       newly-linked files from "links" target.  <*sigh*>
all: links
	$(MAKE) $(patsubst %,%.exe,$(EXES))

links: $(patsubst %,%-link,$(ENZO_SRC))
%-link:
	rm -f $*
	ln -s $(TOP_DIR)/amr_mpi/src/$* .

OBJS_MAIN = \
        $(OBJS_FFT) \
        $(OBJS_PERF) \
	FastFourierTransform.o \
	c_message.o \
	f_message.o \
	CosmologyReadParameters.o  \
	EvaluatePowerSpectrum.o \
	GenerateField.o \
	GenerateRealizationHDF4.o \
	GenerateRealizationHDF5.o \
	InitializePowerSpectrum.o \
	InterpretCommandLine.o \
	ReadParameterFile.o \
	ReadPowerSpectrumParameters.o \
	SetParameterDefaults.o \
	WriteFieldHDF4.o \
	WriteFieldHDF5.o \
	adjfft.o \
	cosmo_functions.o \
	eisenstein_power.o \
	fcol.o \
	fortfft.o \
	main.o \
	make_field.o \
	make_field_kpreserving.o \
        prefort.o \
	psfunc.o \
	qromo.o \
	ran1.o \
	shift.o \
	spline.o

DEFINES_MAIN = $(DEFINES_FFT)

#======================================================================
# Include configuration targets
#======================================================================

include $(TOP_DIR)/config/MAKE.CONFIG

# OVERRIDE: force single-precision!

include $(TOP_DIR)/config/Make.precision.single

# OVERRIDE: turn off MPI

include $(TOP_DIR)/config/Make.mpi.off

# OVERRIDE: turn off performance stuff!

include $(TOP_DIR)/config/Make.perf.off

#======================================================================
# MAIN TARGET
#======================================================================

inits.exe: $(OBJS)
	-$(LD) $(LDFLAGS) -o inits.exe $(OBJS) $(LIBS)
	($(MAKE) show-config) | grep -v directory > inits.config
	($(MAKE) show-flags) | grep -v directory > inits.flags

#======================================================================
# Pattern rules
#======================================================================

.SUFFIXES: .c .C .src .src90 .f .f90

.c.o:
	$(CC) -c $(DEFINES) $(CFLAGS) $*.c
.C.o:
	$(CXX) -c $(DEFINES) $(CXXFLAGS) $*.C
.src.f:
	$(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $*.f
.src90.f90:
	$(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $*.f90
.f.o:
	$(FC) -c $(FFLAGS) $*.f
	rm -f $*.f
.f90.o:
	$(F90) -c $(F90FLAGS) $*.f90
	rm -f $*.f90

#======================================================================
# install, clean, and help targets
#======================================================================

.PHONY: install install-%
install: $(patsubst %,install-%,$(EXES))

install-%:
	-@mv $*.exe    $(TOP_DIR)/bin/$*
	-@mv $*.config $(TOP_DIR)/bin/$*.config
	-@mv $*.flags  $(TOP_DIR)/bin/$*.flags


#----------------------------------------------------------------------

.PHONY: clean distclean
clean:
	-@rm -f *.o *.f *~ work.pc $(ENZO_SRC) singleton.mod
distclean: clean
	-@rm -f $(patsubst %,%.exe,   $(EXES))
	-@rm -f $(patsubst %,%.config,$(EXES))
	-@rm -f $(patsubst %,%.flags, $(EXES))

#----------------------------------------------------------------------

.PHONY: help
help:
	@echo ""
	@echo "INITS HELP: amr_mpi"
	@echo ""
	@echo "   make                Compile inits"
	@echo "   make install        Install inits"
	@echo "   make help           Show this help information"
	@echo "   make clean          Clean *.o *.f *~ files"
	@echo ""
	@$(MAKE) -s help-config

#======================================================================
# Dependencies
#======================================================================

# include DEPEND

#======================================================================
#
# ENZO Makefile: uniglue/Makefile.in
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
# 2003-05-20  jb  created
#
#======================================================================

TOP_DIR = ..

EXES = glop

all: $(patsubst %,%.exe,$(EXES))

OBJS_MAIN = glop.o sub.o

DEFINES_MAIN = -DIO_LOG

#======================================================================
# Include configuration targets
#======================================================================

include $(TOP_DIR)/config/MAKE.CONFIG

#======================================================================
# MAIN TARGET
#======================================================================

%.exe: $(OBJS)
	-$(LD) $(LDFLAGS) -o $*.exe $(OBJS) $(LIBS)
	($(MAKE) show-config) | grep -v directory > $*.config
	($(MAKE) show-flags) | grep -v directory > $*.flags

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
	-@rm -f *.o *.f *~ work.pc 
distclean: clean
	-@rm -f $(patsubst %,%.exe,$(EXES))
	-@rm -f $(patsubst %,%.config,$(EXES))
	-@rm -f $(patsubst %,%.flags,$(EXES))

#----------------------------------------------------------------------

.PHONY: help
help:
	@echo ""
	@echo "ENZO HELP: uniglue"
	@echo ""
	@echo "   make                Compile inits"
	@echo "   make install        Install inits"
	@echo "   make help           Show this help information"
	@echo "   make clean          Clean *.o *.f *~ files"
	@echo "   make distclean      Also clean executables *.exe
	@echo ""
	@$(MAKE) -s help-config

#======================================================================
# Dependencies
#======================================================================

glop.o: macros_and_parameters.h extern_hdf5.h
old_glop.o: macros_and_parameters.h extern_hdf5.h

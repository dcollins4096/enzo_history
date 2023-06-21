#======================================================================
#
# ENZO Makefile: hdf5_tools/Makefile.in
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
# 2003-05-16  jb  created
#
#======================================================================

TOP_DIR = ..

EXES = differ diffx slice extract dumper logger newname

all: $(patsubst %,%.exe,$(EXES))

#======================================================================
# Include configuration targets
#======================================================================

include $(TOP_DIR)/config/MAKE.CONFIG

# OVERRIDE: use single precision

include $(TOP_DIR)/config/Make.precision.single

# OVERRIDE: turn off MPI

include $(TOP_DIR)/config/Make.mpi.off

#======================================================================
# MAIN TARGETS
#======================================================================

.SUFFIXES: .exe .C

.C.exe:
	-$(CXX) -o $*.exe $(DEFINES) $(CXXFLAGS) $*.C $(LDFLAGS) $(LIBS)
	($(MAKE) show-config) | grep -v directory > $*.config
	($(MAKE) show-flags) | grep -v directory > $*.flags

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
	-@rm -f $(patsubst %,%.exe,   $(EXES))
	-@rm -f $(patsubst %,%.config,$(EXES))
	-@rm -f $(patsubst %,%.flags, $(EXES))


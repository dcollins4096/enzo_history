#=======================================================================
#
# FILE:        Makefile.config
#
# SUMMARY:     Configurable Makefile for Enzo
#
# AUTHOR:      James Bordner (jobordner@ucsd.edu)
#
# DATE:        2007-02-21
#
# DESCRIPTION
#              See 'gmake help' for definitive description of targets
#
#              Makefile.config includes the following files:
# 
#              Make.config.settings   default configuration settings
#              Make.config.override   local user settings
#              Make.config.assemble   maps 'config' settings to 'flag' settings
#              Make.config.targets    configure targets
#              Make.mach.*            all machine-dependent settings
#              Make.config.objects    list of object files
#              DEPEND                 Make-generated dependencies
#
#              Make.mach.* should be the only file that one should
#              edit when porting Enzo to a new machine.
#
#              Make.config.override should be the only file that
#              one should edit when defining local user settings.
#              Preferably, this is done implicitly through
#              the available make targets (e.g. "gmake precision-32").
#              Use "gmake help-config" for a list of all configuration
#              settings.  These make targets do error-checking; hand-editing 
#              Make.config.override is more error-prone.
#
#=======================================================================

TOP_DIR  = ../..
EXE      = enzo
THIS     = Makefile.config

ENZO_DIR             = ./
MAKE_CONFIG_OVERRIDE = $(ENZO_DIR)/Make.config.override

#-----------------------------------------------------------------------
# Make.config.settings is used for setting default values to all compile-time 
# configuration settings.
#-----------------------------------------------------------------------

include Make.config.settings

#-----------------------------------------------------------------------
# Make.config.override is used for overriding the default settings in
# Make.config.settings.  This was made separate from the default settings 
# to enable easily interfacing Enzo with a software testing environment 
# like lcatest.
#-----------------------------------------------------------------------

include $(MAKE_CONFIG_OVERRIDE)

#-----------------------------------------------------------------------
# Make.config.assemble takes the settings in the Make.config.settings
# and Make.config.override, and generates the appropriate make variables
# required by this makefile.  E.g. $(CXX), $(CXXFLAGS), etc.
#-----------------------------------------------------------------------

include Make.config.assemble

#-----------------------------------------------------------------------
# Make.mach.<machine-name> defines all machine-dependent settings.
#-----------------------------------------------------------------------

include Make.mach.$(CONFIG_MACHINE)

#=======================================================================
# OBJECT FILES
#=======================================================================

include Make.config.objects

#-----------------------------------------------------------------------
# MAIN TARGET
#-----------------------------------------------------------------------

$(EXE).exe: get-version dep $(OBJS)
	$(LD) $(LDFLAGS) -o $(EXE).exe $(OBJS) $(LIBS)

#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------

.SUFFIXES: .c .C .src .src90 .f .f90

.src.f:
	$(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $*.f
.src90.f90:
	$(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $*.f90
.f.o:
	$(FC) -c -o $*.o $(FFLAGS) $*.f
	rm -f $*.f
.f90.o:
	$(F90) -c -o $*.o $(F90FLAGS) $*.f90
	rm -f $*.f90
.c.o:
	$(CC) -c -o $*.o $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c
.C.o:
	$(CXX) -c -o $*.o $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C

#-----------------------------------------------------------------------
# Determine Enzo version
#-----------------------------------------------------------------------

.PHONY: get-version
get-version:
	@svn info $(TOP_DIR) | awk '/URL/{print "#define ENZO_SVN_BRANCH","\""$$2"\""}; /Revision/{print "#define ENZO_SVN_REVISION",$$2}' > svn_version.def
	@svn info $(TOP_DIR) | awk '/URL/{print $$2}'      > $(TOP_DIR)/BRANCH
	@svn info $(TOP_DIR) | awk '/Revision/{print $$2}' > $(TOP_DIR)/REVISION

#-----------------------------------------------------------------------
# Generate dependency file
#-----------------------------------------------------------------------

.PHONY: dep
dep:
	-@(makedepend $(DEFINES) $(INCLUDES) -fDEPEND -Y -- -- *.C *.src *.h) >& temp90210 ; touch temp90210; rm temp90210
	
include DEPEND

#-----------------------------------------------------------------------
# INSTALL TARGET
#-----------------------------------------------------------------------

install: $(TOP_DIR)/bin
	if [ ! -e $(TOP_DIR)/bin ]; then mkdir $(TOP_DIR)/bin; fi
	-@cp $(EXE).exe $(TOP_DIR)/bin/$(EXE)
	-$(MAKE) -f $(THIS) show-flags >& $(TOP_DIR)/bin/$(EXE).show-flags
	-$(MAKE) -f $(THIS) show-config >& $(TOP_DIR)/bin/$(EXE).show-config
	-$(MAKE) -f $(THIS) show-version >& $(TOP_DIR)/bin/$(EXE).show-version

$(TOP_DIR)/bin:
	mkdir $(TOP_DIR)/bin

#-----------------------------------------------------------------------
# HELP TARGET
#-----------------------------------------------------------------------

help:
	@echo
	@echo "========================================================================"
	@echo "   Enzo Makefile Help"
	@echo "========================================================================"
	@echo
	@echo "   gmake                Compile and generate the executable 'enzo.exe'"
	@echo "   gmake install        Copy the executable to bin/enzo"
	@echo "   gmake help           Display this help information"
	@echo "   gmake clean          Remove object files, executable, etc."
	@echo "   gmake dep            Create make dependencies in DEPEND file"
	@echo
	@echo "   gmake get-version    Generate version-related files"
	@echo "   gmake show-version   Display Subversion branch and revision"
	@echo
	@echo "   gmake help-config    Display detailed help on configuration make targets"
	@echo "   gmake show-config    Display the configuration settings"
	@echo "   gmake show-flags     Display specific compilation flags"
	@echo "   gmake default        Reset the configuration to the default values"
	@echo

#-----------------------------------------------------------------------

clean:
	rm -f *.o *.mod *.f *.f90 *~ $(EXE).exe

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include Make.config.targets

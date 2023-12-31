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

# Use bash since sh on datastar does not recognize ">&" used in dep: target

SHELL    = /bin/bash

TOP_DIR  = ../..
EXE      = enzo
OUTPUT   = out.compile

ENZO_DIR = .
MODULES  = 
VERBOSE  = 0

SVN      = hg

#-----------------------------------------------------------------------
# Make.config.settings is used for setting default values to all compile-time 
# configuration settings.
#-----------------------------------------------------------------------

include $(ENZO_DIR)/Make.config.settings

#-----------------------------------------------------------------------
# Make.config.machine is used for setting which Make.mach.* file to use
#-----------------------------------------------------------------------

MAKE_CONFIG_MACHINE  = $(ENZO_DIR)/Make.config.machine
include $(ENZO_DIR)/Make.config.machine

#-----------------------------------------------------------------------
# Make.config.override is used for overriding the default settings in
# Make.config.settings.  This was made separate from the default settings 
# to enable easily interfacing Enzo with a software testing environment 
# like lcatest.
#-----------------------------------------------------------------------

MAKE_CONFIG_OVERRIDE = $(ENZO_DIR)/Make.config.override

include $(MAKE_CONFIG_OVERRIDE)

# THIS WAY OF DOING THE ABOVE DOES NOT WORK:
#
# include $(ENZO_DIR)/Make.config.override

#-----------------------------------------------------------------------
# Make.config.assemble takes the settings in the Make.config.settings
# and Make.config.override, and generates the appropriate make variables
# required by this makefile.  E.g. $(CXX), $(CXXFLAGS), etc.
#-----------------------------------------------------------------------

include Make.config.assemble

#-----------------------------------------------------------------------
# Make.mach.<machine-name> defines all machine-dependent settings.
#-----------------------------------------------------------------------

-include $(ENZO_DIR)/Make.mach.$(CONFIG_MACHINE)
-include $(HOME)/.enzo/Make.mach.$(CONFIG_MACHINE)

#=======================================================================
# OBJECT FILES
#=======================================================================

include Make.config.objects

#-----------------------------------------------------------------------
# MAKE ENZO BY DEFAULT
#-----------------------------------------------------------------------

all: MACHNOTES $(EXE).exe

#-----------------------------------------------------------------------
# MAKE AN EXECUTABLE
#-----------------------------------------------------------------------

MACHNOTES: 
	@echo -e $(MACHINE_NOTES)

%.exe: $(MODULES) autogen dep %.o $(OBJS_LIB) MACHNOTES
	@rm -f $@
	@echo "Linking"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(LD) $(LDFLAGS) -o $*.exe $*.o $(OBJS_LIB) $(LIBS) >& $(OUTPUT) ; \
	else \
	  $(LD) $(LDFLAGS) -o $*.exe $*.o $(OBJS_LIB) $(LIBS) >> $(OUTPUT) \
	    2>&1 ; \
	fi)
	@(if [ -e $@ ]; then \
	   echo "Success!"; \
           if [ ! -e $(TOP_DIR)/bin ]; then mkdir $(TOP_DIR)/bin; fi; \
           cp $(EXE).exe $(TOP_DIR)/bin/$(EXE); \
	else \
	   echo "$(LD) $(LDFLAGS) -o $*.exe $*.o $(OBJS_LIB) $(LIBS)" >> temp1; \
	   cat temp1 $(OUTPUT) > temp2; \
	   rm -f temp1; \
	   mv -f temp2 $(OUTPUT); \
	   echo "Failed! See $(OUTPUT) for error messages"; \
	fi)

#-----------------------------------------------------------------------
# WRITE ALL COMPILER OUTPUT TO FILE
#-----------------------------------------------------------------------

.PHONY: verbose
verbose: VERBOSE = 1
#verbose:
#	@rm -fv $(OUTPUT)
#	@echo "Writing all compiler output to $(OUTPUT)"
verbose: $(EXE).exe

#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------

.SUFFIXES: .c .C .src .src90 .f .f90 .o .cu

# Inhibit removing any *.o files after compiling

.PRECIOUS: %.o %.f %.f90

.src.f:
	@$(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $*.f
.f.o:
	@echo "Compiling $<"
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(FC) -c -o $@ $(FFLAGS) $*.f >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
             echo "$(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $*.f"; \
             echo "$(FC) -c -o $@ $(FFLAGS) $*.f"; \
             echo; \
             $(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $*.f; \
             $(FC) -c -o $@ $(FFLAGS) $*.f; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(FC) -c -o $@ $(FFLAGS) $*.f >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)
	@(if [ -e $@ ]; then \
             rm -f $*.f; \
          fi)

.src90.f90:
	@echo "Compiling $<"
	@$(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $*.f90
.f90.o:
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(F90) -c -o $@ $(F90FLAGS) $*.f90 >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
             echo "$(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $*.f90"; \
             echo "$(F90) -c -o $@ $(F90FLAGS) $*.f90"; \
             echo; \
             $(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $*.f90; \
             $(F90) -c -o $@ $(F90FLAGS) $*.f90; \
             echo; \
             exit 1; \
	  fi ; \
	else \
	  $(F90) -c -o $@ $(F90FLAGS) $*.f90 >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)
	@(if [ -e $@ ]; then \
             rm -f $*.f90; \
          fi)

.c.o:
	@rm -f $@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(CC) -c -o $@ $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c \
	    >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
             echo "$(CC) -c -o $@ $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c"; \
             echo; \
             $(CC) -c -o $@ $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c;\
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(CC) -c -o $@ $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c \
	    >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

.C.o:
	@rm -f $@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C \
	    >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
             echo "$(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C"; \
             echo; \
             $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C;\
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C \
	    >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

.cu.o: 
	@rm -f $@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(CUDACOMPILER) -c -o $@ $(DEFINES) $(CUDACOMPFLAGS) \
	    $(INCLUDES) $*.cu >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
	     echo; \
             echo "$(CUDACOMPILER) -c -o $@ $(DEFINES) $(CUDACOMPFLAGS) $(INCLUDES) $*.cu"; \
             echo; \
             $(CUDACOMPILER)  -c -o $@ $(DEFINES) $(CUDACOMPFLAGS) $(INCLUDES) $*.cu;\
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(CUDACOMPILER) -c -o $@ $(DEFINES) $(CUDACOMPFLAGS) \
	    $(INCLUDES) $*.cu >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

#-----------------------------------------------------------------------
# Generate all make-generated source files
#-----------------------------------------------------------------------

.PHONY: autogen
autogen: svn_version.def auto_show_config.C auto_show_flags.C auto_show_version.C

# Force update of auto_show_config.C

.PHONY: auto_show_config.C
auto_show_config.C:
	-@$(MAKE) -s show-config  >& temp.show-config
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_config(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-config > auto_show_config.C

# Force update of auto_show_flags.C

.PHONY: auto_show_flags.C
auto_show_flags.C:
	-@$(MAKE) -s show-flags  >& temp.show-flags
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_flags(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-flags > auto_show_flags.C

# Force update of auto_show_version.C

.PHONY: auto_show_version.C
auto_show_version.C:
	-@$(MAKE) -s show-version  >& temp.show-version
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_version(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-version > auto_show_version.C

#-----------------------------------------------------------------------
# Generate SVN version source file 'svn_version.def'
#-----------------------------------------------------------------------

# Force update of svn_version.def
.PHONY: svn_version.def
svn_version.def:
	-@$(SVN) identify -i | \
           awk '{print "#define ENZO_SVN_REVISION","\""$$1"\""}' > svn_version.def
	-@$(SVN) identify -b | awk '{print "#define ENZO_SVN_BRANCH","\""$$1"\""};' >> svn_version.def 
	-@if [ ! -s svn_version.def ]; then \
	printf "#define ENZO_SVN_REVISION 0\n" > svn_version.def; \
	printf "#define ENZO_SVN_BRANCH \"\"\n" >> svn_version.def; \
	fi

#-----------------------------------------------------------------------
# Generate dependency file
#-----------------------------------------------------------------------

.PHONY: dep
dep:
	@echo "Updating DEPEND"
	-@(makedepend $(DEFINES) $(INCLUDES) -fDEPEND -o.o -m -- -- *.C *.c *.src *src90 *.h) >& out.make.DEPEND

include DEPEND

#-----------------------------------------------------------------------
# Radiative transfer module
#-----------------------------------------------------------------------

#include $(ENZO_DIR)/photons/Make.config.objects
#
#.PHONY: photon
#photon: OBJS_LIB += photons/*.o
#photon:
#	@echo "Making radiative transfer module"
#	+(cd photons/ ; make photon)

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
	@echo "   gmake help           Display this help information"
	@echo "   gmake clean          Remove object files, executable, etc."
	@echo "   gmake dep            Create make dependencies in DEPEND file"
	@echo
	@echo "   gmake show-version   Display revision control system branch and revision"
	@echo "   gmake show-diff      Display local file modifications"
	@echo
	@echo "   gmake help-config    Display detailed help on configuration make targets"
	@echo "   gmake show-config    Display the configuration settings"
	@echo "   gmake show-flags     Display specific compilation flags"
	@echo "   gmake default        Reset the configuration to the default values"
	@echo

#-----------------------------------------------------------------------

clean:
	-@rm -f *.o uuid/*.o *.mod *.f *.f90 DEPEND.bak *~ $(OUTPUT) *.exe \
          auto_show*.C svn_version.def hydro_rk/*.o

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include $(ENZO_DIR)/Make.config.targets

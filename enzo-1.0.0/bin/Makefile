#======================================================================
#
# ENZO bin Makefile
#
#---------------------------------------------------------------------
#
# Targets: see help target below for list
#
#---------------------------------------------------------------------
#
# James Bordner (jbordner@cosmos.ucsd.edu)
#
# 2003-05-16  created
#
#======================================================================

SHELL = tcsh

.PHONY:all
all:

ENZO_EXES = enzo emon
INITS_EXES = inits
HDF5_EXES = differ diffx extract slice dumper logger newname
ANYL_EXES = dumpgrids enzo_anyl enzohop enzostats findinit findpeaks
UNIGLUE_EXES = glop
RING_EXES = ring ringforce Xring

EXES = $(HDF5_EXES) \
       $(ENZO_EXES) \
       $(ANYL_EXES) \
       $(INITS_EXES) \
       $(UNIGLUE_EXES) \
       $(RING_EXES)

#----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm -f *~

#----------------------------------------------------------------------

.PHONY: distclean distclean-%

distclean: $(patsubst %,distclean-%,$(EXES)) clean
distclean-%:
	-@rm -f $*
	-@rm -f $*.config
	-@rm -f $*.flags

#----------------------------------------------------------------------

.PHONY: list-all list-all-% 

list-all: $(patsubst %,list-all-%,$(EXES))

list-all-%: 
	@echo $*

#----------------------------------------------------------------------

.PHONY: list-existing list-existing-% 

list-existing: $(patsubst %,list-existing-%,$(EXES))

list-existing-%: 
	@if (-e $*) echo $*

#----------------------------------------------------------------------

.PHONY: list-executable list-executable-%

list-executable: $(patsubst %,list-executable-%,$(EXES))

list-executable-%: 
	@if (-x $*) echo $*

#----------------------------------------------------------------------

.PHONY: list-missing list-missing-%

list-missing: $(patsubst %,list-missing-%,$(EXES))

list-missing-%: 
	@if (! -e $*) echo $*

#----------------------------------------------------------------------

.PHONY: help

#----------------------------------------------------------------------

.PHONY: help
help:
	@echo ""
	@echo "ENZO HELP: bin"
	@echo ""
	@echo "   make help            Show this help information"
	@echo "   make clean           Clean *~ files"
	@echo "   make distclean       Clean *~ and executables"
	@echo "   make list-all        List all potential executables"
	@echo "   make list-existing   List all actual executables"
	@echo "   make list-executable List all executable executables"
	@echo "   make list-missing    List all missing executables"
	@echo ""

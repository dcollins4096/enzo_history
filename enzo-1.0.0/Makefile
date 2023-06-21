#======================================================================
#
# ENZO Makefile
#
#======================================================================
#
# Autoconf variable values:
#
#    top_srcdir             = .
#
#---------------------------------------------------------------------
#
# Targets:
#
#    all:          make all Enzo executables
#    install:      install all Enzo executables
#    help:         help for make targets
#    clean:        clean Enzo directories
#    distclean:    clean Enzo directories and configuration
#
# Also, all configuration targets available in config/MAKE.CONFIG
#
#---------------------------------------------------------------------
#
# James Bordner (jbordner@cosmos.ucsd.edu)
#
# 2003-05-13  added target_* variables
#
#======================================================================

TOP_DIR = .
THIS_DIR = /home/dc/Enzo/enzo-1.0.0

SHELL = csh

#======================================================================
# Global targets
#======================================================================

# COMPS = config amr_mpi hdf5_tools inits ring uniglue bin config

COMPS = amr_mpi hdf5_tools inits ring uniglue

#----------------------------------------------------------------------
# Global make
#----------------------------------------------------------------------

.PHONY: all all-%
all: $(patsubst %,all-%,$(COMPS))
all-%:
	-$(MAKE) -C $*

#----------------------------------------------------------------------
# Global make install
#----------------------------------------------------------------------

.PHONY: install install-%
install: $(patsubst %,install-%,$(COMPS))
install-%:
	-$(MAKE) -C $* install

#----------------------------------------------------------------------
# Global make help
#----------------------------------------------------------------------

.PHONY: help
help:
	@echo ""
	@echo "ENZO HELP"
	@echo ""
	@echo "Components: $(COMPS)"
	@echo ""
	@echo "   make all               Make all component executables"
	@echo "   make install           Install all component executables"
	@echo "   make help              Display this help message"
	@echo "   make clean             Clean all component directories"
	@echo "   make distclean         Really clean all component directories"
	@echo ""
	@echo "   make help-config       Help on global ENZO configuration"
	@echo ""

#----------------------------------------------------------------------
# Global make clean
#----------------------------------------------------------------------

.PHONY: clean clean-%
clean: $(patsubst %,clean-%,$(COMPS) config bin)
	-@rm -f *~
clean-%:
	-$(MAKE) -C $* clean

#----------------------------------------------------------------------
# Global make distclean
#----------------------------------------------------------------------

.PHONY: distclean distclean-%
distclean: clean $(patsubst %,distclean-%,$(COMPS) config bin)
distclean-%:
	-@$(MAKE) -C $* distclean

#======================================================================
# Create tarball of given version
#======================================================================

.PHONY: tar-v% cvs-co-enzo-v%
tar-v%: 
#	Create ../enzo-temp directory
	-@$(MAKE) create-tmp
#	create temporary "version" file containing e.g. 1.0.0 instead of v1_0_0
	sed 's/v//g' VERSION | sed 's/_/./g' > version
# 	Get enzo-code-v<%> in ../enzo-temp
	-@$(MAKE) cvs-co-enzo-v$* 
#	Remove private parts (ahem)
	-@$(MAKE) remove-private-v$*
#	tar-up code
	-@$(MAKE) tar-enzo-v$* 
#	move tarball to parent directory from ../enzo-temp
	-@$(MAKE) move-tar-v$*
#	remove ../enzo-temp
	-@$(MAKE) cleantmp
TMP_DIR = $(TOP_DIR)/../enzo-temp
CVSROOT = :pserver:anonymous@jbpc.ucsd.edu:/usr/local/cvs_lca

create-tmp: cleantmp
	@echo "***DEBUG*** target create-tmp"
	-mkdir $(TMP_DIR)
cvs-co-enzo-v%:
	@echo "***DEBUG*** target cvs-co-enzo-v%"
	-cd $(TMP_DIR); setenv CVSROOT $(CVSROOT); cvs co -r v$* -d enzo-`cat $(THIS_DIR)/version` enzo-code; 
remove-private-v%:
	@echo "***DEBUG*** target remove-private"
	-cd $(TMP_DIR)/enzo-`cat $(THIS_DIR)/version`/amr_mpi/src; \
	mv -f public-DEPEND DEPEND; \
	mv -f public-Grid.h Grid.h; \
	mv -f public-main.C main.C; \
	mv -f public-mcooling.src mcooling.src; \
	mv -f public-fortran.def fortran.def; \
	mv -f public-Grid_SolveRadiativeCooling.C Grid_SolveRadiativeCooling.C; \
	mv -f public-Grid_StarParticleHandler.C Grid_StarParticleHandler.C; \
	mv -f public-Grid_ComputeCoolingTime.C Grid_ComputeCoolingTime.C; \
	mv -f public-Grid_ComputeTemperatureField.C Grid_ComputeTemperatureField.C; \
	mv -f public-ReadParameterFile.C ReadParameterFile.C ; \
	mv -f public-Makefile.in Makefile.in; \
	rm -f star_maker1.src; \
	rm -f star_maker2.src; \
	rm -f star_maker2.src.BUGS; \
	rm -f star_maker3.src; \
	rm -f star_maker4.src; \
	rm -f Gadget.C; \
	rm -f Gadget.h; \
	rm -f GadgetInitCool.C; \
	rm -f GadgetInitCoolMemory.C; \
	rm -f GadgetIonizeParams.C; \
	rm -f GadgetIonizeParamsTable.C; \
	rm -f GadgetMakeCoolingTable.C; \
	rm -f GadgetReadIonizeParams.C; \
	rm -f Grid_GadgetCalculateCooling.C; \
	rm -f Grid_GadgetComputeTemperature.C; \
	rm -f Grid_GadgetComputeTemperatureDEF.C; \
	rm -f Grid_Gadgetconvertutotemp.C; \
	rm -f Grid_GadgetCoolingRate.C; \
	rm -f Grid_GadgetCoolingRateFromU.C; \
	rm -f Grid_GadgetCoolingTime.C; \
	rm -f Grid_GadgetEquilibriumCooling.C; \
	rm -f Grid_Gadgetfindabundancesandrates.C; \
	rm -f Grid_GadgetIonizeParamsFunction.C; \
	rm -f InitializeGadgetEquilibriumCoolData.C; \
	cd ../anyl; \
	mv -f public-Makefile.in Makefile.in

tar-enzo-v%:
	@echo "***DEBUG*** target tar-enzo-v%"
	-cd $(TMP_DIR)/enzo-`cat $(THIS_DIR)/version`; rm -rf `find . -name CVS`
	-cd $(TMP_DIR); tar cf enzo-`cat $(THIS_DIR)/version`.tar enzo-`cat $(THIS_DIR)/version`; 
	-cd $(TMP_DIR); gzip enzo-`cat $(THIS_DIR)/version`.tar
move-tar-v%:
	@echo "***DEBUG*** target move-tar-v%"
	-cd $(TMP_DIR); mv enzo-`cat $(THIS_DIR)/version`.tar.gz ..
cleantmp:
	@echo "***DEBUG*** target cleantmp"
	-rm -rf $(TMP_DIR)

#======================================================================
# Include configuration targets
#======================================================================

include $(TOP_DIR)/config/MAKE.CONFIG



#=======================================================================
#  MAKEFILE FOR AMR
#    written by: Greg Bryan
#    date:       November, 1994
#    modified1:  
#
#
# MACHINE_NAME is one of: sun c3 r4400 r8000 r10000 spp alpha sp2 origin2000
# PRECISION is one of r8 (double) or r4 (single)
#

MACHINE_NAME = glamdring
PRECISION    = r4

#-------------------  object files  ------------------------------------

OBJ =   AdiabaticExpansionInitialize.o			\
	CheckEnergyConservation.o			\
	CheckForOutput.o				\
	CheckForTimeAction.o				\
	CollapseTestInitialize.o			\
	Commandline.o					\
	CommunicationBroadcastValue.o			\
	CommunicationBufferedSend.o			\
	CommunicationCombineGrids.o			\
	CommunicationInitialize.o			\
	CommunicationLoadBalanceGrids.o			\
	CommunicationParallelFFT.o			\
	CommunicationPartitionGrid.o			\
	CommunicationReceiveFluxes.o			\
	CommunicationReceiveHandler.o			\
	CommunicationSendFluxes.o			\
	CommunicationShareGrids.o			\
	CommunicationTransferParticles.o		\
	CommunicationTranspose.o			\
	CommunicationUpdateStarParticleCount.o		\
	CommunicationUtilities.o			\
	ComputePotentialFieldLevelZero.o		\
	ComputeTable.o					\
	CopyOverlappingParticleMassFields.o		\
	CopyOverlappingZones.o				\
	CosmologyComputeExpansionFactor.o		\
	CosmologyComputeExpansionTimestep.o		\
	CosmologyComputeTimeFromRedshift.o		\
	CosmologyGetUnits.o				\
	CosmologyReadParameters.o			\
	CosmologySimulationInitialize.o			\
	CosmologyWriteParameters.o			\
	DeleteFluxes.o					\
	DepositBaryons.o				\
	DepositParticleMassField.o			\
	DoubleMachInitialize.o				\
	EvolveHierarchy.o				\
	EvolveLevel.o					\
	EvolveLevelRoutinesOptimized.o			\
	ExternalBoundary_IdentifyPhysicalQuantities.o	\
	ExternalBoundary_InitializeExternalBoundaryFace.o	\
	ExternalBoundary_Prepare.o			\
	ExternalBoundary_ReadExternalBoundary.o		\
	ExternalBoundary_SetDoubleMachBoundary.o	\
	ExternalBoundary_SetExternalBoundary.o		\
	ExternalBoundary_SetExternalBoundaryParticles.o	\
	ExternalBoundary_SetShockPoolBoundary.o		\
	ExternalBoundary_SetWavePoolBoundary.o		\
	ExternalBoundary_WriteExternalBoundary.o	\
	ExternalBoundary_constructor.o			\
	ExtractSection.o				\
	FastFourierTransform.o				\
	FastFourierTransformPrepareFourn.o		\
	FastFourierTransformSGIMATH.o			\
	FastFourierTransformVECLIB.o			\
	FastSiblingLocatorInitialize.o			\
	FastSiblingLocatorFinalize.o			\
	FindField.o					\
	FindSubgrids.o					\
	GravityEquilibriumTestInitialize.o		\
	Grid_AddBaryonsToGravitatingMassField.o		\
	Grid_AddExternalAcceleration.o			\
	Grid_AddFieldMassToMassFlaggingField.o		\
	Grid_AddOverlappingParticleMassField.o		\
	Grid_AddToBoundaryFluxes.o			\
	Grid_AllocateGrids.o				\
	Grid_AnalyzeTrackPeaks.o			\
	Grid_ApplyTimeAction.o				\
	Grid_CalculateAngularMomentum.o			\
	Grid_CheckForOverlap.o				\
	Grid_CheckForPossibleOverlap.o			\
	Grid_CleanUp.o					\
	Grid_CleanUpMovedParticles.o			\
	Grid_ClearBoundaryFluxes.o			\
	Grid_ClearFlaggingField.o			\
	Grid_ClearGravitatingMassField.o		\
	Grid_ClearGravitatingMassFieldParticles.o	\
	Grid_ClearMassFlaggingField.o			\
	Grid_ClearParticleAccelerations.o		\
	Grid_CollapseTestInitializeGrid.o		\
	Grid_CollectGridInformation.o			\
	Grid_CommunicationMoveGrid.o			\
	Grid_CommunicationSendParticles.o		\
	Grid_CommunicationSendRegion.o			\
	Grid_CommunicationReceiveRegion.o		\
	Grid_CommunicationTransferParticles.o		\
	Grid_ComovingExpansionTerms.o			\
	Grid_ComovingGravitySourceTerm.o		\
	Grid_ComputeAccelerationField.o			\
	Grid_ComputeAccelerationFieldExternal.o		\
	Grid_ComputeAccelerations.o			\
	Grid_ComputeCoolingTime.o			\
	Grid_ComputeElementalDensity.o			\
	Grid_ComputeEnergy.o				\
	Grid_ComputeGammaField.o			\
	Grid_ComputePressure.o				\
	Grid_ComputePressureDualEnergyFormalism.o	\
	Grid_ComputeTemperatureField.o			\
	Grid_ComputeTimeStep.o				\
	Grid_ComputeXrayEmissivity.o			\
	Grid_ConvertTotalEnergyToGasEnergy.o		\
	Grid_CopyOverlappingMassField.o			\
	Grid_CopyParentToGravitatingFieldBoundary.o	\
	Grid_CosmologySimulationInitializeGrid.o	\
	Grid_CopyBaryonFieldToOldBaryonField.o		\
	Grid_CopyPotentialField.o			\
	Grid_CopyZonesFromGrid.o			\
	Grid_CopyZonesFromGridCountOnly.o		\
	Grid_CorrectForRefinedFluxes.o			\
	Grid_DebugCheck.o				\
	Grid_DeleteAllFields.o				\
	Grid_DeleteFlaggingField.o			\
	Grid_DepositBaryons.o				\
	Grid_DepositParticlePositions.o			\
	Grid_DoubleMachInitializeGrid.o			\
	Grid_FastSiblingLocatorAddGrid.o		\
	Grid_FastSiblingLocatorFindSiblings.o		\
	Grid_FinishFFT.o				\
	Grid_FlagBufferZones.o				\
	Grid_FlagCellsToBeRefinedByCoolingTime.o	\
	Grid_FlagCellsToBeRefinedByJeansLength.o	\
	Grid_FlagCellsToBeRefinedByMass.o		\
	Grid_FlagCellsToBeRefinedByMustRefineParticles.o \
	Grid_FlagCellsToBeRefinedBySlope.o		\
	Grid_FlagCellsToBeRefinedByShocks.o		\
	Grid_GetProjectedBoundaryFluxes.o		\
	Grid_GravityEquilibriumTestInitializeGrid.o	\
	Grid_IdentifyNewSubgrids.o			\
	Grid_IdentifyNewSubgridsSmall.o			\
	Grid_IdentifyPhysicalQuantities.o		\
	Grid_IdentifySpeciesFields.o			\
	Grid_InheritProperties.o			\
	Grid_InitializeGravitatingMassField.o		\
	Grid_InitializeGravitatingMassFieldParticles.o	\
	Grid_InitializeUniformGrid.o			\
	Grid_InterpolateAccelerations.o			\
	Grid_InterpolateBoundaryFromParent.o		\
	Grid_InterpolateFieldValues.o			\
	Grid_InterpolateParticlePositions.o		\
	Grid_InterpolatePositions.o			\
	Grid_MoveAllParticles.o				\
	Grid_MoveSubgridParticles.o			\
	Grid_MoveSubgridParticlesFast.o			\
	Grid_OutputAsParticleData.o			\
	Grid_OutputGridMovieData.o			\
	Grid_OutputStarParticleInformation.o		\
	Grid_PrepareFFT.o				\
	Grid_PrepareBoundaryFluxes.o			\
	Grid_PrepareGreensFunction.o			\
	Grid_PreparePeriodicGreensFunction.o		\
	Grid_PrepareGrid.o				\
	Grid_PrepareGridDerivedQuantities.o		\
	Grid_PreparePotentialField.o			\
	Grid_PressurelessCollapseInitialize.o		\
	Grid_ProjectSolutionToParentGrid.o		\
	Grid_ProjectToPlane.o				\
	Grid_RadiationComputeDensities.o		\
	Grid_ReadGrid.o					\
	Grid_ReadPartialGrid.o				\
	Grid_ReturnFluxDims.o				\
	Grid_ReturnGridInfo.o				\
	Grid_RestoreEnergyConsistency.o			\
	Grid_SetExternalBoundaryValues.o		\
	Grid_SetIsolatedPotentialBoundary.o		\
	Grid_SetFlaggingField.o				\
	Grid_SetFlaggingFieldStaticRegions.o		\
	Grid_SetMinimumSupport.o			\
	Grid_ShockTubeInitializeGrid.o			\
	Grid_SolveForPotential.o			\
	Grid_SolveHydroEquations.o			\
	Grid_SolveRadiativeCooling.o			\
	Grid_SolveRateEquations.o			\
	Grid_SortParticlesByNumber.o			\
	Grid_SphericalInfallGetProfile.o		\
	Grid_SphericalInfallInitializeGrid.o		\
	Grid_StarParticleHandler.o			\
	Grid_SupernovaRestartInitialize.o		\
	Grid_TestGravityCheckResults.o			\
	Grid_TestGravitySphereCheckResults.o		\
	Grid_TestGravityInitializeGrid.o		\
	Grid_TestGravityMotionInitializeGrid.o		\
	Grid_TestGravitySphereInitializeGrid.o		\
	Grid_TracerParticleCreateParticles.o		\
	Grid_TracerParticleOutputData.o			\
	Grid_TracerParticleSetVelocity.o		\
	Grid_UpdateParticlePosition.o			\
	Grid_UpdateParticleVelocity.o			\
	Grid_WriteGridDFSD.o				\
	Grid_WriteGridFlexIO.o				\
	Grid_WriteGridFortran.o				\
	Grid_WriteGridInterpolate.o			\
	Grid_WriteGridRaw.o				\
	Grid_WriteGridSD.o				\
	Grid_ZeldovichPancakeInitializeGrid.o		\
	Grid_ZeroSolutionUnderSubgrid.o			\
	Grid_constructor.o				\
	Grid_destructor.o				\
	Hierarchy_DeleteHierarchyEntry.o		\
	IdentifyNewSubgridsBySignature.o		\
	InitializeEquilibriumCoolData.o			\
	InitializeNew.o					\
	InitializeRadiationFieldData.o			\
	InitializeRateData.o				\
	InterpretCommandLine.o				\
	LevelHierarchy_AddLevel.o			\
	ListIO.o					\
	MemoryAllocationRoutines.o			\
	MultigridSolver.o				\
	OutputAsParticleData.o				\
	OutputLevelInformation.o			\
	OutputPerformanceInformation.o			\
	PrepareGravitatingMassField.o			\
	PrepareIsolatedGreensFunction.o			\
	PressurelessCollapseInitialize.o		\
	ProjectToPlane.o				\
	ProtoSubgrid_AcceptableGrid.o			\
	ProtoSubgrid_CleanUp.o				\
	ProtoSubgrid_ComputeSecondDerivative.o		\
	ProtoSubgrid_ComputeSignature.o			\
	ProtoSubgrid_CopyFlaggedZonesFromGrid.o		\
	ProtoSubgrid_CopyToNewSubgrid.o			\
	ProtoSubgrid_FindGridsByZeroSignature.o		\
	ProtoSubgrid_ReturnNthLongestDimension.o	\
	ProtoSubgrid_ShrinkToMinimumSize.o		\
	ProtoSubgrid_constructor.o			\
	QuickSortAndDrag.o				\
	RadiationFieldCalculateRates.o			\
	RadiationFieldUpdate.o				\
	ReadAllData.o					\
	ReadDataHierarchy.o				\
	ReadFlexIOFile.o				\
	ReadHDFFile.o					\
	ReadParameterFile.o				\
	ReadRadiationData.o				\
	ReadStarParticleData.o				\
	RebuildHierarchy.o				\
	ReportMemoryUsage.o				\
	ReturnCPUTime.o					\
	SetDefaultGlobalValues.o			\
	ShockInABoxInitialize.o				\
	ShockPoolInitialize.o				\
	ShockTubeInitialize.o				\
	SphericalInfallInitialize.o			\
	SupernovaRestartInitialize.o			\
	TestGravityCheckResults.o			\
	TestGravityInitialize.o				\
	TestGravityMotion.o				\
	TestGravitySphereInitialize.o			\
	TestGravitySphereCheckResults.o			\
	TracerParticleCreation.o			\
	UpdateParticlePositions.o			\
	WavePoolInitialize.o				\
	WriteAllData.o					\
	WriteDataHierarchy.o				\
	WriteEvolution.o				\
	WriteMovieData.o				\
	WriteParameterFile.o				\
	WriteRadiationData.o				\
	WriteStarParticleData.o				\
	WriteTracerParticleData.o			\
	ZeldovichPancakeInitialize.o			\
	main.o						\
        calc_dt.o                                       \
	calc_photo_rates.o				\
	calc_rad.o					\
	calc_rates.o					\
        calcdiss.o					\
	chtable.o					\
	cicinterp.o					\
	cic_deposit.o					\
	cic_flag.o					\
	cic_interp.o					\
	cie_thin_cooling_rate.o				\
	colh2diss.o					\
	coll_rates.o					\
	comp_accel.o					\
	cool1d.o					\
	cool1d_multi.o					\
	cool_multi_time.o				\
	cool_time.o					\
        euler.o						\
        expand_terms.o                                  \
	fft_utils.o					\
	fortio.o					\
	fourn.o						\
	grid_cic.o					\
	int_lin3d.o					\
	int_spline.o					\
	interp1d.o					\
	interp2d.o					\
	interp3d.o					\
	interpolate.o					\
        inteuler.o					\
        intlgrg.o                                       \
        intrmp.o                                        \
        intvar.o                                        \
        lgrg.o                                          \
	mcooling.o					\
	mg_calc_defect.o				\
	mg_prolong.o					\
	mg_prolong2.o					\
	mg_relax.o					\
	mg_restrict.o					\
	multi_cool.o					\
	ngpinterp.o					\
        pgas2d.o					\
        pgas2d_dual.o					\
        ppm_de.o					\
        ppm_lr.o					\
	prefourn.o					\
	project.o					\
	projplane.o					\
	prolong.o					\
	prolong_tsc.o					\
        remap.o                                         \
	smooth.o					\
	smooth_deposit.o				\
	solve_cool.o					\
	solve_rate.o					\
	star_maker1.o					\
	star_maker2.o					\
	star_maker3.o					\
	star_maker4.o					\
        tscint1d.o					\
        tscint2d.o					\
        tscint3d.o					\
        twoshock.o					\
	utilities.o					\
        xeuler_sweep.o					\
        xlagrange_sweep.o                               \
        yeuler_sweep.o					\
        ylagrange_sweep.o                               \
        zeuler_sweep.o                                  \
	zeus_main.o					\
	zeus_utilities.o				\
	zeus_source.o					\
	zeus_xtransport.o				\
	zeus_ytransport.o				\
	zeus_ztransport.o				\
        zlagrange_sweep.o 

#-------------------  macro definitions -------------------------------

.SUFFIXES:
.SUFFIXES: .src .c .f .o .C
HEADER  = CosmologyParameters.h ExternalBoundary.h Fluxes.h \
	  Grid.h GridList.h Hierarchy.h LevelHierarchy.h ProtoSubgrid.h \
	  TopGridData.h global_data.h macros_and_parameters.h typedefs.h \
	  CoolData.h RateData.h StarParticleData.h Commandline.h \
          RadiationFieldData.h WriteEvolution.h communication.h
EXEDIR  = ../exe/
EXENAME = enzo
#FLEXLIB = -L./FlexIO -lAMR -lhlio -lhdfio -lieeeio
#------------------ include machine-dependent macros -----------------

INCLUDE_MAKEFILE = Makefile.${MACHINE_NAME}
include ${INCLUDE_MAKEFILE}

#--------------------  implicit rules ---------------------------------

.src.o:
	${CPP} ${DEFINES} -P $< > $*.f
	${FC} ${FC_OPTIMIZE} ${DEBUG} ${OPTS} $*.f 
.C.o:
	${CCxx} ${CCxx_OPTIMIZE} -D ${PRECISION} ${DEBUG} ${OPTS} $*.C
.c.o:
	${CC} ${CC_OPTIMIZE} -D ${PRECISION} ${DEBUG} ${OPTS} -c $<

#-------------------------  main target  -------------------------------

${EXENAME}: ${OBJ}
#	(cd FlexIO ; gmake)
	${LDR} ${LDR_OPTIMIZE} ${DEBUG} -o ${EXEDIR}${EXENAME} ${OBJ} \
	${FLEXLIB} ${LIB}

#---------------------  special compilation rules ---------------------


#----------------------------- help ------------------------------------

help:
	@echo Type 'make' to generate executable
	@echo Type 'make clean' to remove object files, etc

#----------------------------- clean ------------------------------------

clean:
	rm *.o *.f *~
#	(cd FlexIO ; gmake clean)

# gets rid of objects *and* libraries *and* executables
distclean: clean
	(cd FlexIO ; gmake distclean)

#-----------------------------------------------------------------------

${OBJ}: ${HEADER}

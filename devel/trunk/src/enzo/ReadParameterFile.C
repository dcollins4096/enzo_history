/***********************************************************************
/
/  READ A PARAMETER FILE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       June 25, 2006
/  modified2:  Robert Harkness
/  date:       February 29th, 2008
/  modified3:  Robert Harkness
/  date:       May, 2008
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine reads the parameter file in the argument and sets parameters
//   based on it.
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "StarParticleData.h"
 
 
/* This variable is only declared and used in Grid_DepositPositions. */
 
extern float DepositPositionsParticleSmoothRadius;
 
/* This variable is declared here and only used in Grid_ReadGrid. */
 
// int ParticleTypeInFile = 0;
 
/* function prototypes */
 
int ReadListOfFloats(FILE *fptr, int N, float floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
int CosmologyReadParameters(FILE *fptr, FLOAT *StopTime, FLOAT *InitTime);
int ReadUnits(FILE *fptr);
int InitializeRateData(FLOAT Time);
int InitializeEquilibriumCoolData(FLOAT Time);
int InitializeGadgetEquilibriumCoolData(FLOAT Time);
int InitializeRadiationFieldData(FLOAT Time);
 
 
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt)
{
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int dim, ret, int_dummy;
  float TempFloat;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
 
  /* read until out of lines */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read MetaData parameters */
 
    ret += sscanf(line, "InitialCycleNumber = %"ISYM, &MetaData.CycleNumber);
    ret += sscanf(line, "InitialTime        = %"PSYM, &MetaData.Time);
    ret += sscanf(line, "InitialCPUTime     = %"FSYM, &MetaData.CPUTime);
    ret += sscanf(line, "Initialdt          = %"FSYM, Initialdt);
 
    ret += sscanf(line, "StopTime    = %"PSYM, &MetaData.StopTime);
    ret += sscanf(line, "StopCycle   = %"ISYM, &MetaData.StopCycle);
    ret += sscanf(line, "StopSteps   = %"ISYM, &MetaData.StopSteps);
    ret += sscanf(line, "StopCPUTime = %"FSYM, &MetaData.StopCPUTime);
 
    ret += sscanf(line, "TimeLastRestartDump = %"PSYM,
		  &MetaData.TimeLastRestartDump);
    ret += sscanf(line, "dtRestartDump       = %"PSYM, &MetaData.dtRestartDump);
    ret += sscanf(line, "TimeLastDataDump    = %"PSYM,
		  &MetaData.TimeLastDataDump);
    ret += sscanf(line, "dtDataDump          = %"PSYM, &MetaData.dtDataDump);
    ret += sscanf(line, "TimeLastHistoryDump = %"PSYM,
		  &MetaData.TimeLastHistoryDump);
    ret += sscanf(line, "dtHistoryDump       = %"PSYM, &MetaData.dtHistoryDump);
    ret += sscanf(line, "TimeLastMovieDump = %"PSYM,
		  &MetaData.TimeLastMovieDump);
    ret += sscanf(line, "dtMovieDump       = %"PSYM, &MetaData.dtMovieDump);
 
    ret += sscanf(line, "TracerParticleOn  = %"ISYM, &TracerParticleOn);
    ret += sscanf(line, "ParticleTypeInFile = %"ISYM, &ParticleTypeInFile);
    ret += sscanf(line, "TimeLastTracerParticleDump = %"PSYM,
                  &MetaData.TimeLastTracerParticleDump);
    ret += sscanf(line, "dtTracerParticleDump       = %"PSYM,
                  &MetaData.dtTracerParticleDump);
 
    ret += sscanf(line, "MovieRegionLeftEdge  = %"PSYM" %"PSYM" %"PSYM,
		  MetaData.MovieRegionLeftEdge,
		  MetaData.MovieRegionLeftEdge+1,
		  MetaData.MovieRegionLeftEdge+2);
    ret += sscanf(line, "MovieRegionRightEdge = %"PSYM" %"PSYM" %"PSYM,
		  MetaData.MovieRegionRightEdge,
		  MetaData.MovieRegionRightEdge+1,
		  MetaData.MovieRegionRightEdge+2);
 

    ret += sscanf(line, "NewMovieLeftEdge  = %"FSYM" %"FSYM" %"FSYM, 
		  MetaData.NewMovieLeftEdge,
		  MetaData.NewMovieLeftEdge+1, 
		  MetaData.NewMovieLeftEdge+2);
    ret += sscanf(line, "NewMovieRightEdge = %"FSYM" %"FSYM" %"FSYM, 
		  MetaData.NewMovieRightEdge, 
		  MetaData.NewMovieRightEdge+1,
		  MetaData.NewMovieRightEdge+2);

    ret += sscanf(line, "CycleLastRestartDump = %"ISYM,
		  &MetaData.CycleLastRestartDump);
    ret += sscanf(line, "CycleSkipRestartDump = %"ISYM,
		  &MetaData.CycleSkipRestartDump);
    ret += sscanf(line, "CycleLastDataDump    = %"ISYM,
		  &MetaData.CycleLastDataDump);
    ret += sscanf(line, "CycleSkipDataDump    = %"ISYM,
		  &MetaData.CycleSkipDataDump);
    ret += sscanf(line, "CycleLastHistoryDump = %"ISYM,
		  &MetaData.CycleLastHistoryDump);
    ret += sscanf(line, "CycleSkipHistoryDump = %"ISYM,
		  &MetaData.CycleSkipHistoryDump);
    ret += sscanf(line, "CycleSkipGlobalDataDump = %"ISYM, //AK
                  &MetaData.CycleSkipGlobalDataDump);
    ret += sscanf(line, "OutputFirstTimeAtLevel = %"ISYM,
		  &MetaData.OutputFirstTimeAtLevel);
    ret += sscanf(line, "StopFirstTimeAtLevel = %"ISYM,
		  &MetaData.StopFirstTimeAtLevel);
 
    ret += sscanf(line, "RestartDumpNumber = %"ISYM, &MetaData.RestartDumpNumber);
    ret += sscanf(line, "DataDumpNumber    = %"ISYM, &MetaData.DataDumpNumber);
    ret += sscanf(line, "HistoryDumpNumber = %"ISYM, &MetaData.HistoryDumpNumber);
    ret += sscanf(line, "MovieDumpNumber   = %"ISYM, &MetaData.MovieDumpNumber);
    ret += sscanf(line, "TracerParticleDumpNumber = %"ISYM, &MetaData.TracerParticleDumpNumber);
 
    if (sscanf(line, "RestartDumpName      = %s", dummy) == 1)
      MetaData.RestartDumpName = dummy;
    if (sscanf(line, "DataDumpName         = %s", dummy) == 1)
      MetaData.DataDumpName = dummy;
    if (sscanf(line, "HistoryDumpName      = %s", dummy) == 1)
      MetaData.HistoryDumpName = dummy;
    if (sscanf(line, "MovieDumpName        = %s", dummy) == 1)
      MetaData.MovieDumpName = dummy;
    if (sscanf(line, "TracerParticleDumpName = %s", dummy) == 1)
      MetaData.TracerParticleDumpName = dummy;
    if (sscanf(line, "RedshiftDumpName     = %s", dummy) == 1)
      MetaData.RedshiftDumpName = dummy;
 
    if (sscanf(line, "RestartDumpDir      = %s", dummy) == 1)
      MetaData.RestartDumpDir = dummy;
    if (sscanf(line, "DataDumpDir         = %s", dummy) == 1)
      MetaData.DataDumpDir = dummy;
    if (sscanf(line, "HistoryDumpDir      = %s", dummy) == 1)
      MetaData.HistoryDumpDir = dummy;
    if (sscanf(line, "MovieDumpDir        = %s", dummy) == 1)
      MetaData.MovieDumpDir = dummy;
    if (sscanf(line, "TracerParticleDumpDir = %s", dummy) == 1)
      MetaData.TracerParticleDumpDir = dummy;
    if (sscanf(line, "RedshiftDumpDir     = %s", dummy) == 1)
      MetaData.RedshiftDumpDir = dummy;
 
    if (sscanf(line, "LocalDir            = %s", dummy) == 1)
      MetaData.LocalDir = dummy;
    if (sscanf(line, "GlobalDir           = %s", dummy) == 1)
      MetaData.GlobalDir = dummy;
 
    if (sscanf(line, "CubeDump[%"ISYM"] = %s", &dim, dummy) == 2) {
      ret++; CubeDumps[dim] = dummy;
      if (dim >= MAX_CUBE_DUMPS) {
        fprintf(stderr, "CubeDump %"ISYM" > maximum allowed.\n", dim);
        return FAIL;
      }
    }
 
    if (sscanf(line, "TimeActionType[%"ISYM"] = %"ISYM, &dim, &int_dummy) == 2) {
      ret++; TimeActionType[dim] = int_dummy;
      if (dim >= MAX_TIME_ACTIONS-1) {
	fprintf(stderr, "Time action %"ISYM" > maximum allowed.\n", dim);
	return FAIL;
      }
    }
    if (sscanf(line, "TimeActionRedshift[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionRedshift[%"ISYM"] = %"PSYM, &dim,
		    TimeActionRedshift+dim);
    if (sscanf(line, "TimeActionTime[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionTime[%"ISYM"] = %"PSYM, &dim,
		    TimeActionTime+dim);
    if (sscanf(line, "TimeActionParameter[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionParameter[%"ISYM"] = %"FSYM, &dim,
		    TimeActionParameter+dim);
 
    ret += sscanf(line, "StaticHierarchy = %"ISYM, &MetaData.StaticHierarchy);
 
    ret += sscanf(line, "TopGridRank       = %"ISYM, &MetaData.TopGridRank);
    ret += sscanf(line, "TopGridDimensions = %"ISYM" %"ISYM" %"ISYM, MetaData.TopGridDims,
		  MetaData.TopGridDims+1, MetaData.TopGridDims+2);
 
    ret += sscanf(line, "TopGridGravityBoundary = %"ISYM,
		  &MetaData.GravityBoundary);
 
    ret += sscanf(line, "ParticleBoundaryType   = %"ISYM,
		  &MetaData.ParticleBoundaryType);
    ret += sscanf(line, "NumberOfParticles      = %"ISYM,
		  &MetaData.NumberOfParticles);
 
    ret += sscanf(line, "CourantSafetyNumber    = %"FSYM,
		  &MetaData.CourantSafetyNumber);
    ret += sscanf(line, "PPMFlatteningParameter = %"ISYM,
		  &MetaData.PPMFlatteningParameter);
    ret += sscanf(line, "PPMDiffusionParameter  = %"ISYM,
		  &MetaData.PPMDiffusionParameter);
    ret += sscanf(line, "PPMSteepeningParameter = %"ISYM,
		  &MetaData.PPMSteepeningParameter);
 
    /* read global Parameters */
 
    ret += sscanf(line, "ProblemType            = %"ISYM, &ProblemType);
    ret += sscanf(line, "HydroMethod            = %"ISYM, &HydroMethod);
    ret += sscanf(line, "huge_number            = %"FSYM, &huge_number);
    ret += sscanf(line, "tiny_number            = %"FSYM, &tiny_number);
    ret += sscanf(line, "Gamma                  = %"FSYM, &Gamma);
    ret += sscanf(line, "PressureFree           = %"ISYM, &PressureFree);
    ret += sscanf(line, "RefineBy               = %"ISYM, &RefineBy);
    ret += sscanf(line, "MaximumRefinementLevel = %"ISYM,
		  &MaximumRefinementLevel);
    ret += sscanf(line, "MaximumGravityRefinementLevel = %"ISYM,
		  &MaximumGravityRefinementLevel);
    ret += sscanf(line, "MaximumParticleRefinementLevel = %"ISYM,
		  &MaximumParticleRefinementLevel);
    ret += sscanf(line, "CellFlaggingMethod     = %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM,
	     CellFlaggingMethod+0, CellFlaggingMethod+1, CellFlaggingMethod+2,
	     CellFlaggingMethod+3, CellFlaggingMethod+4, CellFlaggingMethod+5,
	     CellFlaggingMethod+6);
    ret += sscanf(line, "FluxCorrection         = %"ISYM, &FluxCorrection);
    ret += sscanf(line, "InterpolationMethod    = %"ISYM, &InterpolationMethod);
    ret += sscanf(line, "ConservativeInterpolation = %"ISYM,
		  &ConservativeInterpolation);
    ret += sscanf(line, "MinimumEfficiency      = %"FSYM, &MinimumEfficiency);
    ret += sscanf(line, "MinimumSubgridEdge     = %"ISYM, &MinimumSubgridEdge);
    ret += sscanf(line, "MaximumSubgridSize     = %"ISYM, &MaximumSubgridSize);
    ret += sscanf(line, "NumberOfBufferZones    = %"ISYM, &NumberOfBufferZones);
 
    ret += sscanf(line, "DomainLeftEdge        = %"PSYM" %"PSYM" %"PSYM, DomainLeftEdge,
		  DomainLeftEdge+1, DomainLeftEdge+2);
    ret += sscanf(line, "DomainRightEdge       = %"PSYM" %"PSYM" %"PSYM, DomainRightEdge,
		  DomainRightEdge+1, DomainRightEdge+2);
    ret += sscanf(line, "GridVelocity          = %"FSYM" %"FSYM" %"FSYM, GridVelocity,
		  GridVelocity+1, GridVelocity+2);
    ret += sscanf(line, "RefineRegionLeftEdge  = %"PSYM" %"PSYM" %"PSYM,
		  RefineRegionLeftEdge, RefineRegionLeftEdge+1,
		  RefineRegionLeftEdge+2);
    ret += sscanf(line, "RefineRegionRightEdge = %"PSYM" %"PSYM" %"PSYM,
		  RefineRegionRightEdge, RefineRegionRightEdge+1,
		  RefineRegionRightEdge+2);
 
    if (sscanf(line, "DataLabel[%"ISYM"] = %s\n", &dim, dummy) == 2)
      DataLabel[dim] = dummy;
    if (sscanf(line, "DataUnits[%"ISYM"] = %s\n", &dim, dummy) == 2)
      DataUnits[dim] = dummy;
 
    ret += sscanf(line, "UniformGravity          = %"ISYM, &UniformGravity);
    ret += sscanf(line, "UniformGravityDirection = %"ISYM,
		  &UniformGravityDirection);
    ret += sscanf(line, "UniformGravityConstant  = %"FSYM,
		  &UniformGravityConstant);
 
    ret += sscanf(line, "PointSourceGravity         = %"ISYM,&PointSourceGravity);
    ret += sscanf(line, "PointSourceGravityPosition = %"PSYM" %"PSYM" %"PSYM,
		  PointSourceGravityPosition, PointSourceGravityPosition+1,
		  PointSourceGravityPosition+2);
    ret += sscanf(line, "PointSourceGravityConstant = %"FSYM,
		  &PointSourceGravityConstant);
    ret += sscanf(line, "PointSourceGravityCoreRadius = %"FSYM,
		  &PointSourceGravityCoreRadius);
 
    ret += sscanf(line, "SelfGravity           = %"ISYM, &SelfGravity);
    ret += sscanf(line, "GravitationalConstant = %"FSYM, &GravitationalConstant);
    ret += sscanf(line, "S2ParticleSize        = %"FSYM, &S2ParticleSize);
    ret += sscanf(line, "GravityResolution     = %"FSYM, &GravityResolution);
    ret += sscanf(line, "ComputePotential      = %"ISYM, &ComputePotential);
    ret += sscanf(line, "WritePotential        = %"ISYM, &WritePotential);
    ret += sscanf(line, "BaryonSelfGravityApproximation = %"ISYM,
		  &BaryonSelfGravityApproximation);
 
    ret += sscanf(line, "GreensFunctionMaxNumber   = %"ISYM,
		  &GreensFunctionMaxNumber);
    ret += sscanf(line, "GreensFunctionMaxSize     = %"ISYM,
		  &GreensFunctionMaxSize);
 
    ret += sscanf(line, "DualEnergyFormalism     = %"ISYM, &DualEnergyFormalism);
    ret += sscanf(line, "DualEnergyFormalismEta1 = %"FSYM,
		  &DualEnergyFormalismEta1);
    ret += sscanf(line, "DualEnergyFormalismEta2 = %"FSYM,
		  &DualEnergyFormalismEta2);
    ret += sscanf(line, "ParticleCourantSafetyNumber = %"FSYM,
		  &ParticleCourantSafetyNumber);
    ret += sscanf(line, "RandomForcing = %"ISYM, &RandomForcing); //AK
    ret += sscanf(line, "RandomForcingEdot = %"FSYM, &RandomForcingEdot); //AK
    ret += sscanf(line, "RandomForcingMachNumber = %"FSYM, //AK
                  &RandomForcingMachNumber);
    ret += sscanf(line, "RadiativeCooling = %"ISYM, &RadiativeCooling);
    ret += sscanf(line, "GadgetEquilibriumCooling = %"ISYM, &GadgetEquilibriumCooling);
    ret += sscanf(line, "MultiSpecies = %"ISYM, &MultiSpecies);
    ret += sscanf(line, "RadiationFieldType = %"ISYM, &RadiationFieldType);
    ret += sscanf(line, "AdjustUVBackground = %"ISYM, &AdjustUVBackground);
    ret += sscanf(line, "SetUVBAmplitude = %"FSYM, &SetUVBAmplitude);
    ret += sscanf(line, "SetHeIIHeatingScale = %"FSYM, &SetHeIIHeatingScale);

    ret += sscanf(line, "RadiationFieldLevelRecompute = %"ISYM,
		  &RadiationFieldLevelRecompute);
    ret += sscanf(line, "RadiationSpectrumNormalization = %"FSYM,
		  &CoolData.f3);
    ret += sscanf(line, "RadiationSpectrumSlope = %"FSYM, &CoolData.alpha0);

    if (sscanf(line, "CoolDataParameterFile = %s", dummy) == 1)
      CoolData.ParameterFilename = dummy;

    ret += sscanf(line, "ZEUSQuadraticArtificialViscosity = %"FSYM,
		  &ZEUSQuadraticArtificialViscosity);
    ret += sscanf(line, "ZEUSLinearArtificialViscosity = %"FSYM,
		  &ZEUSLinearArtificialViscosity);
 
    ret += sscanf(line, "UseMinimumPressureSupport = %"ISYM,
		  &UseMinimumPressureSupport);
    ret += sscanf(line, "MinimumPressureSupportParameter = %"FSYM,
		  &MinimumPressureSupportParameter);
    ret += sscanf(line, "RefineByJeansLengthSafetyFactor = %"FSYM,
		  &RefineByJeansLengthSafetyFactor);
    ret += sscanf(line, "MustRefineParticlesRefineToLevel = %"ISYM,
                  &MustRefineParticlesRefineToLevel);
    ret += sscanf(line, "ParticleTypeInFile = %"ISYM,
                  &ParticleTypeInFile);
 
    if (sscanf(line, "StaticRefineRegionLevel[%"ISYM"] = %"ISYM,&dim,&int_dummy) == 2){
      if (dim > MAX_STATIC_REGIONS-1) {
        fprintf(stderr, "StaticRegion number %"ISYM" > MAX allowed\n", dim);
        return FAIL;
      }
      ret++;
      StaticRefineRegionLevel[dim] = int_dummy;
    }
    if (sscanf(line, "StaticRefineRegionLeftEdge[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "StaticRefineRegionLeftEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, StaticRefineRegionLeftEdge[dim],
		    StaticRefineRegionLeftEdge[dim]+1,
		    StaticRefineRegionLeftEdge[dim]+2);
    if (sscanf(line, "StaticRefineRegionRightEdge[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "StaticRefineRegionRightEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, StaticRefineRegionRightEdge[dim],
		    StaticRefineRegionRightEdge[dim]+1,
		    StaticRefineRegionRightEdge[dim]+2);
 
    ret += sscanf(line, "ParallelRootGridIO = %"ISYM, &ParallelRootGridIO);
 
    ret += sscanf(line, "ParallelParticleIO = %"ISYM, &ParallelParticleIO);
 
    ret += sscanf(line, "Unigrid = %"ISYM, &Unigrid);
 
    ret += sscanf(line, "PartitionNestedGrids = %"ISYM, &PartitionNestedGrids);
 
    ret += sscanf(line, "ExtractFieldsOnly = %"ISYM, &ExtractFieldsOnly);
 
    ret += sscanf(line, "CubeDumpEnabled = %"ISYM, &CubeDumpEnabled);
 
    if (sscanf(line, "SRBprefix = %s", dummy) == 1)
      SRBprefix = dummy;

    ret += sscanf(line, "Debug1 = %"ISYM, &debug1);

    ret += sscanf(line, "Debug2 = %"ISYM, &debug2);

    ret += sscanf(line, "MemoryLimit = %"ISYM, &MemoryLimit);

#ifdef STAGE_INPUT
    ret += sscanf(line, "StageInput = %"ISYM, &StageInput);
#endif

#ifdef OOC_BOUNDARY

    ret += sscanf(line, "ExternalBoundaryIO = %"ISYM, &ExternalBoundaryIO);

    ret += sscanf(line, "ExternalBoundaryTypeIO = %"ISYM, &ExternalBoundaryTypeIO);

    ret += sscanf(line, "ExternalBoundaryValueIO = %"ISYM, &ExternalBoundaryValueIO);

    ret += sscanf(line, "SimpleConstantBoundary = %"ISYM, &SimpleConstantBoundary);

#endif
 
    ret += sscanf(line, "MinimumOverDensityForRefinement  = "
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
		  MinimumOverDensityForRefinement+0, 
		  MinimumOverDensityForRefinement+1,
		  MinimumOverDensityForRefinement+2, 
		  MinimumOverDensityForRefinement+3,
		  MinimumOverDensityForRefinement+4, 
		  MinimumOverDensityForRefinement+5,
		  MinimumOverDensityForRefinement+6);
    ret += sscanf(line, "MinimumMassForRefinement  = "
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
		  MinimumMassForRefinement+0, 
		  MinimumMassForRefinement+1,
		  MinimumMassForRefinement+2, 
		  MinimumMassForRefinement+3,
		  MinimumMassForRefinement+4, 
		  MinimumMassForRefinement+5,
		  MinimumMassForRefinement+6);
    ret += sscanf(line, "MinimumMassForRefinementLevelExponent = "
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
		  MinimumMassForRefinementLevelExponent+0,
		  MinimumMassForRefinementLevelExponent+1,
		  MinimumMassForRefinementLevelExponent+2,
		  MinimumMassForRefinementLevelExponent+3,
		  MinimumMassForRefinementLevelExponent+4,
		  MinimumMassForRefinementLevelExponent+5,
		  MinimumMassForRefinementLevelExponent+6);
    ret += sscanf(line, "MinimumSlopeForRefinement = %"FSYM,
		  &MinimumSlopeForRefinement);
    ret += sscanf(line, "MinimumPressureJumpForRefinement = %"FSYM,
		  &MinimumPressureJumpForRefinement);
    ret += sscanf(line, "MinimumShearForRefinement = %"FSYM,
		  &MinimumShearForRefinement);
    ret += sscanf(line, "MinimumEnergyRatioForRefinement = %"FSYM,
		  &MinimumEnergyRatioForRefinement);
    ret += sscanf(line, "ComovingCoordinates = %"ISYM,&ComovingCoordinates);
    ret += sscanf(line, "StarParticleCreation = %"ISYM, &StarParticleCreation);
    ret += sscanf(line, "StarParticleFeedback = %"ISYM, &StarParticleFeedback);
    ret += sscanf(line, "NumberOfParticleAttributes = %"ISYM,
		  &NumberOfParticleAttributes);
 
    /* read data which defines the boundary conditions */
 
    ret += sscanf(line, "LeftFaceBoundaryCondition  = %"ISYM" %"ISYM" %"ISYM,
		  MetaData.LeftFaceBoundaryCondition,
		  MetaData.LeftFaceBoundaryCondition+1,
		  MetaData.LeftFaceBoundaryCondition+2);
    ret += sscanf(line, "RightFaceBoundaryCondition = %"ISYM" %"ISYM" %"ISYM,
		  MetaData.RightFaceBoundaryCondition,
		  MetaData.RightFaceBoundaryCondition+1,
		  MetaData.RightFaceBoundaryCondition+2);
    if (sscanf(line, "BoundaryConditionName         = %s", dummy) == 1)
      MetaData.BoundaryConditionName = dummy;
 
    /* Check version number. */
 
    if (sscanf(line, "VersionNumber = %"FSYM, &TempFloat) == 1) {
      ret++;
      if (fabs(TempFloat - VERSION) >= 1.0e-3)
	fprintf(stderr, "Warning: Incorrect version number.\n");
    }
 
    /* Read star particle parameters. */
 
    ret += sscanf(line, "StarMakerOverDensityThreshold = %"FSYM,
		  &StarMakerOverDensityThreshold);
    ret += sscanf(line, "StarMakerMassEfficiency = %"FSYM,
		  &StarMakerMassEfficiency);
    ret += sscanf(line, "StarMakerMinimumMass = %"FSYM, &StarMakerMinimumMass);
    ret += sscanf(line, "StarMakerMinimumDynamicalTime = %"FSYM,
                  &StarMakerMinimumDynamicalTime);
    ret += sscanf(line, "StarMassEjectionFraction = %"FSYM,
		  &StarMassEjectionFraction);
    ret += sscanf(line, "StarMetalYield = %"FSYM, &StarMetalYield);
    ret += sscanf(line, "StarEnergyToThermalFeedback = %"FSYM,
		  &StarEnergyToThermalFeedback);
    ret += sscanf(line, "StarEnergyToStellarUV = %"FSYM, &StarEnergyToStellarUV);
    ret += sscanf(line, "StarEnergyToQuasarUV = %"FSYM, &StarEnergyToQuasarUV);
 
    /* Read Movie Dump parameters */

    ret += sscanf(line, "MovieSkipTimestep = %"ISYM, &MovieSkipTimestep);
    ret += sscanf(line, "NewMovieParticleOn = %"ISYM, &NewMovieParticleOn);
    ret += sscanf(line, "MovieDataField = %"ISYM" %"ISYM" %"ISYM,
                  MovieDataField+0, MovieDataField+1, MovieDataField+2);
    ret += sscanf(line, "NewMovieDumpNumber = %"ISYM, &NewMovieDumpNumber);
    if (sscanf(line, "NewMovieName = %s", dummy) == 1)
      NewMovieName = dummy;

    ret += sscanf(line, "MultiMetals = %"ISYM, &MultiMetals);

#ifdef STAGE_INPUT
    sscanf(line, "LocalPath = %s\n", LocalPath);
    sscanf(line, "GlobalPath = %s\n", GlobalPath);
#endif
 
    /* If the dummy char space was used, then make another. */
 
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
    /* check to see if the line belongs to one of the test problems */
 
    if (strstr(line, "ShockTube")           ) ret++;
    if (strstr(line, "WavePool" )           ) ret++;
    if (strstr(line, "ShockPool")           ) ret++;
    if (strstr(line, "DoubleMach")          ) ret++;
    if (strstr(line, "Implosion")           ) ret++;
    if (strstr(line, "SedovBlast")          ) ret++;
    if (strstr(line, "Units")               ) ret++;
    if (strstr(line, "KelvinHelmholtz")     ) ret++;
    if (strstr(line, "KH")                  ) ret++;
    if (strstr(line, "Noh")                 ) ret++;
    if (strstr(line, "ZeldovichPancake")    ) ret++;
    if (strstr(line, "PressurelessCollapse")) ret++;
    if (strstr(line, "AdiabaticExpansion" ) ) ret++;
    if (strstr(line, "CosmologySimulation") ) ret++;
    if (strstr(line, "TestGravity"        ) ) ret++;
    if (strstr(line, "SphericalInfall"    ) ) ret++;
    if (strstr(line, "TestGravitySphere"  ) ) ret++;
    if (strstr(line, "CollapseTest"       ) ) ret++;
    if (strstr(line, "Cosmology"          ) ) ret++;
    if (strstr(line, "SupernovaRestart"   ) ) ret++;
    if (strstr(line, "TracerParticleCreation")) ret++;
    if (strstr(line, "TurbulenceSimulation")) ret++;
    if (strstr(line, "ProtostellarCollapse")) ret++;
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#')
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s", line);
 
  }
 
  /* clean up */
 
  delete dummy;
  rewind(fptr);
 
  /* If we have turned on Comoving coordinates, read cosmology parameters. */
 
  if (ComovingCoordinates) {
    if (CosmologyReadParameters(fptr, &MetaData.StopTime, &MetaData.Time)
	== FAIL) {
      fprintf(stderr, "Error in ReadCosmologyParameters.\n");;
      return FAIL;
    }
    rewind(fptr);
  }
  else {
    if (ReadUnits(fptr) == FAIL){
      fprintf(stderr, "Error in ReadUnits. \n");
      return FAIL;
    }
    rewind(fptr);
  }
 
  /* If GadgetEquilibriumCooling == TRUE, we don't want MultiSpecies
     or RadiationFieldType to be on - both are taken care of in
     the Gadget cooling routine.  Therefore, we turn them off!
     Also, initialize the Gadget equilibrium cooling data. */

  if(GadgetEquilibriumCooling == TRUE){

    if(MyProcessorNumber == ROOT_PROCESSOR ) {
      fprintf(stderr, "WARNING:  GadgetEquilibriumCooling = 1.  Forcing\n");
      fprintf(stderr, "WARNING:  RadiationFieldType = 0, MultiSpecies = 0, and\n");
      fprintf(stderr, "WARNING:  RadiativeCooling = 1.\n");
    }

    RadiationFieldType = 0;
    MultiSpecies       = 0;
    RadiativeCooling   = 1;

    // initialize Gadget equilibrium cooling
    if (InitializeGadgetEquilibriumCoolData(MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in InitializeGadgetEquilibriumCoolData.\n");
      return FAIL;
    } 
  }

  /* If set, initialize the RadiativeCooling and RateEquations data. */
 
  if (MultiSpecies > 0)
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in InitializeRateData.\n");
      return FAIL;
    }
 
  if (MultiSpecies             == 0 && 
      RadiativeCooling          > 0 && 
      GadgetEquilibriumCooling == 0) {
    if (InitializeEquilibriumCoolData(MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in InitializeEquilibriumCoolData.\n");
      return FAIL;
    }
  }
 
  /* If using the internal radiation field, initialize it. */
 
  if (RadiationFieldType >= 10 && RadiationFieldType <= 11)
    if (InitializeRadiationFieldData(MetaData.Time) == FAIL) {
	fprintf(stderr, "Error in InitializeRadiationFieldData.\n");
	return FAIL;
      }
 
  /* Turn off DualEnergyFormalism for zeus hydro (and a few other things). */
 
  if (HydroMethod == Zeus_Hydro) {
    ConservativeInterpolation = FALSE;
    DualEnergyFormalism       = FALSE;
    //    FluxCorrection            = FALSE;
  }
 
  /* Set the number of particle attributes, if left unset. */
 
  if (NumberOfParticleAttributes == INT_UNDEFINED)
    if (StarParticleCreation || StarParticleFeedback)
      NumberOfParticleAttributes = 3;
    else
      NumberOfParticleAttributes = 0;
 
#ifdef UNUSED
  if (MaximumGravityRefinementLevel == INT_UNDEFINED)
    MaximumGravityRefinementLevel = (RadiativeCooling && SelfGravity
				     && HydroMethod == Zeus_Hydro) ?
       max(MaximumRefinementLevel-2, 5) : MaximumRefinementLevel;
#else
  if (MaximumGravityRefinementLevel == INT_UNDEFINED)
    MaximumGravityRefinementLevel = MaximumRefinementLevel;
#endif
 
  MaximumGravityRefinementLevel =
    min(MaximumGravityRefinementLevel, MaximumRefinementLevel);
 
  /* Use the value in MaximumParticleRefinementLevel to set the smoothing
     radius for the particles, to be used to Grid_DepositPositions. */
 
  if (MaximumParticleRefinementLevel >= 0)
    DepositPositionsParticleSmoothRadius =
      (DomainRightEdge[0] - DomainLeftEdge[0])/
      (float(MetaData.TopGridDims[0])*
       POW(float(RefineBy), float(MaximumParticleRefinementLevel)));
 
//  PPMDiffusion causes an out-of-bounds condition as currently written
//  The following is an over-ride to force PPMDiffusion OFF. This has
//  been fixed in this latest version (AK).

  if (MetaData.PPMDiffusionParameter != 0 && ProblemType != 60 // Turbulence
                                          && ProblemType != 4  // Double Mach Reflection test
                                          && ProblemType != 6  // Implosion test
                                          && ProblemType != 7  // SedovBlast test
                                          && ProblemType != 8  // KH test
                                          && ProblemType != 9  // Noh test
                                          ) {
    printf("WARNING! Setting MetaData.PPMDiffusionParameter = 0\n");
    MetaData.PPMDiffusionParameter = 0;
  }
 
  if (PartitionNestedGrids == 1) {
    printf("WARNING! PartitionNestedGrids = 1 forces Parallel IO = 1\n");
    ParallelRootGridIO = 1;
    ParallelParticleIO = 1;
  }

  if (TracerParticleOn) {
    ParticleTypeInFile = TRUE;
  }
 
  if (WritePotential && ComovingCoordinates && SelfGravity) {
    CopyGravPotential = TRUE;
  }

  if (MyProcessorNumber == ROOT_PROCESSOR) {
 
    if ( MetaData.GlobalDir != NULL ) {
      fprintf(stderr, "Output to Global Dir %s\n", MetaData.GlobalDir);
    }
 
    if ( MetaData.LocalDir != NULL ) {
      fprintf(stderr, "Output to Local Dir %s\n", MetaData.LocalDir);
    }

  }
 
  if ( (MetaData.GlobalDir != NULL) && (MetaData.LocalDir != NULL) ) {
    fprintf(stderr, "Cannot select GlobalDir AND LocalDir!\n");
    return FAIL;
  }
 
  char *cwd_buffer = new char[MAX_LINE_LENGTH];
  size_t cwd_buffer_len = MAX_LINE_LENGTH;
 
  if ( (MetaData.GlobalDir == NULL) && (MetaData.LocalDir == NULL) ) {
    if(getcwd(cwd_buffer, cwd_buffer_len) == NULL) {
      fprintf(stderr, "GETCWD call FAILED\n");
    }
    fprintf(stderr,"CWD %s\n", cwd_buffer);
    MetaData.GlobalDir = cwd_buffer;
    fprintf(stderr,"Global Dir set to %s\n", cwd_buffer);
  }
 
  char *srb_prefix_buffer = new char[MAX_LINE_LENGTH];
 
  if ( SRBprefix == NULL ) {
    strcpy(srb_prefix_buffer, "/NONE");
    SRBprefix = srb_prefix_buffer;
  }
 
  return SUCCESS;
}

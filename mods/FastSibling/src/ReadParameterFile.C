/***********************************************************************
/
/  READ A PARAMETER FILE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
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

int ParticleTypeInFile = 0;

/* function prototypes */

int ReadListOfFloats(FILE *fptr, int N, float floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
int CosmologyReadParameters(FILE *fptr, FLOAT *StopTime, FLOAT *InitTime);
int InitializeRateData(FLOAT Time);
int InitializeEquilibriumCoolData(FLOAT Time);
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

    ret += sscanf(line, "InitialCycleNumber = %d", &MetaData.CycleNumber);
    ret += sscanf(line, "InitialTime        = %"FSYM, &MetaData.Time);
    ret += sscanf(line, "InitialCPUTime     = %f", &MetaData.CPUTime);
    ret += sscanf(line, "Initialdt          = %f", Initialdt);

    ret += sscanf(line, "StopTime    = %"FSYM, &MetaData.StopTime);
    ret += sscanf(line, "StopCycle   = %d", &MetaData.StopCycle);
    ret += sscanf(line, "StopCPUTime = %f", &MetaData.StopCPUTime);

    ret += sscanf(line, "TimeLastRestartDump = %"FSYM, 
		  &MetaData.TimeLastRestartDump);
    ret += sscanf(line, "dtRestartDump       = %"FSYM, &MetaData.dtRestartDump);
    ret += sscanf(line, "TimeLastDataDump    = %"FSYM, 
		  &MetaData.TimeLastDataDump);
    ret += sscanf(line, "dtDataDump          = %"FSYM, &MetaData.dtDataDump);
    ret += sscanf(line, "TimeLastHistoryDump = %"FSYM, 
		  &MetaData.TimeLastHistoryDump);
    ret += sscanf(line, "dtHistoryDump       = %"FSYM, &MetaData.dtHistoryDump);
    ret += sscanf(line, "TimeLastMovieDump = %"FSYM, 
		  &MetaData.TimeLastMovieDump);
    ret += sscanf(line, "dtMovieDump       = %"FSYM, &MetaData.dtMovieDump);

    ret += sscanf(line, "TracerParticleOn  = %d", &TracerParticleOn);
    ret += sscanf(line, "TimeLastTracerParticleDump = %"FSYM, 
		  &MetaData.TimeLastTracerParticleDump);
    ret += sscanf(line, "dtTracerParticleDump       = %"FSYM, 
		  &MetaData.dtTracerParticleDump);

    ret += sscanf(line, "MovieRegionLeftEdge  = %"FSYM" %"FSYM" %"FSYM, 
		  MetaData.MovieRegionLeftEdge,
		  MetaData.MovieRegionLeftEdge+1, 
		  MetaData.MovieRegionLeftEdge+2);
    ret += sscanf(line, "MovieRegionRightEdge = %"FSYM" %"FSYM" %"FSYM, 
		  MetaData.MovieRegionRightEdge, 
		  MetaData.MovieRegionRightEdge+1,
		  MetaData.MovieRegionRightEdge+2);

    ret += sscanf(line, "CycleLastRestartDump = %d", 
		  &MetaData.CycleLastRestartDump);
    ret += sscanf(line, "CycleSkipRestartDump = %d", 
		  &MetaData.CycleSkipRestartDump);
    ret += sscanf(line, "CycleLastDataDump    = %d", 
		  &MetaData.CycleLastDataDump);
    ret += sscanf(line, "CycleSkipDataDump    = %d", 
		  &MetaData.CycleSkipDataDump);
    ret += sscanf(line, "CycleLastHistoryDump = %d", 
		  &MetaData.CycleLastHistoryDump);
    ret += sscanf(line, "CycleSkipHistoryDump = %d", 
		  &MetaData.CycleSkipHistoryDump);
    ret += sscanf(line, "OutputFirstTimeAtLevel = %d", 
		  &MetaData.OutputFirstTimeAtLevel);
    ret += sscanf(line, "StopFirstTimeAtLevel = %d", 
		  &MetaData.StopFirstTimeAtLevel);

    ret += sscanf(line, "RestartDumpNumber = %d", &MetaData.RestartDumpNumber);
    ret += sscanf(line, "DataDumpNumber    = %d", &MetaData.DataDumpNumber);
    ret += sscanf(line, "HistoryDumpNumber = %d", &MetaData.HistoryDumpNumber);
    ret += sscanf(line, "MovieDumpNumber   = %d", &MetaData.MovieDumpNumber);
    ret += sscanf(line, "TracerParticleDumpNumber = %d", 
		  &MetaData.TracerParticleDumpNumber);

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

    if (sscanf(line, "TimeActionType[%d] = %d", &dim, &int_dummy) == 2) {
      ret++; TimeActionType[dim] = int_dummy;
      if (dim >= MAX_TIME_ACTIONS-1) {
	fprintf(stderr, "Time action %d > maximum allowed.\n", dim);
	return FAIL;
      }
    }
    if (sscanf(line, "TimeActionRedshift[%d] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionRedshift[%d] = %"FSYM, &dim, 
		    TimeActionRedshift+dim);
    if (sscanf(line, "TimeActionTime[%d] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionTime[%d] = %"FSYM, &dim, 
		    TimeActionTime+dim);
    if (sscanf(line, "TimeActionParameter[%d] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionParameter[%d] = %f", &dim, 
		    TimeActionParameter+dim);
    
    ret += sscanf(line, "StaticHierarchy = %d", &MetaData.StaticHierarchy);

    ret += sscanf(line, "TopGridRank       = %d", &MetaData.TopGridRank);
    ret += sscanf(line, "TopGridDimensions = %d %d %d", MetaData.TopGridDims, 
		  MetaData.TopGridDims+1, MetaData.TopGridDims+2);

    ret += sscanf(line, "TopGridGravityBoundary = %d", 
		  &MetaData.GravityBoundary);

    ret += sscanf(line, "ParticleBoundaryType   = %d", 
		  &MetaData.ParticleBoundaryType);
    ret += sscanf(line, "NumberOfParticles      = %d", 
		  &MetaData.NumberOfParticles);

    ret += sscanf(line, "CourantSafetyNumber    = %f",
		  &MetaData.CourantSafetyNumber);
    ret += sscanf(line, "PPMFlatteningParameter = %d", 
		  &MetaData.PPMFlatteningParameter);
    ret += sscanf(line, "PPMDiffusionParameter  = %d", 
		  &MetaData.PPMDiffusionParameter);
    ret += sscanf(line, "PPMSteepeningParameter = %d", 
		  &MetaData.PPMSteepeningParameter);

    /* read global Parameters */

    ret += sscanf(line, "ProblemType            = %d", &ProblemType);
    ret += sscanf(line, "HydroMethod            = %d", &HydroMethod);
    ret += sscanf(line, "huge_number            = %f", &huge_number);
    ret += sscanf(line, "tiny_number            = %f", &tiny_number);
    ret += sscanf(line, "Gamma                  = %f", &Gamma);
    ret += sscanf(line, "PressureFree           = %d", &PressureFree);
    ret += sscanf(line, "RefineBy               = %d", &RefineBy);
    ret += sscanf(line, "MaximumRefinementLevel = %d", 
		  &MaximumRefinementLevel);
    ret += sscanf(line, "MaximumGravityRefinementLevel = %d", 
		  &MaximumGravityRefinementLevel);
    ret += sscanf(line, "MaximumParticleRefinementLevel = %d", 
		  &MaximumParticleRefinementLevel);
    ret += sscanf(line, "CellFlaggingMethod     = %d %d %d %d %d", 
	     CellFlaggingMethod+0, CellFlaggingMethod+1, CellFlaggingMethod+2,
	     CellFlaggingMethod+3, CellFlaggingMethod+4);
    ret += sscanf(line, "FluxCorrection         = %d", &FluxCorrection);
    ret += sscanf(line, "InterpolationMethod    = %d", &InterpolationMethod);
    ret += sscanf(line, "ConservativeInterpolation = %d", 
		  &ConservativeInterpolation);
    ret += sscanf(line, "MinimumEfficiency      = %f", &MinimumEfficiency);
    ret += sscanf(line, "NumberOfBufferZones    = %d", &NumberOfBufferZones);

    ret += sscanf(line, "DomainLeftEdge        = %"FSYM" %"FSYM" %"FSYM, DomainLeftEdge,
		  DomainLeftEdge+1, DomainLeftEdge+2);
    ret += sscanf(line, "DomainRightEdge       = %"FSYM" %"FSYM" %"FSYM, DomainRightEdge,
		  DomainRightEdge+1, DomainRightEdge+2);
    ret += sscanf(line, "GridVelocity          = %f %f %f", GridVelocity,
		  GridVelocity+1, GridVelocity+2);
    ret += sscanf(line, "RefineRegionLeftEdge  = %"FSYM" %"FSYM" %"FSYM, 
		  RefineRegionLeftEdge, RefineRegionLeftEdge+1, 
		  RefineRegionLeftEdge+2);
    ret += sscanf(line, "RefineRegionRightEdge = %"FSYM" %"FSYM" %"FSYM, 
		  RefineRegionRightEdge, RefineRegionRightEdge+1,
		  RefineRegionRightEdge+2);

    if (sscanf(line, "DataLabel[%d] = %s\n", &dim, dummy) == 2)
      DataLabel[dim] = dummy;
    if (sscanf(line, "DataUnits[%d] = %s\n", &dim, dummy) == 2)
      DataUnits[dim] = dummy;

    ret += sscanf(line, "UniformGravity          = %d", &UniformGravity);
    ret += sscanf(line, "UniformGravityDirection = %d", 
		  &UniformGravityDirection);
    ret += sscanf(line, "UniformGravityConstant  = %f", 
		  &UniformGravityConstant);

    ret += sscanf(line, "PointSourceGravity         = %d",&PointSourceGravity);
    ret += sscanf(line, "PointSourceGravityPosition = %"FSYM" %"FSYM" %"FSYM, 
		  PointSourceGravityPosition, PointSourceGravityPosition+1, 
		  PointSourceGravityPosition+2);
    ret += sscanf(line, "PointSourceGravityConstant = %lf", 
		  &PointSourceGravityConstant);
    ret += sscanf(line, "PointSourceGravityCoreRadius = %lf", 
		  &PointSourceGravityCoreRadius);

    ret += sscanf(line, "SelfGravity           = %d", &SelfGravity);
    ret += sscanf(line, "GravitationalConstant = %f", &GravitationalConstant);
    ret += sscanf(line, "S2ParticleSize        = %f", &S2ParticleSize);
    ret += sscanf(line, "GravityResolution     = %f", &GravityResolution);
    ret += sscanf(line, "ComputePotential      = %d", &ComputePotential);
    ret += sscanf(line, "BaryonSelfGravityApproximation = %d",
		  &BaryonSelfGravityApproximation);

    ret += sscanf(line, "GreensFunctionMaxNumber   = %d", 
		  &GreensFunctionMaxNumber);
    ret += sscanf(line, "GreensFunctionMaxSize     = %d",
		  &GreensFunctionMaxSize);

    ret += sscanf(line, "DualEnergyFormalism     = %d", &DualEnergyFormalism);
    ret += sscanf(line, "DualEnergyFormalismEta1 = %f", 
		  &DualEnergyFormalismEta1);
    ret += sscanf(line, "DualEnergyFormalismEta2 = %f", 
		  &DualEnergyFormalismEta2);
    ret += sscanf(line, "ParticleCourantSafetyNumber = %f",
		  &ParticleCourantSafetyNumber);
    ret += sscanf(line, "RadiativeCooling = %d", &RadiativeCooling);
    ret += sscanf(line, "MultiSpecies = %d", &MultiSpecies);
    ret += sscanf(line, "RadiationFieldType = %d", &RadiationFieldType);
    ret += sscanf(line, "RadiationFieldLevelRecompute = %d", 
		  &RadiationFieldLevelRecompute);
    ret += sscanf(line, "RadiationSpectrumNormalization = %f",
		  &CoolData.f3);
    ret += sscanf(line, "RadiationSpectrumSlope = %f", &CoolData.alpha0);

    ret += sscanf(line, "ZEUSQuadraticArtificialViscosity = %f",
		  &ZEUSQuadraticArtificialViscosity);
    ret += sscanf(line, "ZEUSLinearArtificialViscosity = %f",
		  &ZEUSLinearArtificialViscosity);

    ret += sscanf(line, "UseMinimumPressureSupport = %d",
		  &UseMinimumPressureSupport);
    ret += sscanf(line, "MinimumPressureSupportParameter = %f",
		  &MinimumPressureSupportParameter);
    ret += sscanf(line, "RefineByJeansLengthSafetyFactor = %f",
		  &RefineByJeansLengthSafetyFactor);
    ret += sscanf(line, "MustRefineParticlesRefineToLevel = %d",
		  &MustRefineParticlesRefineToLevel);
    ret += sscanf(line, "ParticleTypeInFile = %d",
		  &ParticleTypeInFile);
    
    if (sscanf(line, "StaticRefineRegionLevel[%d] = %d",&dim,&int_dummy) == 2){
      if (dim > MAX_STATIC_REGIONS-1) {
	fprintf(stderr, "StaticRegion number %d > MAX allowed\n", dim);
	return FAIL;
      }
      ret++;
      StaticRefineRegionLevel[dim] = int_dummy;
    }
    if (sscanf(line, "StaticRefineRegionLeftEdge[%d] = ", &dim) == 1)
      ret += sscanf(line, 
		    "StaticRefineRegionLeftEdge[%d] = %"FSYM" %"FSYM" %"FSYM,
		    &dim, StaticRefineRegionLeftEdge[dim],
		    StaticRefineRegionLeftEdge[dim]+1,
		    StaticRefineRegionLeftEdge[dim]+2);
    if (sscanf(line, "StaticRefineRegionRightEdge[%d] = ", &dim) == 1)
      ret += sscanf(line, 
		    "StaticRefineRegionRightEdge[%d] = %"FSYM" %"FSYM" %"FSYM,
		    &dim, StaticRefineRegionRightEdge[dim],
		    StaticRefineRegionRightEdge[dim]+1,
		    StaticRefineRegionRightEdge[dim]+2);

    ret += sscanf(line, "ParallelRootGridIO = %d", &ParallelRootGridIO);

    ret += sscanf(line, "MinimumOverDensityForRefinement  = %f %f %f %f %f", 
	  MinimumOverDensityForRefinement+0, MinimumOverDensityForRefinement+1,
	  MinimumOverDensityForRefinement+2, MinimumOverDensityForRefinement+3,
	  MinimumOverDensityForRefinement+4);
    ret += sscanf(line, "MinimumMassForRefinement  = %f %f %f %f %f", 
	  MinimumMassForRefinement+0, MinimumMassForRefinement+1,
	  MinimumMassForRefinement+2, MinimumMassForRefinement+3,
	  MinimumMassForRefinement+4);
    ret += sscanf(line, "MinimumMassForRefinementLevelExponent = %f %f %f %f %f", 
		  MinimumMassForRefinementLevelExponent+0, 
		  MinimumMassForRefinementLevelExponent+1,
		  MinimumMassForRefinementLevelExponent+2, 
		  MinimumMassForRefinementLevelExponent+3,
		  MinimumMassForRefinementLevelExponent+4);
    ret += sscanf(line, "MinimumSlopeForRefinement = %f", 
		  &MinimumSlopeForRefinement);
    ret += sscanf(line, "MinimumPressureJumpForRefinement = %f", 
		  &MinimumPressureJumpForRefinement);
    ret += sscanf(line, "MinimumEnergyRatioForRefinement = %f", 
		  &MinimumEnergyRatioForRefinement);
    ret += sscanf(line, "ComovingCoordinates = %d",&ComovingCoordinates);
    ret += sscanf(line, "StarParticleCreation = %d", &StarParticleCreation);
    ret += sscanf(line, "StarParticleFeedback = %d", &StarParticleFeedback);
    ret += sscanf(line, "NumberOfParticleAttributes = %d", 
		  &NumberOfParticleAttributes);

    /* read data which defines the boundary conditions */

    ret += sscanf(line, "LeftFaceBoundaryCondition  = %d %d %d", 
		  MetaData.LeftFaceBoundaryCondition,
		  MetaData.LeftFaceBoundaryCondition+1,
		  MetaData.LeftFaceBoundaryCondition+2);
    ret += sscanf(line, "RightFaceBoundaryCondition = %d %d %d",
		  MetaData.RightFaceBoundaryCondition,
		  MetaData.RightFaceBoundaryCondition+1,
		  MetaData.RightFaceBoundaryCondition+2);
    if (sscanf(line, "BoundaryConditionName         = %s", dummy) == 1)
      MetaData.BoundaryConditionName = dummy;

    /* Check version number. */

    if (sscanf(line, "VersionNumber = %f", &TempFloat) == 1) {
      ret++;
      if (fabs(TempFloat - VERSION) >= 1.0e-3)
	fprintf(stderr, "Warning: Incorrect version number.\n");
    }

    /* Read star particle parameters. */

    ret += sscanf(line, "StarMakerOverDensityThreshold = %f",
		  &StarMakerOverDensityThreshold);
    ret += sscanf(line, "StarMakerMassEfficiency = %f",
		  &StarMakerMassEfficiency);
    ret += sscanf(line, "StarMakerMinimumMass = %f", &StarMakerMinimumMass);
    ret += sscanf(line, "StarMakerMinimumDynamicalTime = %f",
                  &StarMakerMinimumDynamicalTime);
    ret += sscanf(line, "StarMassEjectionFraction = %f", 
		  &StarMassEjectionFraction);
    ret += sscanf(line, "StarMetalYield = %f", &StarMetalYield);
    ret += sscanf(line, "StarEnergyToThermalFeedback = %f", 
		  &StarEnergyToThermalFeedback);
    ret += sscanf(line, "StarEnergyToStellarUV = %f", &StarEnergyToStellarUV);
    ret += sscanf(line, "StarEnergyToQuasarUV = %f", &StarEnergyToQuasarUV);

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
    if (strstr(line, "ZeldovichPancake")    ) ret++;
    if (strstr(line, "PressurelessCollapse")) ret++;
    if (strstr(line, "AdiabaticExpansion")  ) ret++;
    if (strstr(line, "CosmologySimulation") ) ret++;
    if (strstr(line, "TestGravity"        ) ) ret++;
    if (strstr(line, "SphericalInfall"    ) ) ret++;
    if (strstr(line, "TestGravitySphere"  ) ) ret++;
    if (strstr(line, "CollapseTest"       ) ) ret++;
    if (strstr(line, "Cosmology")           ) ret++;
    if (strstr(line, "SupernovaRestart")    ) ret++;
    if (strstr(line, "TracerParticleCreation")) ret++;
    if (strstr(line, "GalaxySimulation")    ) ret++;

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

  /* If set, initialize the RadiativeCooling and RateEquations data. */

  if (MultiSpecies > 0)
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in InitializeRateData.\n");
      return FAIL;
    }

  if (MultiSpecies == 0 && RadiativeCooling > 0)
    if (InitializeEquilibriumCoolData(MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in InitializeEquilibriumCoolData.\n");
      return FAIL;
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
       pow(float(RefineBy), float(MaximumParticleRefinementLevel)));

  return SUCCESS;
}

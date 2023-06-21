/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  WRITES A PARAMETER FILE (i.e. TopGrid data)
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

// This routine writes the parameter file in the argument and sets parameters
//   based on it.

#include <stdio.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "StarParticleData.h"

/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int  CosmologyWriteParameters(FILE *fptr, FLOAT StopTime, FLOAT CurrentTime);
int  CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		       float *TemperatureUnits, float *TimeUnits,
		       float *VelocityUnits, FLOAT Time);

int WriteParameterFile(FILE *fptr, TopGridData &MetaData)
{

  int dim;

  /* Compute Units. */

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
        VelocityUnits = 1;

  if (ComovingCoordinates)
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }

  /* write data to Parameter output file */

  /* write MetaData parameters */

  fprintf(fptr, "InitialCycleNumber  = %d\n", MetaData.CycleNumber);
  fprintf(fptr, "InitialTime         = %"GOUTSYM"\n", MetaData.Time);
  fprintf(fptr, "InitialCPUTime      = %g\n\n", MetaData.CPUTime);

  fprintf(fptr, "StopTime            = %"GOUTSYM"\n", MetaData.StopTime);
  fprintf(fptr, "StopCycle           = %d\n", MetaData.StopCycle);
  fprintf(fptr, "StopCPUTime         = %g\n\n", MetaData.StopCPUTime);

  fprintf(fptr, "TimeLastRestartDump = %"GOUTSYM"\n", MetaData.TimeLastRestartDump);
  fprintf(fptr, "dtRestartDump       = %"GOUTSYM"\n", MetaData.dtRestartDump);
  fprintf(fptr, "TimeLastDataDump    = %"GOUTSYM"\n", MetaData.TimeLastDataDump);
  fprintf(fptr, "dtDataDump          = %"GOUTSYM"\n", MetaData.dtDataDump);
  fprintf(fptr, "TimeLastHistoryDump = %"GOUTSYM"\n", MetaData.TimeLastHistoryDump);
  fprintf(fptr, "dtHistoryDump       = %"GOUTSYM"\n\n", MetaData.dtHistoryDump);
  fprintf(fptr, "TimeLastMovieDump     = %"GOUTSYM"\n", MetaData.TimeLastMovieDump);
  fprintf(fptr, "dtMovieDump           = %"GOUTSYM"\n", MetaData.dtMovieDump);

  fprintf(fptr, "MovieRegionLeftEdge   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.MovieRegionLeftEdge);
  fprintf(fptr, "MovieRegionRightEdge  = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.MovieRegionRightEdge);
  fprintf(fptr, "\n");

  fprintf(fptr, "CycleLastRestartDump = %d\n", MetaData.CycleLastRestartDump);
  fprintf(fptr, "CycleSkipRestartDump = %d\n", MetaData.CycleSkipRestartDump);
  fprintf(fptr, "CycleLastDataDump    = %d\n", MetaData.CycleLastDataDump);
  fprintf(fptr, "CycleSkipDataDump    = %d\n", MetaData.CycleSkipDataDump);
  fprintf(fptr, "CycleLastHistoryDump = %d\n", MetaData.CycleLastHistoryDump);
  fprintf(fptr, "CycleSkipHistoryDump = %d\n\n", 
	  MetaData.CycleSkipHistoryDump);

  fprintf(fptr, "OutputFirstTimeAtLevel = %d\n", 
	  MetaData.OutputFirstTimeAtLevel);
  fprintf(fptr, "StopFirstTimeAtLevel = %d\n\n", 
	  MetaData.StopFirstTimeAtLevel);

  fprintf(fptr, "RestartDumpNumber   = %d\n", MetaData.RestartDumpNumber);
  fprintf(fptr, "DataDumpNumber      = %d\n", MetaData.DataDumpNumber);
  fprintf(fptr, "HistoryDumpNumber   = %d\n", MetaData.HistoryDumpNumber);
  fprintf(fptr, "MovieDumpNumber     = %d\n", MetaData.MovieDumpNumber);

  fprintf(fptr, "RestartDumpName     = %s\n", MetaData.RestartDumpName);
  fprintf(fptr, "DataDumpName        = %s\n", MetaData.DataDumpName);
  fprintf(fptr, "HistoryDumpName     = %s\n", MetaData.HistoryDumpName);
  fprintf(fptr, "MovieDumpName       = %s\n", MetaData.MovieDumpName);
  fprintf(fptr, "RedshiftDumpName    = %s\n\n", MetaData.RedshiftDumpName);

  for (dim = 0; dim < MAX_TIME_ACTIONS; dim++)
    if (TimeActionType[dim] > 0) {
      fprintf(fptr, "TimeActionType[%d]      = %d\n", dim,TimeActionType[dim]);
      if (ComovingCoordinates)
	fprintf(fptr, "TimeActionRedshift[%d]  = %"GOUTSYM"\n", dim, 
		TimeActionRedshift[dim]);
      else
	fprintf(fptr, "TimeActionTime[%d]      = %"GOUTSYM"\n", dim, 
		TimeActionRedshift[dim]);
      fprintf(fptr, "TimeActionParameter[%d] = %g\n", dim,
	      TimeActionParameter[dim]);
    }

  fprintf(fptr, "StaticHierarchy     = %d\n", MetaData.StaticHierarchy);

  fprintf(fptr, "TopGridRank         = %d\n", MetaData.TopGridRank);
  fprintf(fptr, "TopGridDimensions   = ");
  WriteListOfInts(fptr, MetaData.TopGridRank, MetaData.TopGridDims);
  fprintf(fptr, "\n");

  fprintf(fptr, "TopGridGravityBoundary = %d\n", MetaData.GravityBoundary);

  fprintf(fptr, "ParticleBoundaryType   = %d\n",MetaData.ParticleBoundaryType);
  fprintf(fptr, "NumberOfParticles      = %d (do not modify)\n", 
	  MetaData.NumberOfParticles);

  fprintf(fptr, "CourantSafetyNumber    = %g\n",
	  MetaData.CourantSafetyNumber);
  fprintf(fptr, "PPMFlatteningParameter = %d\n", 
	  MetaData.PPMFlatteningParameter);
  fprintf(fptr, "PPMDiffusionParameter  = %d\n", 
	  MetaData.PPMDiffusionParameter);
  fprintf(fptr, "PPMSteepeningParameter = %d\n\n", 
	  MetaData.PPMSteepeningParameter);

  /* write global Parameters */

  fprintf(fptr, "ProblemType            = %d\n", ProblemType);
  fprintf(fptr, "HydroMethod            = %d\n", HydroMethod);
  fprintf(fptr, "huge_number            = %e\n", huge_number);
  fprintf(fptr, "tiny_number            = %e\n", tiny_number);
  fprintf(fptr, "Gamma                  = %g\n", Gamma);
  fprintf(fptr, "PressureFree           = %d\n", PressureFree);
  fprintf(fptr, "RefineBy               = %d\n", RefineBy);
  fprintf(fptr, "MaximumRefinementLevel = %d\n", MaximumRefinementLevel);
  fprintf(fptr, "MaximumGravityRefinementLevel = %d\n", 
	  MaximumGravityRefinementLevel);
  fprintf(fptr, "MaximumParticleRefinementLevel = %d\n", 
	  MaximumParticleRefinementLevel);
  fprintf(fptr, "CellFlaggingMethod     = ");
  WriteListOfInts(fptr, MAX_FLAGGING_METHODS, CellFlaggingMethod);
  fprintf(fptr, "FluxCorrection         = %d\n", FluxCorrection);
  fprintf(fptr, "InterpolationMethod    = %d\n", InterpolationMethod);
  fprintf(fptr, "ConservativeInterpolation = %d\n", ConservativeInterpolation);
  fprintf(fptr, "MinimumEfficiency      = %g\n", MinimumEfficiency);
  fprintf(fptr, "NumberOfBufferZones    = %d\n\n", NumberOfBufferZones);

  fprintf(fptr, "DomainLeftEdge         = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, DomainLeftEdge);
  fprintf(fptr, "DomainRightEdge        = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, DomainRightEdge);
  fprintf(fptr, "GridVelocity           = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, GridVelocity);
  fprintf(fptr, "RefineRegionLeftEdge   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, RefineRegionLeftEdge);
  fprintf(fptr, "RefineRegionRightEdge  = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, RefineRegionRightEdge);
  fprintf(fptr, "\n");

  for (dim = 0; dim < MAX_NUMBER_OF_BARYON_FIELDS; dim++) {
    if (DataLabel[dim])
      fprintf(fptr, "DataLabel[%d]              = %s\n", dim, DataLabel[dim]);
    if (DataUnits[dim])
      fprintf(fptr, "DataUnits[%d]              = %s\n", dim, DataUnits[dim]);
    if (DataLabel[dim]) {
      if (strstr(DataLabel[dim], "Density") != NULL)
	fprintf(fptr, "#DataCGSConversionFactor[%d] = %g\n", dim,DensityUnits);
      if (strstr(DataLabel[dim], "velocity") != NULL)
	fprintf(fptr, "#DataCGSConversionFactor[%d] = %g\n",dim,VelocityUnits);
    }
  }
  fprintf(fptr, "\n");

  fprintf(fptr, "UniformGravity             = %d\n", UniformGravity);
  fprintf(fptr, "UniformGravityDirection    = %d\n", UniformGravityDirection);
  fprintf(fptr, "UniformGravityConstant     = %g\n", UniformGravityConstant);

  fprintf(fptr, "PointSourceGravity           = %d\n",PointSourceGravity);
  fprintf(fptr, "PointSourceGravityPosition   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, PointSourceGravityPosition);
  fprintf(fptr, "PointSourceGravityConstant   = %g\n", 
	  PointSourceGravityConstant);
  fprintf(fptr, "PointSourceGravityCoreRadius = %g\n\n", 
	  PointSourceGravityCoreRadius);

  fprintf(fptr, "SelfGravity                    = %d\n", SelfGravity);
  fprintf(fptr, "GravitationalConstant          = %e\n", 
	  GravitationalConstant);
  fprintf(fptr, "S2ParticleSize                 = %g\n", S2ParticleSize);
  fprintf(fptr, "GravityResolution              = %g\n", GravityResolution);
  fprintf(fptr, "ComputePotential               = %d\n", ComputePotential);
  fprintf(fptr, "BaryonSelfGravityApproximation = %d\n\n",
	  BaryonSelfGravityApproximation);

  fprintf(fptr, "GreensFunctionMaxNumber     = %d\n", GreensFunctionMaxNumber);
  fprintf(fptr, "GreensFunctionMaxSize       = %d\n", GreensFunctionMaxSize);

  fprintf(fptr, "DualEnergyFormalism         = %d\n", DualEnergyFormalism);
  fprintf(fptr, "DualEnergyFormalismEta1     = %e\n", DualEnergyFormalismEta1);
  fprintf(fptr, "DualEnergyFormalismEta2     = %e\n", DualEnergyFormalismEta2);
  fprintf(fptr, "ParticleCourantSafetyNumber = %g\n\n", 
	  ParticleCourantSafetyNumber);
  fprintf(fptr, "RadiativeCooling               = %d\n", RadiativeCooling);
  fprintf(fptr, "GadgetEquilibriumCooling       = %d\n", 
	  GadgetEquilibriumCooling);
  fprintf(fptr, "MultiSpecies                   = %d\n", MultiSpecies);
  fprintf(fptr, "RadiationFieldType             = %d\n", RadiationFieldType);
  fprintf(fptr, "RadiationFieldLevelRecompute   = %d\n", 
	  RadiationFieldLevelRecompute);
  fprintf(fptr, "RadiationSpectrumNormalization = %g\n", CoolData.f3);
  fprintf(fptr, "RadiationSpectrumSlope         = %g\n", CoolData.alpha0);

  fprintf(fptr, "ZEUSLinearArtificialViscosity    = %g\n", 
	  ZEUSLinearArtificialViscosity);
  fprintf(fptr, "ZEUSQuadraticArtificialViscosity = %g\n", 
	  ZEUSQuadraticArtificialViscosity);
  fprintf(fptr, "UseMinimumPressureSupport        = %d\n", 
	  UseMinimumPressureSupport);
  fprintf(fptr, "MinimumPressureSupportParameter  = %f\n", 
	  MinimumPressureSupportParameter);
  fprintf(fptr, "RefineByJeansLengthSafetyFactor  = %f\n\n", 
	  RefineByJeansLengthSafetyFactor);

  for (dim = 0; dim < MAX_STATIC_REGIONS; dim++)
    if (StaticRefineRegionLevel[dim] != INT_UNDEFINED) {
      fprintf(fptr, "StaticRefineRegionLevel[%d] = %d\n", dim,
	      StaticRefineRegionLevel[dim]);
      fprintf(fptr, "StaticRefineRegionLeftEdge[%d] = ", dim);
      WriteListOfFloats(fptr, MAX_DIMENSION, StaticRefineRegionLeftEdge[dim]);
      fprintf(fptr, "StaticRefineRegionRightEdge[%d] = ", dim);
      WriteListOfFloats(fptr, MAX_DIMENSION, StaticRefineRegionRightEdge[dim]);
    }

  fprintf(fptr, "ParallelRootGridIO              = %d\n", ParallelRootGridIO);
  fprintf(fptr, "ParallelParticleIO              = %d\n", ParallelParticleIO);
  fprintf(fptr, "Unigrid                         = %d\n", Unigrid);
  fprintf(fptr, "ExtractFieldsOnly               = %d\n", ExtractFieldsOnly);

  fprintf(fptr, "MinimumOverDensityForRefinement = %g %g %g %g %g\n",
       MinimumOverDensityForRefinement[0], MinimumOverDensityForRefinement[1],
       MinimumOverDensityForRefinement[2], MinimumOverDensityForRefinement[3],
       MinimumOverDensityForRefinement[4]);
  fprintf(fptr, "MinimumMassForRefinement = %.9g %.9g %.9g %.9g %.9g\n",
       MinimumMassForRefinement[0], MinimumMassForRefinement[1],
       MinimumMassForRefinement[2], MinimumMassForRefinement[3],
       MinimumMassForRefinement[4]);
  fprintf(fptr, "MinimumMassForRefinementLevelExponent = %f %f %f %f %f\n",
	  MinimumMassForRefinementLevelExponent[0], 
	  MinimumMassForRefinementLevelExponent[1],
	  MinimumMassForRefinementLevelExponent[2], 
	  MinimumMassForRefinementLevelExponent[3],
	  MinimumMassForRefinementLevelExponent[4]);
  fprintf(fptr, "MinimumSlopeForRefinement             = %e\n",
	  MinimumSlopeForRefinement);
  fprintf(fptr, "MinimumPressureJumpForRefinement      = %e\n",
	  MinimumPressureJumpForRefinement);
  fprintf(fptr, "MinimumEnergyRatioForRefinement       = %e\n",
	  MinimumEnergyRatioForRefinement);
  fprintf(fptr, "ComovingCoordinates                   = %d\n", 
	  ComovingCoordinates);
  fprintf(fptr, "StarParticleCreation                  = %d\n", 
	  StarParticleCreation);
  fprintf(fptr, "StarParticleFeedback                  = %d\n", 
	  StarParticleFeedback);
  fprintf(fptr, "NumberOfParticleAttributes            = %d\n", 
	  NumberOfParticleAttributes);
  fprintf(fptr, "StarMakerOverDensityThreshold         = %g\n",
	  StarMakerOverDensityThreshold);
  fprintf(fptr, "StarMakerMassEfficiency               = %g\n",
	  StarMakerMassEfficiency);
  fprintf(fptr, "StarMakerMinimumMass                  = %g\n",
          StarMakerMinimumMass);
  fprintf(fptr, "StarMakerMinimumDynamicalTime         = %g\n",
          StarMakerMinimumDynamicalTime);
  fprintf(fptr, "StarMassEjectionFraction              = %g\n",
          StarMassEjectionFraction);
  fprintf(fptr, "StarMetalYield                        = %g\n",
          StarMetalYield);
  fprintf(fptr, "StarEnergyToThermalFeedback           = %g\n",
          StarEnergyToThermalFeedback);
  fprintf(fptr, "StarEnergyToStellarUV                 = %g\n",
          StarEnergyToStellarUV);
  fprintf(fptr, "StarEnergyToQuasarUV                  = %g\n\n",
          StarEnergyToQuasarUV);
  fprintf(fptr, "MultiMetals                           = %d\n",
          MultiMetals);

  /* write data which defines the boundary conditions */
  
  fprintf(fptr, "LeftFaceBoundaryCondition  = ");
  WriteListOfInts(fptr, MetaData.TopGridRank, 
		  (int*) MetaData.LeftFaceBoundaryCondition);
  fprintf(fptr, "RightFaceBoundaryCondition = ");
  WriteListOfInts(fptr, MetaData.TopGridRank, 
		  (int*) MetaData.RightFaceBoundaryCondition);
  if (MetaData.BoundaryConditionName)
    fprintf(fptr, "BoundaryConditionName      = %s\n\n", 
	    MetaData.BoundaryConditionName);

  /* If appropriate, write Cosmology data. */

  if (ComovingCoordinates)
    if (CosmologyWriteParameters(fptr, MetaData.StopTime, MetaData.Time) == 
	FAIL) {
      fprintf(stderr, "Error in CosmologyWriteParameters.\n");
      return FAIL;
    }

  /* write version info */

  fprintf(fptr, "VersionNumber              = %f\n\n", VERSION);

  return SUCCESS;
}

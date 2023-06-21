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
 
  fprintf(fptr, "InitialCycleNumber  = %"ISYM"\n", MetaData.CycleNumber);
  fprintf(fptr, "InitialTime         = %"GOUTSYM"\n", MetaData.Time);
  fprintf(fptr, "InitialCPUTime      = %"GSYM"\n\n", MetaData.CPUTime);
 
  fprintf(fptr, "StopTime            = %"GOUTSYM"\n", MetaData.StopTime);
  fprintf(fptr, "StopCycle           = %"ISYM"\n", MetaData.StopCycle);
  fprintf(fptr, "StopCPUTime         = %"GSYM"\n\n", MetaData.StopCPUTime);
 
  fprintf(fptr, "TimeLastRestartDump = %"GOUTSYM"\n", MetaData.TimeLastRestartDump);
  fprintf(fptr, "dtRestartDump       = %"GOUTSYM"\n", MetaData.dtRestartDump);
  fprintf(fptr, "TimeLastDataDump    = %"GOUTSYM"\n", MetaData.TimeLastDataDump);
  fprintf(fptr, "dtDataDump          = %"GOUTSYM"\n", MetaData.dtDataDump);
  fprintf(fptr, "TimeLastHistoryDump = %"GOUTSYM"\n", MetaData.TimeLastHistoryDump);
  fprintf(fptr, "dtHistoryDump       = %"GOUTSYM"\n\n", MetaData.dtHistoryDump);
  fprintf(fptr, "TimeLastMovieDump     = %"GOUTSYM"\n", MetaData.TimeLastMovieDump);
  fprintf(fptr, "dtMovieDump           = %"GOUTSYM"\n", MetaData.dtMovieDump);
 
  fprintf(fptr, "TracerParticleOn           = %"ISYM"\n", TracerParticleOn);
  fprintf(fptr, "TimeLastTracerParticleDump = %"GOUTSYM"\n",
          MetaData.TimeLastTracerParticleDump);
  fprintf(fptr, "dtTracerParticleDump       = %"GOUTSYM"\n",
          MetaData.dtTracerParticleDump);
 
  fprintf(fptr, "MovieRegionLeftEdge   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.MovieRegionLeftEdge);
  fprintf(fptr, "MovieRegionRightEdge  = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.MovieRegionRightEdge);
  fprintf(fptr, "\n");
 
  fprintf(fptr, "NewMovieLeftEdge   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.NewMovieLeftEdge);
  fprintf(fptr, "NewMovieRightEdge  = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.NewMovieRightEdge);
  fprintf(fptr, "MovieSkipTimestep = %"ISYM"\n", MovieSkipTimestep);
  fprintf(fptr, "NewMovieParticleOn = %"ISYM"\n", NewMovieParticleOn);
  fprintf(fptr, "MovieDataField = ");
  WriteListOfInts(fptr, MAX_MOVIE_FIELDS, MovieDataField);
  fprintf(fptr, "NewMovieDumpNumber = %"ISYM"\n", NewMovieDumpNumber);
  fprintf(fptr, "NewMovieName = %s\n", NewMovieName);
  fprintf(fptr, "\n");

  fprintf(fptr, "CycleLastRestartDump = %"ISYM"\n", MetaData.CycleLastRestartDump);
  fprintf(fptr, "CycleSkipRestartDump = %"ISYM"\n", MetaData.CycleSkipRestartDump);
  fprintf(fptr, "CycleLastDataDump    = %"ISYM"\n", MetaData.CycleLastDataDump);
  fprintf(fptr, "CycleSkipDataDump    = %"ISYM"\n", MetaData.CycleSkipDataDump);
  fprintf(fptr, "CycleLastHistoryDump = %"ISYM"\n", MetaData.CycleLastHistoryDump);
  fprintf(fptr, "CycleSkipHistoryDump = %"ISYM"\n\n",
	  MetaData.CycleSkipHistoryDump);
  fprintf(fptr, "CycleSkipGlobalDataDump = %"ISYM"\n\n", //AK
          MetaData.CycleSkipGlobalDataDump);
 
  fprintf(fptr, "OutputFirstTimeAtLevel = %"ISYM"\n",
	  MetaData.OutputFirstTimeAtLevel);
  fprintf(fptr, "StopFirstTimeAtLevel = %"ISYM"\n\n",
	  MetaData.StopFirstTimeAtLevel);
 
  fprintf(fptr, "RestartDumpNumber   = %"ISYM"\n", MetaData.RestartDumpNumber);
  fprintf(fptr, "DataDumpNumber      = %"ISYM"\n", MetaData.DataDumpNumber);
  fprintf(fptr, "HistoryDumpNumber   = %"ISYM"\n", MetaData.HistoryDumpNumber);
  fprintf(fptr, "MovieDumpNumber     = %"ISYM"\n", MetaData.MovieDumpNumber);
  fprintf(fptr, "TracerParticleDumpNumber = %"ISYM"\n",
          MetaData.TracerParticleDumpNumber);
 
  fprintf(fptr, "RestartDumpName     = %s\n", MetaData.RestartDumpName);
  fprintf(fptr, "DataDumpName        = %s\n", MetaData.DataDumpName);
  fprintf(fptr, "HistoryDumpName     = %s\n", MetaData.HistoryDumpName);
  fprintf(fptr, "MovieDumpName       = %s\n", MetaData.MovieDumpName);
  fprintf(fptr, "TracerParticleDumpName = %s\n",
          MetaData.TracerParticleDumpName);
  fprintf(fptr, "RedshiftDumpName    = %s\n\n", MetaData.RedshiftDumpName);
 
  if (MetaData.RestartDumpDir != NULL)
    fprintf(fptr, "RestartDumpDir      = %s\n", MetaData.RestartDumpDir);
  if (MetaData.DataDumpDir != NULL)
    fprintf(fptr, "DataDumpDir         = %s\n", MetaData.DataDumpDir);
  if (MetaData.HistoryDumpDir != NULL)
    fprintf(fptr, "HistoryDumpDir      = %s\n", MetaData.HistoryDumpDir);
  if (MetaData.MovieDumpDir != NULL)
    fprintf(fptr, "MovieDumpDir        = %s\n", MetaData.MovieDumpDir);
  if (MetaData.TracerParticleDumpDir != NULL)
    fprintf(fptr, "TracerParticleDumpDir = %s\n", MetaData.TracerParticleDumpDir);
  if (MetaData.RedshiftDumpDir != NULL)
    fprintf(fptr, "RedshiftDumpDir     = %s\n\n", MetaData.RedshiftDumpDir);
 
  if (MetaData.LocalDir != NULL)
    fprintf(fptr, "LocalDir            = %s\n", MetaData.LocalDir);
  if (MetaData.GlobalDir != NULL)
    fprintf(fptr, "GlobalDir           = %s\n", MetaData.GlobalDir);
 
  for (dim = 0; dim < MAX_CUBE_DUMPS; dim++)
    if (CubeDumps[dim] != NULL)
      fprintf(fptr, "CubeDump[%"ISYM"]            = %s\n", dim, CubeDumps[dim]);
 
  for (dim = 0; dim < MAX_TIME_ACTIONS; dim++)
    if (TimeActionType[dim] > 0) {
      fprintf(fptr, "TimeActionType[%"ISYM"]      = %"ISYM"\n", dim,TimeActionType[dim]);
      if (ComovingCoordinates)
	fprintf(fptr, "TimeActionRedshift[%"ISYM"]  = %"GOUTSYM"\n", dim,
		TimeActionRedshift[dim]);
      else
	fprintf(fptr, "TimeActionTime[%"ISYM"]      = %"GOUTSYM"\n", dim,
		TimeActionRedshift[dim]);
      fprintf(fptr, "TimeActionParameter[%"ISYM"] = %"GSYM"\n", dim,
	      TimeActionParameter[dim]);
    }
 
  fprintf(fptr, "StaticHierarchy     = %"ISYM"\n", MetaData.StaticHierarchy);
 
  fprintf(fptr, "TopGridRank         = %"ISYM"\n", MetaData.TopGridRank);
  fprintf(fptr, "TopGridDimensions   = ");
  WriteListOfInts(fptr, MetaData.TopGridRank, MetaData.TopGridDims);
  fprintf(fptr, "\n");
 
  fprintf(fptr, "TopGridGravityBoundary = %"ISYM"\n", MetaData.GravityBoundary);
 
  fprintf(fptr, "ParticleBoundaryType   = %"ISYM"\n",MetaData.ParticleBoundaryType);
  fprintf(fptr, "NumberOfParticles      = %"ISYM" (do not modify)\n",
	  MetaData.NumberOfParticles);
 
  fprintf(fptr, "CourantSafetyNumber    = %"GSYM"\n",
	  MetaData.CourantSafetyNumber);
  fprintf(fptr, "PPMFlatteningParameter = %"ISYM"\n",
	  MetaData.PPMFlatteningParameter);
  fprintf(fptr, "PPMDiffusionParameter  = %"ISYM"\n",
	  MetaData.PPMDiffusionParameter);
  fprintf(fptr, "PPMSteepeningParameter = %"ISYM"\n\n",
	  MetaData.PPMSteepeningParameter);
 
  /* write global Parameters */
 
  fprintf(fptr, "ProblemType            = %"ISYM"\n", ProblemType);
  fprintf(fptr, "HydroMethod            = %"ISYM"\n", HydroMethod);
  fprintf(fptr, "huge_number            = %e\n", huge_number);
  fprintf(fptr, "tiny_number            = %e\n", tiny_number);
  fprintf(fptr, "Gamma                  = %"GSYM"\n", Gamma);
  fprintf(fptr, "PressureFree           = %"ISYM"\n", PressureFree);
  fprintf(fptr, "RefineBy               = %"ISYM"\n", RefineBy);
  fprintf(fptr, "MaximumRefinementLevel = %"ISYM"\n", MaximumRefinementLevel);
  fprintf(fptr, "MaximumGravityRefinementLevel = %"ISYM"\n",
	  MaximumGravityRefinementLevel);
  fprintf(fptr, "MaximumParticleRefinementLevel = %"ISYM"\n",
	  MaximumParticleRefinementLevel);
  fprintf(fptr, "CellFlaggingMethod     = ");
  WriteListOfInts(fptr, MAX_FLAGGING_METHODS, CellFlaggingMethod);
  fprintf(fptr, "FluxCorrection         = %"ISYM"\n", FluxCorrection);
  fprintf(fptr, "InterpolationMethod    = %"ISYM"\n", InterpolationMethod);
  fprintf(fptr, "ConservativeInterpolation = %"ISYM"\n", ConservativeInterpolation);
  fprintf(fptr, "MinimumEfficiency      = %"GSYM"\n", MinimumEfficiency);
  fprintf(fptr, "MinimumSubgridEdge     = %"ISYM"\n", MinimumSubgridEdge);
  fprintf(fptr, "MaximumSubgridSize     = %"ISYM"\n", MaximumSubgridSize);
  fprintf(fptr, "NumberOfBufferZones    = %"ISYM"\n\n", NumberOfBufferZones);
  fprintf(fptr, "MustRefineRegionMinRefinementLevel = %"ISYM"\n", MustRefineRegionMinRefinementLevel);
  fprintf(fptr, "MetallicityRefinementMinLevel = %"ISYM"\n", MetallicityRefinementMinLevel);
  fprintf(fptr, "MetallicityRefinementMinMetallicity      = %"GSYM"\n", 
	  MetallicityRefinementMinMetallicity);

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
  fprintf(fptr, "MustRefineRegionLeftEdge   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MustRefineRegionLeftEdge);
  fprintf(fptr, "MustRefineRegionRightEdge  = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MustRefineRegionRightEdge);
  fprintf(fptr, "\n");
 
  for (dim = 0; dim < MAX_NUMBER_OF_BARYON_FIELDS; dim++) {
    if (DataLabel[dim])
      fprintf(fptr, "DataLabel[%"ISYM"]              = %s\n", dim, DataLabel[dim]);
    if (DataUnits[dim])
      fprintf(fptr, "DataUnits[%"ISYM"]              = %s\n", dim, DataUnits[dim]);
    if (DataLabel[dim]) {
      if (strstr(DataLabel[dim], "Density") != NULL)
	fprintf(fptr, "#DataCGSConversionFactor[%"ISYM"] = %"GSYM"\n", dim,DensityUnits);
      if (strstr(DataLabel[dim], "velocity") != NULL)
	fprintf(fptr, "#DataCGSConversionFactor[%"ISYM"] = %"GSYM"\n",dim,VelocityUnits);
    }
  }
  fprintf(fptr, "\n");
 
  fprintf(fptr, "UniformGravity             = %"ISYM"\n", UniformGravity);
  fprintf(fptr, "UniformGravityDirection    = %"ISYM"\n", UniformGravityDirection);
  fprintf(fptr, "UniformGravityConstant     = %"GSYM"\n", UniformGravityConstant);
 
  fprintf(fptr, "PointSourceGravity           = %"ISYM"\n",PointSourceGravity);
  fprintf(fptr, "PointSourceGravityPosition   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, PointSourceGravityPosition);
  fprintf(fptr, "PointSourceGravityConstant   = %"GSYM"\n",
	  PointSourceGravityConstant);
  fprintf(fptr, "PointSourceGravityCoreRadius = %"GSYM"\n\n",
	  PointSourceGravityCoreRadius);
 
  fprintf(fptr, "SelfGravity                    = %"ISYM"\n", SelfGravity);
  fprintf(fptr, "GravitationalConstant          = %e\n",
	  GravitationalConstant);
  fprintf(fptr, "S2ParticleSize                 = %"GSYM"\n", S2ParticleSize);
  fprintf(fptr, "GravityResolution              = %"GSYM"\n", GravityResolution);
  fprintf(fptr, "ComputePotential               = %"ISYM"\n", ComputePotential);
  fprintf(fptr, "WritePotential                 = %"ISYM"\n", WritePotential);
  fprintf(fptr, "BaryonSelfGravityApproximation = %"ISYM"\n\n",
	  BaryonSelfGravityApproximation);
 
  fprintf(fptr, "GreensFunctionMaxNumber     = %"ISYM"\n", GreensFunctionMaxNumber);
  fprintf(fptr, "GreensFunctionMaxSize       = %"ISYM"\n", GreensFunctionMaxSize);
 
  fprintf(fptr, "DualEnergyFormalism         = %"ISYM"\n", DualEnergyFormalism);
  fprintf(fptr, "DualEnergyFormalismEta1     = %e\n", DualEnergyFormalismEta1);
  fprintf(fptr, "DualEnergyFormalismEta2     = %e\n", DualEnergyFormalismEta2);
  fprintf(fptr, "ParticleCourantSafetyNumber = %"GSYM"\n\n",
	  ParticleCourantSafetyNumber);
  fprintf(fptr, "RandomForcing               = %"ISYM"\n", RandomForcing);     //AK
  fprintf(fptr, "RandomForcingEdot           = %"GSYM"\n", RandomForcingEdot); //AK
  fprintf(fptr, "RadiativeCooling               = %"ISYM"\n", RadiativeCooling);
  fprintf(fptr, "GadgetEquilibriumCooling       = %"ISYM"\n", 
	  GadgetEquilibriumCooling);
  fprintf(fptr, "MultiSpecies                   = %"ISYM"\n", MultiSpecies);
  fprintf(fptr, "CloudyCooling                  = %"ISYM"\n", CloudyCooling);
  fprintf(fptr, "CloudyCoolingGridRank          = %"ISYM"\n", 
	  CloudyCoolingData.CloudyCoolingGridRank);
  fprintf(fptr, "CloudyCoolingGridRunFile       = %s\n",
	  CloudyCoolingData.CloudyCoolingGridRunFile);
  fprintf(fptr, "IncludeCloudyHeating           = %"ISYM"\n",
	  CloudyCoolingData.IncludeCloudyHeating);
  fprintf(fptr, "CMBTemperatureFloor            = %"ISYM"\n", CMBTemperatureFloor);
  fprintf(fptr, "RadiationFieldType             = %"ISYM"\n", RadiationFieldType);
  fprintf(fptr, "AdjustUVBackground             = %"ISYM"\n", AdjustUVBackground);
  fprintf(fptr, "SetUVBAmplitude                = %"GSYM"\n", SetUVBAmplitude);
  fprintf(fptr, "SetHeIIHeatingScale            = %"GSYM"\n", SetHeIIHeatingScale);
  fprintf(fptr, "RadiationFieldLevelRecompute   = %"ISYM"\n",
	  RadiationFieldLevelRecompute);
  fprintf(fptr, "RadiationSpectrumNormalization = %"GSYM"\n", CoolData.f3);
  fprintf(fptr, "RadiationSpectrumSlope         = %"GSYM"\n", CoolData.alpha0);
 
  fprintf(fptr, "ZEUSLinearArtificialViscosity    = %"GSYM"\n",
	  ZEUSLinearArtificialViscosity);
  fprintf(fptr, "ZEUSQuadraticArtificialViscosity = %"GSYM"\n",
	  ZEUSQuadraticArtificialViscosity);
  fprintf(fptr, "UseMinimumPressureSupport        = %"ISYM"\n",
	  UseMinimumPressureSupport);
  fprintf(fptr, "MinimumPressureSupportParameter  = %"FSYM"\n",
	  MinimumPressureSupportParameter);
  fprintf(fptr, "RefineByJeansLengthSafetyFactor  = %"FSYM"\n",
	  RefineByJeansLengthSafetyFactor);
  fprintf(fptr, "MustRefineParticlesRefineToLevel = %"ISYM"\n",
          MustRefineParticlesRefineToLevel);
  fprintf(fptr, "ParticleTypeInFile               = %"ISYM"\n",
          ParticleTypeInFile);
 
  for (dim = 0; dim < MAX_STATIC_REGIONS; dim++)
    if (StaticRefineRegionLevel[dim] != INT_UNDEFINED) {
      fprintf(fptr, "StaticRefineRegionLevel[%"ISYM"] = %"ISYM"\n", dim,
	      StaticRefineRegionLevel[dim]);
      fprintf(fptr, "StaticRefineRegionLeftEdge[%"ISYM"] = ", dim);
      WriteListOfFloats(fptr, MAX_DIMENSION, StaticRefineRegionLeftEdge[dim]);
      fprintf(fptr, "StaticRefineRegionRightEdge[%"ISYM"] = ", dim);
      WriteListOfFloats(fptr, MAX_DIMENSION, StaticRefineRegionRightEdge[dim]);
    }
 
  fprintf(fptr, "ParallelRootGridIO              = %"ISYM"\n", ParallelRootGridIO);
  fprintf(fptr, "ParallelParticleIO              = %"ISYM"\n", ParallelParticleIO);
  fprintf(fptr, "Unigrid                         = %"ISYM"\n", Unigrid);
  fprintf(fptr, "PartitionNestedGrids            = %"ISYM"\n", PartitionNestedGrids);
  fprintf(fptr, "ExtractFieldsOnly               = %"ISYM"\n", ExtractFieldsOnly);
  fprintf(fptr, "CubeDumpEnabled                 = %"ISYM"\n", CubeDumpEnabled);
 
  if (SRBprefix != NULL)
    fprintf(fptr, "SRBprefix           = %s\n", SRBprefix);

#ifdef OOC_BOUNDARY

  fprintf(fptr, "ExternalBoundaryIO              = %"ISYM"\n", ExternalBoundaryIO);
  fprintf(fptr, "ExternalBoundaryTypeIO          = %"ISYM"\n", ExternalBoundaryTypeIO);
  fprintf(fptr, "ExternalBoundaryValueIO         = %"ISYM"\n", ExternalBoundaryValueIO);

#endif
 
  fprintf(fptr, "MinimumOverDensityForRefinement = %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
       MinimumOverDensityForRefinement[0], MinimumOverDensityForRefinement[1],
       MinimumOverDensityForRefinement[2], MinimumOverDensityForRefinement[3],
       MinimumOverDensityForRefinement[4], MinimumOverDensityForRefinement[5],
       MinimumOverDensityForRefinement[6]);

  fprintf(fptr, "MinimumMassForRefinement = %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM"\n",
       MinimumMassForRefinement[0], MinimumMassForRefinement[1],
       MinimumMassForRefinement[2], MinimumMassForRefinement[3],
       MinimumMassForRefinement[4], MinimumMassForRefinement[5],
       MinimumMassForRefinement[6]);

  fprintf(fptr, "MinimumMassForRefinementLevelExponent = %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
	  MinimumMassForRefinementLevelExponent[0],
	  MinimumMassForRefinementLevelExponent[1],
	  MinimumMassForRefinementLevelExponent[2],
	  MinimumMassForRefinementLevelExponent[3],
	  MinimumMassForRefinementLevelExponent[4],
	  MinimumMassForRefinementLevelExponent[5],
	  MinimumMassForRefinementLevelExponent[6]);

  fprintf(fptr, "MinimumSlopeForRefinement             = %e\n",
	  MinimumSlopeForRefinement);
  fprintf(fptr, "MinimumShearForRefinement             = %e\n",
	  MinimumShearForRefinement);
  fprintf(fptr, "MinimumPressureJumpForRefinement      = %e\n",
	  MinimumPressureJumpForRefinement);
  fprintf(fptr, "MinimumEnergyRatioForRefinement       = %e\n",
	  MinimumEnergyRatioForRefinement);
  fprintf(fptr, "ComovingCoordinates                   = %"ISYM"\n",
	  ComovingCoordinates);
  fprintf(fptr, "StarParticleCreation                  = %"ISYM"\n",
	  StarParticleCreation);
  fprintf(fptr, "StarParticleFeedback                  = %"ISYM"\n",
	  StarParticleFeedback);
  fprintf(fptr, "NumberOfParticleAttributes            = %"ISYM"\n",
	  NumberOfParticleAttributes);
  fprintf(fptr, "StarMakerOverDensityThreshold         = %"GSYM"\n",
	  StarMakerOverDensityThreshold);
  fprintf(fptr, "StarMakerMassEfficiency               = %"GSYM"\n",
	  StarMakerMassEfficiency);
  fprintf(fptr, "StarMakerMinimumMass                  = %"GSYM"\n",
          StarMakerMinimumMass);
  fprintf(fptr, "StarMakerMinimumDynamicalTime         = %"GSYM"\n",
          StarMakerMinimumDynamicalTime);
  fprintf(fptr, "StarMassEjectionFraction              = %"GSYM"\n",
          StarMassEjectionFraction);
  fprintf(fptr, "StarMetalYield                        = %"GSYM"\n",
          StarMetalYield);
  fprintf(fptr, "StarEnergyToThermalFeedback           = %"GSYM"\n",
          StarEnergyToThermalFeedback);
  fprintf(fptr, "StarEnergyToStellarUV                 = %"GSYM"\n",
          StarEnergyToStellarUV);
  fprintf(fptr, "StarEnergyToQuasarUV                  = %"GSYM"\n\n",
          StarEnergyToQuasarUV);
  fprintf(fptr, "MultiMetals                           = %"ISYM"\n",
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
 
  fprintf(fptr, "VersionNumber              = %"FSYM"\n\n", VERSION);
 
  return SUCCESS;
}

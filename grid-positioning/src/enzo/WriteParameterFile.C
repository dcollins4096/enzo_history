/***********************************************************************
/
/  WRITES A PARAMETER FILE (i.e. TopGrid data)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       February 29th, 2008
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
#include <time.h>
#include "ErrorExceptions.h"
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
int  WriteUnits(FILE *fptr);
int  GetUnits(float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MAssUnits, FLOAT Time);
#ifdef TRANSFER
int RadiativeTransferWriteParameters(FILE *fptr);
int WritePhotonSources(FILE *fptr, FLOAT CurrentTime);
#endif /* TRANSFER */
 
int WriteParameterFile(FILE *fptr, TopGridData &MetaData)
{
 
  int dim;
 
  /* Compute Units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits,  MetaData.Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    ENZO_FAIL("");
  }
 
  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0, velu = 1.0, presu = 1.0;
  double massu = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&rhou, &lenu, &tempu, &tu, &velu, &massu, MetaData.Time);
    presu = rhou*lenu*lenu/tu/tu;
  }
  double mh = 1.6726e-24;
  double uheat = VelocityUnits*VelocityUnits*2.0*mh/TimeUnits;

  /* write data to Parameter output file */
 
  /* write MetaData parameters */
 
  fprintf(fptr, "InitialCycleNumber  = %"ISYM"\n", MetaData.CycleNumber);
  fprintf(fptr, "InitialTime         = %"GOUTSYM"\n", MetaData.Time);
  fprintf(fptr, "InitialCPUTime      = %"GSYM"\n\n", MetaData.CPUTime);
 
  fprintf(fptr, "CheckpointRestart   = %"ISYM"\n", CheckpointRestart);
  fprintf(fptr, "StopTime            = %"GOUTSYM"\n", MetaData.StopTime);
  fprintf(fptr, "StopCycle           = %"ISYM"\n", MetaData.StopCycle);
  fprintf(fptr, "StopSteps           = %"ISYM"\n", MetaData.StopSteps);
  fprintf(fptr, "StopCPUTime         = %lg\n", MetaData.StopCPUTime);
  fprintf(fptr, "ResubmitOn          = %"ISYM"\n", MetaData.ResubmitOn);
  fprintf(fptr, "ResubmitCommand     = %s\n\n", MetaData.ResubmitCommand);
 
  fprintf(fptr, "TimeLastRestartDump = %"GOUTSYM"\n", MetaData.TimeLastRestartDump);
  fprintf(fptr, "dtRestartDump       = %"GOUTSYM"\n", MetaData.dtRestartDump);
  fprintf(fptr, "TimeLastDataDump    = %"GOUTSYM"\n", MetaData.TimeLastDataDump);
  fprintf(fptr, "dtDataDump          = %"GOUTSYM"\n", MetaData.dtDataDump);
  fprintf(fptr, "TimeLastHistoryDump = %"GOUTSYM"\n", MetaData.TimeLastHistoryDump);
  fprintf(fptr, "dtHistoryDump       = %"GOUTSYM"\n\n", MetaData.dtHistoryDump);
 
  fprintf(fptr, "TracerParticleOn           = %"ISYM"\n", TracerParticleOn);
  fprintf(fptr, "TimeLastTracerParticleDump = %"GOUTSYM"\n",
          MetaData.TimeLastTracerParticleDump);
  fprintf(fptr, "dtTracerParticleDump       = %"GOUTSYM"\n",
          MetaData.dtTracerParticleDump);
 
  fprintf(fptr, "NewMovieLeftEdge     = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.NewMovieLeftEdge);
  fprintf(fptr, "NewMovieRightEdge    = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.NewMovieRightEdge);
  fprintf(fptr, "MovieSkipTimestep    = %"ISYM"\n", MovieSkipTimestep);
  fprintf(fptr, "Movie3DVolumes       = %"ISYM"\n", Movie3DVolumes);
  fprintf(fptr, "MovieVertexCentered  = %"ISYM"\n", MovieVertexCentered);
  fprintf(fptr, "NewMovieParticleOn   = %"ISYM"\n", NewMovieParticleOn);
  fprintf(fptr, "MovieDataField       = ");
  WriteListOfInts(fptr, MAX_MOVIE_FIELDS, MovieDataField);
  fprintf(fptr, "NewMovieDumpNumber   = %"ISYM"\n", NewMovieDumpNumber);
  fprintf(fptr, "NewMovieName         = %s\n", NewMovieName);
  fprintf(fptr, "MovieTimestepCounter = %"ISYM"\n", MetaData.TimestepCounter);
  fprintf(fptr, "\n");

  fprintf(fptr, "CycleLastRestartDump = %"ISYM"\n", MetaData.CycleLastRestartDump);
  fprintf(fptr, "CycleSkipRestartDump = %"ISYM"\n", MetaData.CycleSkipRestartDump);
  fprintf(fptr, "CycleLastDataDump    = %"ISYM"\n", MetaData.CycleLastDataDump);
  fprintf(fptr, "CycleSkipDataDump    = %"ISYM"\n", MetaData.CycleSkipDataDump);
  fprintf(fptr, "CycleLastHistoryDump = %"ISYM"\n", MetaData.CycleLastHistoryDump);
  fprintf(fptr, "CycleSkipHistoryDump = %"ISYM"\n\n",
	  MetaData.CycleSkipHistoryDump);


  fprintf(fptr, "PythonSubcycleSkip      = %"ISYM"\n", PythonSubcycleSkip);
  fprintf(fptr, "CycleSkipGlobalDataDump = %"ISYM"\n\n", //AK
          MetaData.CycleSkipGlobalDataDump);

  fprintf(fptr, "SubcycleNumber          = %"ISYM"\n", MetaData.SubcycleNumber);
  fprintf(fptr, "SubcycleSkipDataDump    = %"ISYM"\n", MetaData.SubcycleSkipDataDump);
  fprintf(fptr, "SubcycleLastDataDump    = %"ISYM"\n", MetaData.SubcycleLastDataDump);
   fprintf(fptr, "OutputFirstTimeAtLevel = %"ISYM"\n",
	  MetaData.OutputFirstTimeAtLevel);
  fprintf(fptr, "StopFirstTimeAtLevel    = %"ISYM"\n\n",
	  MetaData.StopFirstTimeAtLevel);
 
  fprintf(fptr, "RestartDumpNumber   = %"ISYM"\n", MetaData.RestartDumpNumber);
  fprintf(fptr, "DataDumpNumber      = %"ISYM"\n", MetaData.DataDumpNumber);
  fprintf(fptr, "HistoryDumpNumber   = %"ISYM"\n", MetaData.HistoryDumpNumber);
  fprintf(fptr, "TracerParticleDumpNumber = %"ISYM"\n",
          MetaData.TracerParticleDumpNumber);
 
  fprintf(fptr, "RestartDumpName     = %s\n", MetaData.RestartDumpName);
  fprintf(fptr, "DataDumpName        = %s\n", MetaData.DataDumpName);
  fprintf(fptr, "HistoryDumpName     = %s\n", MetaData.HistoryDumpName);
  fprintf(fptr, "TracerParticleDumpName = %s\n",
          MetaData.TracerParticleDumpName);
  fprintf(fptr, "RedshiftDumpName    = %s\n\n", MetaData.RedshiftDumpName);
 
  if (MetaData.RestartDumpDir != NULL)
    fprintf(fptr, "RestartDumpDir        = %s\n", MetaData.RestartDumpDir);
  if (MetaData.DataDumpDir != NULL)
    fprintf(fptr, "DataDumpDir           = %s\n", MetaData.DataDumpDir);
  if (MetaData.HistoryDumpDir != NULL)
    fprintf(fptr, "HistoryDumpDir        = %s\n", MetaData.HistoryDumpDir);
  if (MetaData.TracerParticleDumpDir != NULL)
    fprintf(fptr, "TracerParticleDumpDir = %s\n", MetaData.TracerParticleDumpDir);
  if (MetaData.RedshiftDumpDir != NULL)
    fprintf(fptr, "RedshiftDumpDir       = %s\n\n", MetaData.RedshiftDumpDir);
 
  if (MetaData.LocalDir != NULL)
    fprintf(fptr, "LocalDir            = %s\n", MetaData.LocalDir);
  if (MetaData.GlobalDir != NULL)
    fprintf(fptr, "GlobalDir           = %s\n", MetaData.GlobalDir);
 
  for (dim = 0; dim < MAX_CUBE_DUMPS; dim++)
    if (CubeDumps[dim] != NULL)
      fprintf(fptr, "CubeDump[%"ISYM"]            = %s\n", dim, CubeDumps[dim]);

  fprintf(fptr, "LoadBalancing          = %"ISYM"\n", LoadBalancing);
  fprintf(fptr, "LoadBalancingCycleSkip = %"ISYM"\n", LoadBalancingCycleSkip);
 
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
 
  fprintf(fptr, "CourantSafetyNumber    = %"FSYM"\n",
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
  fprintf(fptr, "RefineRegionAutoAdjust = %"ISYM"\n", RefineRegionAutoAdjust);
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
  fprintf(fptr, "#TimeUnits                 = %"GSYM"\n", TimeUnits);
  fprintf(fptr, "#TemperatureUnits          = %"GSYM"\n", TemperatureUnits);
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
  fprintf(fptr, "PotentialIterations            = %"ISYM"\n", PotentialIterations);
  fprintf(fptr, "WritePotential                 = %"ISYM"\n", WritePotential);
  fprintf(fptr, "BaryonSelfGravityApproximation = %"ISYM"\n\n",
	  BaryonSelfGravityApproximation);

  fprintf(fptr, "InlineHaloFinder               = %"ISYM"\n", InlineHaloFinder);
  fprintf(fptr, "HaloFinderSubfind              = %"ISYM"\n", HaloFinderSubfind);
  fprintf(fptr, "HaloFinderCycleSkip            = %"ISYM"\n", 
	  HaloFinderCycleSkip);
  fprintf(fptr, "HaloFinderOutputParticleList   = %"ISYM"\n", 
	  HaloFinderOutputParticleList);
  fprintf(fptr, "HaloFinderMinimumSize          = %"ISYM"\n", 
	  HaloFinderMinimumSize);
  fprintf(fptr, "HaloFinderLinkingLength        = %"FSYM"\n",
	  HaloFinderLinkingLength);
  fprintf(fptr, "HaloFinderTimestep             = %"FSYM"\n",
	  HaloFinderTimestep);
  fprintf(fptr, "HaloFinderLastTime             = %"PSYM"\n\n", 
	  HaloFinderLastTime);
 
  fprintf(fptr, "GreensFunctionMaxNumber     = %"ISYM"\n", GreensFunctionMaxNumber);
  fprintf(fptr, "GreensFunctionMaxSize       = %"ISYM"\n", GreensFunctionMaxSize);
 
  fprintf(fptr, "DualEnergyFormalism         = %"ISYM"\n", DualEnergyFormalism);
  fprintf(fptr, "DualEnergyFormalismEta1     = %e\n", DualEnergyFormalismEta1);
  fprintf(fptr, "DualEnergyFormalismEta2     = %e\n", DualEnergyFormalismEta2);
  fprintf(fptr, "ParticleCourantSafetyNumber = %"FSYM"\n\n", ParticleCourantSafetyNumber);
  fprintf(fptr, "RootGridCourantSafetyNumber = %"FSYM"\n\n", RootGridCourantSafetyNumber);
  fprintf(fptr, "RandomForcing                  = %"ISYM"\n", RandomForcing);
  fprintf(fptr, "RandomForcingEdot              = %"GSYM"\n", RandomForcingEdot);
  fprintf(fptr, "RadiativeCooling               = %"ISYM"\n", RadiativeCooling);
  fprintf(fptr, "GadgetEquilibriumCooling       = %"ISYM"\n", GadgetEquilibriumCooling);
  fprintf(fptr, "MultiSpecies                   = %"ISYM"\n", MultiSpecies);
  fprintf(fptr, "PrimordialChemistrySolver      = %"ISYM"\n", PrimordialChemistrySolver);
  fprintf(fptr, "CIECooling                     = %"ISYM"\n", CIECooling);
  fprintf(fptr, "H2OpticalDepthApproximation    = %"ISYM"\n", H2OpticalDepthApproximation);
  fprintf(fptr, "ThreeBodyRate                  = %"ISYM"\n", ThreeBodyRate);
  fprintf(fptr, "CloudyCoolingGridFile          = %s\n", CloudyCoolingData.CloudyCoolingGridFile);
  fprintf(fptr, "IncludeCloudyHeating           = %"ISYM"\n", CloudyCoolingData.IncludeCloudyHeating);
  fprintf(fptr, "IncludeCloudyMMW               = %"ISYM"\n", CloudyCoolingData.IncludeCloudyMMW);
  fprintf(fptr, "CMBTemperatureFloor            = %"ISYM"\n", CloudyCoolingData.CMBTemperatureFloor);
  fprintf(fptr, "ConstantTemperatureFloor       = %"FSYM"\n", CloudyCoolingData.ConstantTemperatureFloor);
  fprintf(fptr, "CloudyMetallicityNormalization = %"FSYM"\n", CloudyCoolingData.CloudyMetallicityNormalization);
  fprintf(fptr, "CloudyElectronFractionFactor   = %"FSYM"\n", CloudyCoolingData.CloudyElectronFractionFactor);
  fprintf(fptr, "MetalCooling                   = %"ISYM"\n", MetalCooling);
  fprintf(fptr, "MetalCoolingTable              = %s\n", MetalCoolingTable);
  fprintf(fptr, "RadiativeTransfer              = %"ISYM"\n", RadiativeTransfer);
  fprintf(fptr, "RadiationXRaySecondaryIon      = %"ISYM"\n", 
	  RadiationXRaySecondaryIon);
  fprintf(fptr, "RadiationFieldType             = %"ISYM"\n", RadiationFieldType);
  fprintf(fptr, "AdjustUVBackground             = %"ISYM"\n", AdjustUVBackground);
  fprintf(fptr, "SetUVBAmplitude                = %"GSYM"\n", SetUVBAmplitude);
  fprintf(fptr, "SetHeIIHeatingScale            = %"GSYM"\n", SetHeIIHeatingScale);
  fprintf(fptr, "RadiationFieldLevelRecompute   = %"ISYM"\n", RadiationFieldLevelRecompute);
  fprintf(fptr, "RadiationSpectrumNormalization = %"GSYM"\n", CoolData.f3);
  fprintf(fptr, "RadiationSpectrumSlope         = %"GSYM"\n", CoolData.alpha0);

  if (CoolData.ParameterFilename != NULL)
    fprintf(fptr, "CoolDataParameterFile = %s\n\n", CoolData.ParameterFilename);

  fprintf(fptr, "OutputCoolingTime              = %"ISYM"\n", OutputCoolingTime);
  fprintf(fptr, "OutputTemperature              = %"ISYM"\n", OutputTemperature);
  fprintf(fptr, "OutputSmoothedDarkMatter       = %"ISYM"\n", 
	  OutputSmoothedDarkMatter);
  fprintf(fptr, "SmoothedDarkMatterNeighbors    = %"ISYM"\n", 
	  SmoothedDarkMatterNeighbors);
 
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
  fprintf(fptr, "RefineByResistiveLengthSafetyFactor  = %"FSYM"\n", 
	  RefineByResistiveLengthSafetyFactor);
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
  fprintf(fptr, "UnigridTranspose                = %"ISYM"\n", UnigridTranspose);
  fprintf(fptr, "PartitionNestedGrids            = %"ISYM"\n", PartitionNestedGrids);
  fprintf(fptr, "ExtractFieldsOnly               = %"ISYM"\n", ExtractFieldsOnly);
  fprintf(fptr, "CubeDumpEnabled                 = %"ISYM"\n", CubeDumpEnabled);
 
  fprintf(fptr, "Debug1                          = %"ISYM"\n", debug1);
  fprintf(fptr, "Debug2                          = %"ISYM"\n", debug2);

  fprintf(fptr, "MemoryLimit                     = %"ISYM"\n", MemoryLimit);

#ifdef STAGE_INPUT
  fprintf(fptr, "StageInput                      = %"ISYM"\n", StageInput);
  fprintf(fptr, "LocalPath                       = %s\n", LocalPath);
  fprintf(fptr, "GlobalPath                      = %s\n", GlobalPath);
#endif

#ifdef OOC_BOUNDARY

  fprintf(fptr, "ExternalBoundaryIO              = %"ISYM"\n", ExternalBoundaryIO);
  fprintf(fptr, "ExternalBoundaryTypeIO          = %"ISYM"\n", ExternalBoundaryTypeIO);
  fprintf(fptr, "ExternalBoundaryValueIO         = %"ISYM"\n", ExternalBoundaryValueIO);
  fprintf(fptr, "SimpleConstantBoundary          = %"ISYM"\n", SimpleConstantBoundary);

#endif
 

  fprintf(fptr, "SlopeFlaggingFields ="
	  " %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",
	  SlopeFlaggingFields[0], 
	  SlopeFlaggingFields[1],
	  SlopeFlaggingFields[2], 
	  SlopeFlaggingFields[3],
	  SlopeFlaggingFields[4]);

  fprintf(fptr, "MinimumSlopeForRefinement ="
	  " %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
	  MinimumSlopeForRefinement[0],
	  MinimumSlopeForRefinement[1],
	  MinimumSlopeForRefinement[2],
	  MinimumSlopeForRefinement[3],
	  MinimumSlopeForRefinement[4],
	  MinimumSlopeForRefinement[5],
	  MinimumSlopeForRefinement[6]);


  fprintf(fptr, "MinimumOverDensityForRefinement ="
	  " %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
	  MinimumOverDensityForRefinement[0],
	  MinimumOverDensityForRefinement[1],
	  MinimumOverDensityForRefinement[2],
	  MinimumOverDensityForRefinement[3],
	  MinimumOverDensityForRefinement[4],
	  MinimumOverDensityForRefinement[5],
	  MinimumOverDensityForRefinement[6]);

  fprintf(fptr, "MinimumMassForRefinement ="
	  " %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM"\n",
	  MinimumMassForRefinement[0],
	  MinimumMassForRefinement[1],
	  MinimumMassForRefinement[2],
	  MinimumMassForRefinement[3],
	  MinimumMassForRefinement[4],
	  MinimumMassForRefinement[5],
	  MinimumMassForRefinement[6]);

  fprintf(fptr, "MinimumMassForRefinementLevelExponent ="
	  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
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
  fprintf(fptr, "AddParticleAttributes                 = %"ISYM"\n", 
	  AddParticleAttributes);

    /* Sink particles (for present day star formation) & winds */
  fprintf(fptr, "SinkMergeDistance                     = %"FSYM"\n", 
	  SinkMergeDistance);
  fprintf(fptr, "SinkMergeMass                         = %"FSYM"\n", 
	  SinkMergeMass);
  fprintf(fptr, "StellarWindFeedback                   = %"ISYM"\n", 
	  StellarWindFeedback);
  fprintf(fptr, "StellarWindTurnOnMass                 = %"FSYM"\n", 
	  StellarWindTurnOnMass);


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

  fprintf(fptr, "RefineByJeansLengthUnits              = %"ISYM"\n",RefineByJeansLengthUnits);
  fprintf(fptr, "IsothermalSoundSpeed                  = %"GSYM"\n",IsothermalSoundSpeed);
          
  fprintf(fptr, "StarClusterUseMetalField              = %"ISYM"\n",
	  StarClusterUseMetalField);
  fprintf(fptr, "StarClusterMinDynamicalTime           = %"GSYM"\n",
          StarClusterMinDynamicalTime);
  fprintf(fptr, "StarClusterIonizingLuminosity         = %lg\n",
          StarClusterIonizingLuminosity);
  fprintf(fptr, "StarClusterSNEnergy                   = %lg\n",
          StarClusterSNEnergy);
  fprintf(fptr, "StarClusterSNRadius                   = %"GSYM"\n",
          StarClusterSNRadius);
  fprintf(fptr, "StarClusterFormEfficiency             = %"GSYM"\n",
          StarClusterFormEfficiency);
  fprintf(fptr, "StarClusterMinimumMass                = %"GSYM"\n",
          StarClusterMinimumMass);
  fprintf(fptr, "StarClusterCombineRadius              = %"GSYM"\n",
          StarClusterCombineRadius);
  fprintf(fptr, "StarClusterRegionLeftEdge             = %"FSYM" %"FSYM" %"FSYM"\n",
          StarClusterRegionLeftEdge[0], StarClusterRegionLeftEdge[1],
	  StarClusterRegionLeftEdge[2]);
  fprintf(fptr, "StarClusterRegionRightEdge            = %"FSYM" %"FSYM" %"FSYM"\n",
          StarClusterRegionRightEdge[0], StarClusterRegionRightEdge[1],
	  StarClusterRegionRightEdge[2]);
  fprintf(fptr, "PopIIIStarMass                        = %"GSYM"\n",
          PopIIIStarMass);
  fprintf(fptr, "PopIIIBlackHoles                      = %"ISYM"\n",
          PopIIIBlackHoles);
  fprintf(fptr, "PopIIIBHLuminosityEfficiency          = %"FSYM"\n",
          PopIIIBHLuminosityEfficiency);
  fprintf(fptr, "PopIIIOverDensityThreshold            = %"GSYM"\n",
          PopIIIOverDensityThreshold);
  fprintf(fptr, "PopIIIH2CriticalFraction              = %"GSYM"\n",
          PopIIIH2CriticalFraction);
  fprintf(fptr, "PopIIIMetalCriticalFraction           = %"GSYM"\n",
          PopIIIMetalCriticalFraction);
  fprintf(fptr, "PopIIISupernovaRadius                 = %"GSYM"\n",
          PopIIISupernovaRadius);
  fprintf(fptr, "PopIIISupernovaUseColour              = %"ISYM"\n\n",
          PopIIISupernovaUseColour);
  fprintf(fptr, "MBHMinDynamicalTime            = %"GSYM"\n",
          MBHMinDynamicalTime);
  fprintf(fptr, "MBHMinimumMass                 = %"GSYM"\n",
          MBHMinimumMass);
  fprintf(fptr, "MBHFeedbackThermal             = %"ISYM"\n",
	  MBHFeedbackThermal);
  fprintf(fptr, "MBHFeedbackRadius              = %"GSYM"\n",
          MBHFeedbackRadius);
  fprintf(fptr, "MBHFeedbackRadiativeEfficiency = %"GSYM"\n",
          MBHFeedbackRadiativeEfficiency);
  fprintf(fptr, "MBHFeedbackThermalCoupling     = %"GSYM"\n",
          MBHFeedbackThermalCoupling);
  fprintf(fptr, "MBHCombineRadius               = %"GSYM"\n",
          MBHCombineRadius);
  fprintf(fptr, "MBHIonizingLuminosity          = %lg\n\n",
          MBHIonizingLuminosity);

  fprintf(fptr, "PopIIIColorDensityThreshold            = %"GSYM"\n",
          PopIIIColorDensityThreshold);
  fprintf(fptr, "PopIIIColorMass                        = %"GSYM"\n",
          PopIIIColorMass);

  /* Most Stanford additions: */

  fprintf(fptr, "Theta_Limiter = %f\n", Theta_Limiter);
  fprintf(fptr, "RiemannSolver = %d\n", RiemannSolver);
  fprintf(fptr, "ReconstructionMethod = %d\n", ReconstructionMethod);
  fprintf(fptr, "RKOrder = %d\n", RKOrder);
  fprintf(fptr, "UsePhysicalUnit = %d\n", UsePhysicalUnit);
  fprintf(fptr, "UseFloor = %d\n", UseFloor);
  fprintf(fptr, "UseViscosity = %d\n", UseViscosity);
  fprintf(fptr, "UseAmbipolarDiffusion = %d\n", UseAmbipolarDiffusion);
  fprintf(fptr, "UseResistivity = %d\n", UseResistivity);
  fprintf(fptr, "SmallRho = %g\n", SmallRho*rhou);
  fprintf(fptr, "SmallP = %g\n", SmallP*presu);
  fprintf(fptr, "SmallT = %g\n", SmallT*tempu);
  fprintf(fptr, "MaximumAlvenSpeed = %g\n", MaximumAlvenSpeed*velu);
  fprintf(fptr, "Coordinate = %d\n", Coordinate);
  fprintf(fptr, "EOSType = %d\n", EOSType);
  fprintf(fptr, "EOSSoundSpeed = %g\n", EOSSoundSpeed);
  fprintf(fptr, "EOSCriticalDensity = %g\n", EOSCriticalDensity);
  fprintf(fptr, "EOSGamma = %g\n", EOSGamma); 
  fprintf(fptr, "Mu = %g\n", Mu);
  fprintf(fptr, "CoolingCutOffDensity1 = %g\n", CoolingCutOffDensity1);
  fprintf(fptr, "CoolingCutOffDensity2 = %g\n", CoolingCutOffDensity2);
  fprintf(fptr, "CoolingCutOffTemperature = %g\n", CoolingCutOffTemperature);
  fprintf(fptr, "CoolingPowerCutOffDensity1 = %g\n", CoolingPowerCutOffDensity1);
  fprintf(fptr, "CoolingPowerCutOffDensity2 = %g\n", CoolingPowerCutOffDensity2);
  fprintf(fptr, "UseConstantAcceleration = %d\n", UseConstantAcceleration);
  fprintf(fptr, "ConstantAcceleration = %g %g %g\n", ConstantAcceleration[0],
	  ConstantAcceleration[1], ConstantAcceleration[2]);


  fprintf(fptr, "AngularVelocity = %g\n", AngularVelocity);
  fprintf(fptr, "VelocityGradient = %g\n", VelocityGradient);
  fprintf(fptr, "UseDrivingField = %d\n", UseDrivingField);
  fprintf(fptr, "DrivingEfficiency = %f\n", DrivingEfficiency);
#ifdef ECUDA
  fprintf(fptr, "UseCUDA = %f\n", UseCUDA);
#endif

  /* Poisson Solver */

  fprintf(fptr, "DivergenceCleaningBoundaryBuffer = %"ISYM"\n",
	  DivergenceCleaningBoundaryBuffer);
  fprintf(fptr, "UseDivergenceCleaning            = %d\n", UseDivergenceCleaning);
  fprintf(fptr, "DivergenceCleaningThreshold      = %g\n", 
	  DivergenceCleaningThreshold);
  fprintf(fptr, "PoissonApproximationThreshold    = %g\n", 
	  PoissonApproximationThreshold);

  /* Shearing Box Boundary parameters */
  fprintf(fptr, "AngularVelocity              = %"FSYM"\n",AngularVelocity);
  fprintf(fptr, "VelocityGradient             = %"FSYM"\n",VelocityGradient);
  fprintf(fptr, "ShearingVelocityDirection    = %"ISYM"\n\n",ShearingVelocityDirection);
  fprintf(fptr, "ShearingBoxProblemType    = %"ISYM"\n\n", ShearingBoxProblemType);

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
 
  if (ComovingCoordinates) {
    if (CosmologyWriteParameters(fptr, MetaData.StopTime, MetaData.Time) ==
	FAIL) {
      fprintf(stderr, "Error in CosmologyWriteParameters.\n");
      ENZO_FAIL("");
    }
  }
  else {
    if (WriteUnits(fptr) == FAIL) {
      fprintf(stderr, "Error in WriteUnits.\n");
      ENZO_FAIL("");
    }
  }

  /* If radiative transfer, write parameters.  If photon test, write
     source data. */

#ifdef TRANSFER
  if (RadiativeTransferWriteParameters(fptr) == FAIL) {
    fprintf(stderr, "Error in RadiativeTransferWriteParameters.\n");
    ENZO_FAIL("");
  }

  if (ProblemType == 50)
    if (WritePhotonSources(fptr, MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in WritePhotonSources.\n");
      ENZO_FAIL("");
    }
#endif

  /* Output current time */
  time_t ID;
  ID = time(NULL);
  fprintf(fptr, "CurrentTimeIdentifier = %"ISYM"\n", int(ID));

  /* write version info */
 
  fprintf(fptr, "VersionNumber              = %"FSYM"\n\n", VERSION);
 
  return SUCCESS;
}

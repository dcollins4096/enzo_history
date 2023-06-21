/***********************************************************************
/
/  SETS THE DEFAULT GLOBAL VALUES
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

// This routine intializes a new simulation based on the parameter file.
//

#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "TopGridData.h"
#include "StarParticleData.h"

/* character strings */

char DefaultDimUnits[] = "cm";
char *DefaultDimLabel[] = {"x", "y", "z"};
char DefaultRestartName[] = "restart";
char DefaultDataName[] = "data";
char DefaultHistoryName[] = "history";
char DefaultRedshiftName[] = "RedshiftOutput";
char DefaultMovieName[] = "MovieOutput";
char DefaultTracerParticleName[] = "TracerOutput";

int SetDefaultGlobalValues(TopGridData &MetaData)
{

  /* declarations */

  const float Pi = 3.14159;
  int dim, i;

  /* set the default MetaData values. */

  MetaData.CycleNumber     = 0;
  MetaData.Time            = 0.0;
  MetaData.CPUTime         = 0.0;
  
  MetaData.StopTime        = FLOAT_UNDEFINED;  // This must be set be the user
  MetaData.StopCycle       = 10000;            // 10000 timesteps
  MetaData.StopCPUTime     = 100.0*3600.0;     // 100 hours
  
  MetaData.TimeLastRestartDump = 0.0;
  MetaData.dtRestartDump       = 5.0*3600.0;   // every 5 hours
  MetaData.TimeLastDataDump    = FLOAT_UNDEFINED;
  MetaData.dtDataDump          = 0.0;
  MetaData.TimeLastHistoryDump = FLOAT_UNDEFINED;
  MetaData.dtHistoryDump       = 0.0;
  MetaData.TimeLastMovieDump   = FLOAT_UNDEFINED;
  MetaData.dtMovieDump         = 0.0;
  MetaData.TimeLastTracerParticleDump = FLOAT_UNDEFINED;
  MetaData.dtTracerParticleDump       = 0.0;

  MetaData.CycleLastRestartDump = 0;
  MetaData.CycleSkipRestartDump = 0;
  MetaData.CycleLastDataDump    = INT_UNDEFINED;
  MetaData.CycleSkipDataDump    = 0;
  MetaData.CycleLastHistoryDump = INT_UNDEFINED;
  MetaData.CycleSkipHistoryDump = 0;

  MetaData.OutputFirstTimeAtLevel = 0; // zero is off
  MetaData.StopFirstTimeAtLevel   = 0; // zero is off

  MetaData.RestartDumpNumber   = 0;            // starting restart id number
  MetaData.RestartDumpName     = DefaultRestartName;
  MetaData.DataDumpNumber      = 0;
  MetaData.DataDumpName        = DefaultDataName;
  MetaData.HistoryDumpNumber   = 0;
  MetaData.HistoryDumpName     = DefaultHistoryName;
  MetaData.MovieDumpNumber     = 0;
  MetaData.MovieDumpName       = DefaultMovieName;
  MetaData.TracerParticleDumpNumber = 0;
  MetaData.TracerParticleDumpName   = DefaultTracerParticleName;
  MetaData.RedshiftDumpName    = DefaultRedshiftName;

  for (i = 0; i < MAX_TIME_ACTIONS; i++) {
    TimeActionType[i]      = 0;
    TimeActionRedshift[i]  = -1;
    TimeActionTime[i]      = 0;
    TimeActionParameter[i] = FLOAT_UNDEFINED;
  }

  MetaData.StaticHierarchy     = TRUE;

  MetaData.TopGridRank = INT_UNDEFINED;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    MetaData.TopGridDims[dim]                = INT_UNDEFINED;  
    MetaData.LeftFaceBoundaryCondition[dim]  = reflecting;
    MetaData.RightFaceBoundaryCondition[dim] = reflecting;
  }
  MetaData.BoundaryConditionName = NULL;

  MetaData.GravityBoundary        = TopGridPeriodic;

  MetaData.ParticleBoundaryType   = periodic;  // only one implemented!
  MetaData.NumberOfParticles      = 0;         // no particles

  MetaData.CourantSafetyNumber    = 0.6;
  MetaData.PPMFlatteningParameter = 0;    // off
  MetaData.PPMDiffusionParameter  = 1;    // on
  MetaData.PPMSteepeningParameter = 0;    // off

  /* set the default global data. */
                                                 // Debug flag set in main
  ProblemType               = 0;                 // None
  HydroMethod               = PPM_DirectEuler;   //
  huge_number               = 1.0e+20;
  tiny_number               = 1.0e-20;
  Gamma                     = 5.0/3.0;           // 5/3
  PressureFree              = FALSE;             // use pressure (duh)
  RefineBy                  = 4;                 // Refinement factor
  MaximumRefinementLevel    = 2;                 // three levels (w/ topgrid)
  MaximumGravityRefinementLevel = INT_UNDEFINED;
  MaximumParticleRefinementLevel = -1;            // unused if negative
  FluxCorrection            = TRUE;
  InterpolationMethod       = SecondOrderA;      // ?
  ConservativeInterpolation = TRUE;              // true for ppm
  MinimumEfficiency         = 0.2;               // between 0-1, usually ~0.1
  NumberOfBufferZones       = 1;

  for (i = 0; i < MAX_FLAGGING_METHODS; i++) {
    CellFlaggingMethod[i]       = INT_UNDEFINED;
    MinimumMassForRefinement[i] = FLOAT_UNDEFINED;   // usually set by:
    MinimumOverDensityForRefinement[i]       = 1.5; 
    MinimumMassForRefinementLevelExponent[i] = 0;
  }

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    DomainLeftEdge[dim]             = 0.0;
    DomainRightEdge[dim]            = 1.0;
    GridVelocity[dim]               = 0.0;
    DimLabels[dim]                  = DefaultDimLabel[dim];
    DimUnits[dim]                   = DefaultDimUnits;
    RefineRegionLeftEdge[dim]       = FLOAT_UNDEFINED;
    RefineRegionRightEdge[dim]      = FLOAT_UNDEFINED;
    MetaData.MovieRegionLeftEdge[dim]  = FLOAT_UNDEFINED;
    MetaData.MovieRegionRightEdge[dim] = FLOAT_UNDEFINED;
    PointSourceGravityPosition[dim] = 0.0;
  }

  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    StaticRefineRegionLevel[i] = INT_UNDEFINED;

  ParallelRootGridIO          = FALSE;

  UniformGravity              = FALSE;             // off
  UniformGravityDirection     = 0;                 // x-direction
  UniformGravityConstant      = 1.0;

  PointSourceGravity           = FALSE;             // off
  PointSourceGravityConstant   = 1.0;
  PointSourceGravityCoreRadius = 0.0;

  SelfGravity                 = FALSE;             // off
  GravitationalConstant       = 4*Pi;              // G = 1
  S2ParticleSize              = 3.0;               // ~3 is reasonable
  GravityResolution           = 1.0;               // equivalent to grid
  ComputePotential            = FALSE;
  BaryonSelfGravityApproximation = TRUE;           // less accurate but faster

  GreensFunctionMaxNumber     = 1;                 // only one at a time
  GreensFunctionMaxSize       = 1;                 // not used yet

  DualEnergyFormalism         = FALSE;             // off
  DualEnergyFormalismEta1     = 0.001;             // typical 0.001
  DualEnergyFormalismEta2     = 0.1;               // 0.08-0.1
  ParticleCourantSafetyNumber = 0.5;
  RadiativeCooling            = FALSE;             // off
  MultiSpecies                = FALSE;             // off
  RadiationFieldType          = 0;
  RadiationFieldLevelRecompute = 0;
  CoolData.alpha0             = 1.5;               // radiation spectral slope
  CoolData.f3                 = 1.0e-21;           // radiation normalization

  ZEUSLinearArtificialViscosity    = 0.0;
  ZEUSQuadraticArtificialViscosity = 2.0;
  UseMinimumPressureSupport        = FALSE;
  MinimumPressureSupportParameter  = 100.0;
  
  MinimumSlopeForRefinement        = 0.3;          // 30% change in value
  MinimumPressureJumpForRefinement = 0.33;         // As in PPM method paper
  MinimumEnergyRatioForRefinement  = 0.1;          // conservative!
  RefineByJeansLengthSafetyFactor  = 4.0;
  MustRefineParticlesRefineToLevel = 0;
  ComovingCoordinates              = FALSE;        // No comoving coordinates
  StarParticleCreation             = FALSE;
  StarParticleFeedback             = FALSE;
  StarMakerOverDensityThreshold    = 100;          // times mean total density
  StarMakerMassEfficiency          = 1;
  StarMakerMinimumMass             = 1.0e9;        // in solar masses
  StarMakerMinimumDynamicalTime    = 1.0e6;        // in years
  StarMassEjectionFraction         = 0.25;
  StarMetalYield                   = 0.02;
  StarEnergyToThermalFeedback      = 1.0e-5;
  StarEnergyToStellarUV            = 3.0e-6;
  StarEnergyToQuasarUV             = 5.0e-6;
  NumberOfParticleAttributes       = INT_UNDEFINED;
  TracerParticleOn                 = FALSE;

  return SUCCESS;
}

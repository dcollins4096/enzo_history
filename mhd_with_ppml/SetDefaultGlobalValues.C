/***********************************************************************
/
/  SETS THE DEFAULT GLOBAL VALUES
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       October, 2004
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
char DefaultNewMovieName[] = "MoviePack";
char DefaultTracerParticleName[] = "TracerOutput";
 
char DefaultRestartDir[] = "RS";
char DefaultDataDir[] = "DD";
char DefaultHistoryDir[] = "HD";
char DefaultRedshiftDir[] = "RD";
char DefaultMovieDir[] = "MD";
char DefaultTracerParticleDir[] = "TD";
 
 
 
 
int SetDefaultGlobalValues(TopGridData &MetaData)
{
 
  /* declarations */
 
  const float Pi = 3.14159;
  int dim, i;

#ifdef MHDF
  MHD_DivB = 4;
  MHD_DivBparam = 0.0;
  MHD_WriteElectric = 0;
  MHD_CenteringMethod = 6;
#endif //MHDF 
#ifdef PPML
  /* PPML stuff. */

  // ATH || PPML
  for( i=0; i< NUMBER_OF_INTEGRATION_STEPS; i++){
    MHD_Recon[i]              = 3; // PPML by default.  The other option is 0, for Piecewise Constant.
    MHD_Riemann[i]            = -1; // 
    MHD_DiffusionMethod[i]    = 0;  
  }
  MHD_DiffusionParameter    = 0.0;
  MHD_DirectionalSplitting  = 0;  //Turn directional Splitting on or off. (0 for unsplit.)
  //All default to on.
  for( i=0; i<3; i++)
    PPML_SlopeLimiter[i] = 1;
  PPML_EigenSystem       = 0; //0 for Powell, 1 for Ryu and Jones

  // Physics Parameters
  MHD_Used                  = 0;
  EquationOfState           = 0;                 //0 for adiabatic, 1 for isothermal.
  IsothermalSoundSpeed      = 1.0;               //For isothermal

  // Problem parameters
  for( i=0 ;i<4; i++)
    DiscontNormal[i]        = 0.0;               //For oblique shocks

  PPML_InitInterfaceMethod  = 1;

  // Counters for debugging.
  dccCounter01 = 0;  // Flux output, Grid_PPML_Update
  dccCounter02 = 0;  // Forced Timestep Read, Grid_ComputeTimestep
  dccCounter03 = 0;  // unused
  dccCounter04 = 0;  // unused
  dccCounter05 = 0;  // unused
  dccCounter06 = 0;  // unused
  dccCounter07 = 0;  // unused
  dccCounter08 = 0;  // unused
  dccCounter09 = 0;  // unused
  dccCounter10 = 0;  // unused
  dccCounter11 = 0;  // unused

  //Extraneous Debug Controls.
  for( int mid=0; mid<N_MidWayDumps; mid++) MidWayDumpList[mid] = 0;
  WriteBoundary = 0;
  WriteAcceleration = 0;
  FixedTimestep = -1.0;
  ProcessorTopology[0]      = INT_UNDEFINED;
  ProcessorTopology[1]      = INT_UNDEFINED;
  ProcessorTopology[2]      = INT_UNDEFINED;



#ifdef DCC_EXTERNAL_CYCLE_COUNT
  for(i=0; i<MAX_DEPTH_OF_HIERARCHY; i++)
    LevelCycleCount[i] = 0;
#endif //!DCC_EXTERNAL_CYCLE_COUNT

#endif //PPML

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
  MetaData.CycleSkipGlobalDataDump = 0; //AK
 
  MetaData.OutputFirstTimeAtLevel = 0; // zero is off
  MetaData.StopFirstTimeAtLevel   = 0; // zero is off
 
  MetaData.RestartDumpNumber   = 0;            // starting restart id number
  MetaData.RestartDumpName     = DefaultRestartName;
  MetaData.RestartDumpDir      = DefaultRestartDir;
  MetaData.DataDumpNumber      = 0;
  MetaData.DataDumpName        = DefaultDataName;
  MetaData.DataDumpDir         = DefaultDataDir;
  MetaData.HistoryDumpNumber   = 0;
  MetaData.HistoryDumpName     = DefaultHistoryName;
  MetaData.HistoryDumpDir      = DefaultHistoryDir;
  MetaData.MovieDumpNumber     = 0;
  MetaData.MovieDumpName       = DefaultMovieName;
  MetaData.MovieDumpDir        = DefaultMovieDir;
  MetaData.TracerParticleDumpNumber = 0;
  MetaData.TracerParticleDumpName   = DefaultTracerParticleName;
  MetaData.TracerParticleDumpDir    = DefaultTracerParticleDir;
//MetaData.RedshiftDumpNumber  = 0;
  MetaData.RedshiftDumpName    = DefaultRedshiftName;
  MetaData.RedshiftDumpDir     = DefaultRedshiftDir;
 
  MetaData.LocalDir            = NULL;
  MetaData.GlobalDir           = NULL;
 
  for (i = 0; i < MAX_TIME_ACTIONS; i++) {
    TimeActionType[i]      = 0;
    TimeActionRedshift[i]  = -1;
    TimeActionTime[i]      = 0;
    TimeActionParameter[i] = FLOAT_UNDEFINED;
  }
 
  for (i = 0; i < MAX_CUBE_DUMPS; i++) {
    CubeDumps[i] = NULL;
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
#ifdef ISO_GRAV
  MetaData.GravityBoundaryRestart  = 0;
  MetaData.GravityBoundaryName     = "potential_bdry";
  for (dim=0; dim < MAX_DIMENSION; dim++) 
    MetaData.GravityBoundaryFaces[dim] = 0;
#endif

#ifdef RAD_HYDRO
  MetaData.RadHydroParameterFname = NULL;
#endif

  MetaData.ParticleBoundaryType   = periodic;  // only one implemented!
  MetaData.NumberOfParticles      = 0;         // no particles
 
  MetaData.CourantSafetyNumber    = 0.6;
  MetaData.PPMFlatteningParameter = 0;    // off
  MetaData.PPMDiffusionParameter  = 0;    // off
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
  MinimumSubgridEdge        = 4;                 // min for acceptable subgrid
  MaximumSubgridSize        = 2000;              // max for acceptable subgrid
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
    MetaData.NewMovieLeftEdge[dim]  = 0.0;
    MetaData.NewMovieRightEdge[dim] = 1.0;
    PointSourceGravityPosition[dim] = 0.0;
  }
 
  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    StaticRefineRegionLevel[i] = INT_UNDEFINED;
 
  ParallelRootGridIO          = FALSE;
  ParallelParticleIO          = FALSE;
  Unigrid                     = FALSE;
  PartitionNestedGrids        = FALSE;
  ExtractFieldsOnly           = TRUE;
  SRBprefix                   = NULL;
  SRBcwd                      = new char[MAX_LINE_LENGTH];

  SimpleConstantBoundary      = FALSE;

  TracerParticleOn            = 0;
 
  CubeDumpEnabled             = 0;
 
  First_Pass = 0;
 
  UniformGravity              = FALSE;             // off
  UniformGravityDirection     = 0;                 // x-direction
  UniformGravityConstant      = 1.0;
 
  PointSourceGravity           = FALSE;             // off
  PointSourceGravityConstant   = 1.0;
  PointSourceGravityCoreRadius = 0.0;

#ifdef ISO_GRAV
  AverageDensity              = 0.0;               // off
#endif
 
  SelfGravity                 = FALSE;             // off
  CopyGravPotential           = FALSE;             // off
  GravitationalConstant       = 4*Pi;              // G = 1
  S2ParticleSize              = 3.0;               // ~3 is reasonable
  GravityResolution           = 1.0;               // equivalent to grid
  ComputePotential            = FALSE;
  WritePotential              = FALSE;
  BaryonSelfGravityApproximation = TRUE;           // less accurate but faster
 
  GreensFunctionMaxNumber     = 1;                 // only one at a time
  GreensFunctionMaxSize       = 1;                 // not used yet
 
  DualEnergyFormalism         = FALSE;             // off
  DualEnergyFormalismEta1     = 0.001;             // typical 0.001
  DualEnergyFormalismEta2     = 0.1;               // 0.08-0.1
  ParticleCourantSafetyNumber = 0.5;
  RandomForcing               = FALSE;             // off //AK
  RandomForcingEdot           = -1.0;              //AK
  RandomForcingMachNumber     = 0.0;               //AK
#ifdef PPML
  RandomForcingNumberOfFields = 0;                 //dcc.  MUST be set, if RandomForcing is on.
#endif //PPML
  RadiativeCooling            = FALSE;             // off
  RadiationHydrodynamics      = FALSE;             // off
  GadgetEquilibriumCooling    = FALSE;             // off
  MultiSpecies                = FALSE;             // off
  RadiationFieldType          = 0;
  RadiationFieldLevelRecompute = 0;
  AdjustUVBackground          = 1;
  SetUVBAmplitude             = 1.0;
  SetHeIIHeatingScale         = 1.8;
  CoolData.alpha0             = 1.5;               // radiation spectral slope
  CoolData.f3                 = 1.0e-21;           // radiation normalization
  CoolData.ParameterFilename = NULL;
 
  ZEUSLinearArtificialViscosity    = 0.0;
  ZEUSQuadraticArtificialViscosity = 2.0;
  UseMinimumPressureSupport        = FALSE;
  MinimumPressureSupportParameter  = 100.0;
 
  MinimumSlopeForRefinement        = 0.3;          // 30% change in value
  MinimumShearForRefinement        = 1.0;          //AK
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
  MultiMetals                      = FALSE;
  NumberOfParticleAttributes       = INT_UNDEFINED;
 
  for (i = 0; i<MAX_MOVIE_FIELDS; i++)
    MovieDataField[i] = INT_UNDEFINED;
  MovieSkipTimestep = INT_UNDEFINED;
  NewMovieName = DefaultNewMovieName;
  NewMovieDumpNumber = 0;
  NewMovieEntries = 0;
  MovieEntriesPP = new long[NumberOfProcessors];
  MaxMovieFilenum = 0;
  for (i = 0; i<NumberOfProcessors; i++)
    MovieEntriesPP[i] = 0;
  NewMovieParticleOn = FALSE;

  ran1_init = 0;

  return SUCCESS;
}

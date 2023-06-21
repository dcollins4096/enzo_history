/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/    This is the global data, which should be held to a minimum.  Any changes
/    in this file require changes in: WriteGlobalData,
/    ReadGlobalData and InitializeNew.  
/    This file is dual-purposed:
/        1) read with    DEFINE_STORAGE defined for the (single) definition
/        2) read without DEFINE_STORAGE defined for external linkage
/
************************************************************************/
#include <stdio.h>
#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* debugging flag */

EXTERN int debug;

/* Problem: 00 = None                    01 = ShockTube
            02 = WavePool                03 = ShockPool  
	    04 = Double-Mach reflection  05 = ShockInABox
	    20 = 1D Zeldovich Pancake    21 = 1D pressureless collapse
	    22 = Adiabatic expansion     23 = TestGravity
            24 = Spherical infall        25 = TestGravitySphere
	    26 = GravityEquilibriumTest  27 = CollapseTest
	    28 = TestGravityMotion
	    30 = Cosmology simulation
	                                                                  */
EXTERN int ProblemType;

/* Hydrodynamics method:
       0 - PPM_DE      1 - PPM_LR (not working)    2 - ZEUS        */

EXTERN hydro_method HydroMethod;

/* Large and small numbers (i.e. compared to any real quantity).  This may
   be machine and problem dependant. */

EXTERN float huge_number, tiny_number;

/* Gamma: Ideal gas law constant. */

EXTERN float Gamma;

/* Flag indicating if the gas is pressureless. */

EXTERN int PressureFree;

/* Factor to refine by */

EXTERN int RefineBy;

/* Maximum refinement level (0 = topgrid). */

EXTERN int MaximumRefinementLevel;
EXTERN int MaximumGravityRefinementLevel;
EXTERN int MaximumParticleRefinementLevel;

/* Cell Flagging method:  0 = None
                          1 = FlagCellsToBeRefinedBySlope
			  2 = FlagCellsToBeRefinedByMass (baryon only)
			  3 = FlagCellsToBeRefinedByShocks
			  4 = FlagCellsToBeRefinedByMass (particles only)
	     (disabled)	  5 = FlagCellsToBeRefinedByOverdensity (baryon only)
			  6 = FlagCellsToBeRefinedByJeansLength
			  7 = FlagCellsToBeRefinedByCoolingTime
			  8 = FlagCellsToBeRefinedByMustRefineParticles */

EXTERN int CellFlaggingMethod[MAX_FLAGGING_METHODS];

/* Flag indicating if the flux correction should be applied. */

EXTERN int FluxCorrection;

/* This specifies the interpolation method (see typedefs.h). */

EXTERN interpolation_type InterpolationMethod;
EXTERN int ConservativeInterpolation;

/* This is the minimum efficiency of combined grid needs to achieve in
   order to be considered better than the two grids from which it formed. */

EXTERN float MinimumEfficiency;

/* The number of zones that will be refined around each flagged zone. */

EXTERN int NumberOfBufferZones;

/* The left and right boundaries of the entire computational domain. */

EXTERN FLOAT DomainLeftEdge[MAX_DIMENSION], DomainRightEdge[MAX_DIMENSION];

/* Velocity of entire computational domain. */

EXTERN float GridVelocity[MAX_DIMENSION];

/* HDF names for labels and scales. */

EXTERN char *DimUnits[MAX_DIMENSION], *DimLabels[MAX_DIMENSION];
EXTERN char *DataLabel[MAX_NUMBER_OF_BARYON_FIELDS];
EXTERN char *DataUnits[MAX_NUMBER_OF_BARYON_FIELDS];

/* Region in which refinement is allowed (in problem space). */

EXTERN FLOAT RefineRegionLeftEdge[MAX_DIMENSION], 
             RefineRegionRightEdge[MAX_DIMENSION];

/* Uniform gravity: on/off flag, direction, and strength. */

EXTERN int UniformGravity, UniformGravityDirection;
EXTERN float UniformGravityConstant;

/* point source gravity: on/off flag position, and strength. */

EXTERN int PointSourceGravity;
EXTERN FLOAT PointSourceGravityPosition[MAX_DIMENSION];
EXTERN double PointSourceGravityConstant;
EXTERN double PointSourceGravityCoreRadius;

/* SelfGravity (TRUE or FALSE) */

EXTERN int SelfGravity;

/* Flag indicating whether or not to use the baryon self-gravity approximation
   (subgrid cells influence are approximated by their projection to the
   current grid). */

EXTERN int BaryonSelfGravityApproximation;

/* Coefficient in front of source term in Poisson's equations.
   (i.e. Del^phi = GravitationConstant * density, usually 4*Pi*G). */

EXTERN float GravitationalConstant;

/* S2 Particle size in top grid cell units (usually around 3).  The S2
   particle is S(r) = A*(a/2-r) (if r < a/2, 0 otherwise).  The constant
   A depends on the dimension: 1D) 4/a^2,  2D) 24/(Pi*a^3)  3D) 48/(Pi*a^3). */

EXTERN float S2ParticleSize;

/* Gravity resolution factor is a float indicating the comparative resolution
   of the gravitational computation compared to the grid (1-2 or so). */

EXTERN float GravityResolution;

/* Flag to indicate if gravitational potential field should be computed
   and stored. */

EXTERN int ComputePotential;

/* Maximum number of GreensFunctions that will be stored in any time.
   This number must be less than MAX_NUMBER_OF_GREENS_FUNCTIONS. */

EXTERN int GreensFunctionMaxNumber;

/* Maximum number of words associated with GreensFunction storage
   (Not currently implemented). */

EXTERN int GreensFunctionMaxSize;

/* Dual energy formalism (TRUE or FALSE). */

EXTERN int DualEnergyFormalism;

/* Two parameters for the dual energy formalism. */

EXTERN float DualEnergyFormalismEta1;
EXTERN float DualEnergyFormalismEta2;

/* This is the particle equivalent of the Courant factor.  It is the maximum
   number of cells a particle is allowed to travel in a single timestep. */

EXTERN float ParticleCourantSafetyNumber;

/* Radiative cooling on/off flag and associated data. */

EXTERN int RadiativeCooling;
EXTERN CoolDataType CoolData;

/* Multi-species rate equation flag and associated data. */

EXTERN int MultiSpecies;
EXTERN RateDataType RateData;

/* Type of radiation field. 
   0 - none,                    1 - Haardt & Madau alpha=-1.5
   2 - H&M alpha = -1.8       
   10 - homogenous internal radiation field (a la Renyue's work) */

EXTERN int RadiationFieldType;
EXTERN RadiationFieldDataType RadiationData;
EXTERN int RadiationFieldLevelRecompute;

/* ZEUS Hydro artificial viscosity parameters (C1, C2 of Stone & Norman). */

EXTERN float ZEUSLinearArtificialViscosity;
EXTERN float ZEUSQuadraticArtificialViscosity;

/* Parameters for MinimumPressureSupport. */

EXTERN int UseMinimumPressureSupport;
EXTERN float MinimumPressureSupportParameter;

/* Parameters for statically refined regions. */

EXTERN FLOAT StaticRefineRegionLeftEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
EXTERN FLOAT StaticRefineRegionRightEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
EXTERN int   StaticRefineRegionLevel[MAX_STATIC_REGIONS];

/* Processor identifier for this thread/processor */

EXTERN int MyProcessorNumber;
EXTERN int NumberOfProcessors;
EXTERN float CommunicationTime;
EXTERN double PerformanceTimers[MAX_PERFORMANCE_TIMERS];

/* Parameter to indicate if top grid should do parallel IO
   (currently only works for ProblemType == 30). */

EXTERN int ParallelRootGridIO;

/************************************************/
/* Global data for specific problems or methods */
/************************************************/

/* For CellFlaggingMethod = 1,
   The minimum relative slope (da/dx over a) required for refinement. */

EXTERN float MinimumSlopeForRefinement;

/* For CellFlaggingMethod = 2,
   The minimum refined mass for the ByMass refining scheme
   (Usually, user sets OverDensity and code sets MinimumMass but this can be
    overridden by directely setting MinimumMass). 
   The LevelExponent is used to change the minimum mass with level,
   the formula is MinimumMassForRefinement*pow(RefineBy, level*LevelExponent)*/

EXTERN float MinimumOverDensityForRefinement[MAX_FLAGGING_METHODS];
EXTERN float MinimumMassForRefinement[MAX_FLAGGING_METHODS];
EXTERN float MinimumMassForRefinementLevelExponent[MAX_FLAGGING_METHODS];

/* For CellFlaggingMethod = 3,
   The minimum pressure jump required to be a shock.
   The minimum internal/total energy ratio for a shock. */

EXTERN float MinimumPressureJumpForRefinement, MinimumEnergyRatioForRefinement;

/* For CellFlaggingMethod = 6,
   The number of cells by which the Jeans length should be resolved. */

EXTERN float RefineByJeansLengthSafetyFactor;

/* For CellFlaggingMethod = 8,
   The level to which the must refine particles apply */

EXTERN int   MustRefineParticlesRefineToLevel;

/* A boolean flag indicating if we are using coordinate comoving with the
   expansion of the universe. */

EXTERN int   ComovingCoordinates;

/* A flag indicating if we are using star particles. */

EXTERN int   StarParticleCreation;
EXTERN int   StarParticleFeedback;
EXTERN int   NumberOfParticleAttributes;

/* Parameters governing certain time or redshift-dependent actions. */

EXTERN int   TimeActionType[MAX_TIME_ACTIONS];
EXTERN FLOAT TimeActionTime[MAX_TIME_ACTIONS];
EXTERN FLOAT TimeActionRedshift[MAX_TIME_ACTIONS];
EXTERN float TimeActionParameter[MAX_TIME_ACTIONS];

/* Parameters governing whether tracer particles are on or off. */

EXTERN int   TracerParticleOn;

/* Zhiling Lan's modified code */

#ifdef MPI_INSTRUMENTATION
EXTERN float GlobalCommunication;
EXTERN float RecvComm;
EXTERN double timer[40];
EXTERN int counter[40];
EXTERN char name[20];
EXTERN FILE *filePtr;
EXTERN double starttime, endtime;
EXTERN double Start_Wall_Time, End_Wall_Time, WallTime;
EXTERN int flagging_count, in_count, out_count, moving_count;
EXTERN float flagging_pct, moving_pct;
#endif /* MPI_INSTRUMENTATION */




#ifndef __macros_and_parameters_h_
#define __macros_and_parameters_h_
/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

/* Modifiable Parameters */

#define MAX_NUMBER_OF_BARYON_FIELDS        20  /* must be at least 7 */

#define MAX_NUMBER_OF_SUBGRIDS          50000

#define MAX_DEPTH_OF_HIERARCHY             50

#define MAX_LINE_LENGTH                   512

#define DEFAULT_GHOST_ZONES                 3  /* at least 3 */

#define MAX_NUMBER_OF_OUTPUT_REDSHIFTS   1000

#define GRAVITY_BUFFER_SIZE                 3

#define MAX_FLAGGING_METHODS                5

#define MAX_STATIC_REGIONS               1000

#define MAX_NUMBER_OF_PARTICLE_ATTRIBUTES   5

#define MAX_TIME_ACTIONS                   10

#define ROOT_PROCESSOR                      0

#define MAX_PERFORMANCE_TIMERS             32

#define VERSION                           1.3  /* current version number */

/* File interface to use in Grid_WriteGrid.C
   WRITEGRID_HDF_DFSD is a bit slow and inefficient, but is more widely used,
   WRITEGRID_HDF_SD is more stable(?) but may confuse some readers 
   WRITEGRID_RAW is probably stable but highly unportable.
   WRITEGRID_FLEXIO is John Shalf's FlexIO format. 
   WRITEGRID_FORTRAN is Cen's fortran format (ReadGrid doesn't read this!). */

/* #define WRITEGRID_HDF_DFSD */
#define WRITEGRID_HDF_SD  /* standard */
/* #define WRITEGRID_RAW  */
/* #define WRITEGRID_FLEXIO */
/* #define WRITEGRID_FORTRAN */

/* Unmodifiable Parameters */

#define MAX_DIMENSION                       3  /* must be 3! */


/* Fortran name generator (cpp blues) */

#if defined(SUN_OLD)
#define FORTRAN_NAME(NAME) NAME/**/_
#endif

#if defined(IRIS4) || defined(CONVEX) || defined(ALPHA) || defined(SUN) || defined(LINUX)
#define FORTRAN_NAME(NAME) NAME##_
#endif

#if defined(SPP) || defined(SP2)
#define FORTRAN_NAME(NAME) NAME
#endif

/* Precision-related definitions. */

//typedef long long long_int;  /* now done below */

/* Define the long double type. */

typedef long double long_double;

#ifdef r4
#define FLOAT float
#define FSYM "f"
#define GOUTSYM ".7g"
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_FLOAT
#endif /* USE_MPI */
typedef int long_int;
#endif /* r4 */

#ifdef r8
#define FLOAT double
#define FSYM "lf"
#define GOUTSYM ".14g"
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_DOUBLE
#endif /* USE_MPI */
typedef long long long_int;
#endif /* r8 */

#ifdef r16

#ifdef FSYM_lf
#define FSYM "lf"
#define GOUTSYM ".21G"
#else
#define FSYM "Lf"
#define GOUTSYM ".21Lg"
#endif

#define FLOAT long_double
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_LONG_DOUBLE
#endif /* USE_MPI */
typedef long long long_int;
#endif /* r16 */

/* Standard definitions (well, fairly standard) */

#ifndef NULL
#define NULL      0
#endif

#ifdef FAIL
#undef FAIL
#endif
#define FAIL      0
#define SUCCESS   1

#ifndef FALSE
#define FALSE     0
#define TRUE      1
#endif

/* Not-so standard definitions */
#ifndef HDF_FAIL
#define HDF_FAIL -1
#endif

#define FLOAT_UNDEFINED  -99999.0
#define INT_UNDEFINED    -99999

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define nint(A) int((A) + 0.5*sign(A))
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )

/* Definitions for FastFourierTransform and related routines */

#define FFT_FORWARD  +1
#define FFT_INVERSE  -1
#define REAL_TO_COMPLEX    0
#define COMPLEX_TO_COMPLEX 1

/* Definitions for grid::RestoreEnergyConsistency */

#define ENTIRE_REGION  0
#define ONLY_BOUNDARY  1

/* Definitions for grid::ZeroSolutionUnderSubgrid */

#define ZERO_ALL_FIELDS          0
#define ZERO_UNDER_SUBGRID_FIELD 1

/* Definitions for grid::CommunicationSend/ReceiveRegion and 
   grid::DepositPositions */

#define MASS_FLAGGING_FIELD              -6
#define ACCELERATION_FIELDS              -5
#define POTENTIAL_FIELD                  -4
#define GRAVITATING_MASS_FIELD           -3
#define GRAVITATING_MASS_FIELD_PARTICLES -2
#define ALL_FIELDS   -1

#define NEW_AND_OLD   0
#define NEW_ONLY      1
#define OLD_ONLY      2

/* Definitions for grid::ComputeAccelerationField */

#define PARTICLES  0
#define GRIDS      1
#define ZEUS_GRIDS 2

/* Definitions for CommunicationTranspose */

#define NORMAL_ORDER      0
#define TRANSPOSE_FORWARD 1
#define TRANSPOSE_REVERSE 2

/* Definitions for CommunicationTransferParticles */

#define COPY_IN   0
#define COPY_OUT  1

/* Definitions for CommunicationDirection */

#define COMMUNICATION_SEND_RECEIVE 0
#define COMMUNICATION_POST_RECEIVE 1
#define COMMUNICATION_RECEIVE      2
#define COMMUNICATION_SEND         3

/* MPI Tags */

#define MPI_TRANSPOSE_TAG 10
#define MPI_SENDREGION_TAG 11
#define MPI_FLUX_TAG 12
#define MPI_TRANSFERPARTICLE_TAG 13
#define MPI_SENDPARTFLOAT_TAG 14
#define MPI_SENDPARTINT_TAG 15

/* Definitions for CommunicationBufferedSend. */

#define BUFFER_IN_PLACE -1

/* Particle types (note: gas is a conceptual type) */

#define NUM_PARTICLE_TYPES 5

#define PARTICLE_TYPE_GAS          0
#define PARTICLE_TYPE_DARK_MATTER  1
#define PARTICLE_TYPE_STAR         2
#define PARTICLE_TYPE_TRACER       3
#define PARTICLE_TYPE_MUST_REFINE  4

#ifdef USE_MPI
#define NO_MPI_INSTRUMENTATION
#endif /* USE_MPI */

#endif

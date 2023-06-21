#ifndef __macros_and_parameters_h_
#define __macros_and_parameters_h_
/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

/* Modifiable Parameters */

#define MAX_NUMBER_OF_BARYON_FIELDS        20  /* must be at least 6 */

#define MAX_NUMBER_OF_SUBGRIDS          50000

#define MAX_DEPTH_OF_HIERARCHY             40

#define MAX_LINE_LENGTH                   512

#define DEFAULT_GHOST_ZONES                 3  /* at least 3 */

#define MAX_NUMBER_OF_OUTPUT_REDSHIFTS    500

#define GRAVITY_BUFFER_SIZE                 3

#define MAX_FLAGGING_METHODS                5

#define MAX_STATIC_REGIONS                 20

#define MAX_NUMBER_OF_PARTICLE_ATTRIBUTES  10

#define MAX_TIME_ACTIONS                    5

#define ROOT_PROCESSOR                      0

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

#if defined(IRIS4) || defined(CONVEX) || defined(COMPAQ) || defined(SUN) || defined(LINUX) || defined(IA64)
#define FORTRAN_NAME(NAME) NAME##_
#endif

#if defined(SPP) || defined(SP2) || defined(IA65)
#define FORTRAN_NAME(NAME) NAME
#endif

/* Precision-related definitions. */

/* Define the long int type. */

typedef long long long_int;

/* Define the long double type. */

typedef long double long_double;

#define io_float float

#ifdef p4
#define FLOAT float
#define PSYM "f"
#define FSYM "f"
#define GOUTSYM ".7g"
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_FLOAT
#endif /* USE_MPI */
#endif /* p4 */

#ifdef p8
#define FLOAT double
#define PSYM "lf"
#define GOUTSYM ".14g"
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_DOUBLE
#endif /* USE_MPI */
#ifdef r4
#define FSYM "f"
#endif /* r4 */
#ifdef r8
#define float32 TEMP_HOLD_NAME
#define float double
#define TEMP_HOLD_NAME float32
#define FSYM "lf"
#endif /* r8 */
#endif /* p8 */

#ifdef p16
#define FLOAT long_double
#define PSYM "Lf"
#define GOUTSYM ".21Lg"
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_LONG_DOUBLE
#endif /* USE_MPI */
#ifdef r4
#define FSYM "f"
#endif /* r4 */
#ifdef r8
#define float32 TEMP_HOLD_NAME
#define float double
#define TEMP_HOLD_NAME float32
#define FSYM "lf"
#endif /* r8 */
#endif /* p16 */

#ifdef HDF5_BE
#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_R4 H5T_IEEE_F32BE
#define HDF5_FILE_R8 H5T_IEEE_F64BE
#define HDF5_FILE_B8 H5T_STD_B8BE
#endif /* HDF5_BE */

#ifdef HDF5_LE
#define HDF5_FILE_I4 H5T_STD_I32LE
#define HDF5_FILE_R4 H5T_IEEE_F32LE
#define HDF5_FILE_R8 H5T_IEEE_F64LE
#define HDF5_FILE_B8 H5T_STD_B8LE
#endif /* HDF5_LE */

#define HDF5_I4 H5T_NATIVE_INT
#define HDF5_R4 H5T_NATIVE_FLOAT
#define HDF5_R8 H5T_NATIVE_DOUBLE
#define HDF5_R16 H5T_NATIVE_LDOUBLE

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
#define POW(X,Y) pow((double) (X), (double) (Y))

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
#define COMMUNICATION_RECEIVE      1
#define COMMUNICATION_SEND         2

/* MPI Tags */

#define MPI_TRANSPOSE_TAG 10
#define MPI_SENDREGION_TAG 11
#define MPI_FLUX_TAG 12
#define MPI_TRANSFERPARTICLE_TAG 13
#define MPI_SENDPARTFLOAT_TAG 14
#define MPI_SENDPARTINT_TAG 15

/* Definitions for CommunicationBufferedSend. */

#define BUFFER_IN_PLACE -1

#ifdef USE_MPI
#define MPI_INSTRUMENTATION
#endif /* USE_MPI */

#endif

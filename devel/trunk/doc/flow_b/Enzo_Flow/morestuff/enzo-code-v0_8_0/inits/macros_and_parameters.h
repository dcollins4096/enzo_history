/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

/* Modifiable Parameters */

#define MAX_LINE_LENGTH                   512

#define VERSION                           2.0  /* current version number */

#define MAX_SPECIES                        10

#define MAX_PS_NAMES                        2

#define NO_SHIFT_FOR_LARS

#define MAX_DIMENSION 3 /* don't change this */

/* Fortran name generator (cpp blues) */

#if defined(SUN_OLD)
#define FORTRAN_NAME(NAME) NAME/**/_
#endif

#if defined(IRIS4) || defined(CONVEX) || defined(COMPAQ) || defined(SUN) || defined(I686) || defined(IA64) || defined(GNU)
#define FORTRAN_NAME(NAME) NAME##_
#endif

#if defined(SPP) || defined(SP2)
#define FORTRAN_NAME(NAME) NAME
#endif

/* Precision-related definitions. */

/* Define the long int type. */

typedef long long long_int;

/* Define the long double type. */

typedef long double long_double;

#ifdef r4
#define FLOAT float
#define FSYM "f"
#define GOUTSYM ".7g"
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_FLOAT
#endif /* USE_MPI */
#endif /* r4 */

#ifdef r8
#define FLOAT double
#define FSYM "lf"
#define GOUTSYM ".14g"
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_DOUBLE
#endif /* USE_MPI */
#endif /* r8 */

#ifdef r16
#define FLOAT long_double
#define FSYM "Lf"
#define GOUTSYM ".21Lg"
#ifdef USE_MPI
#define MY_MPIFLOAT MPI_LONG_DOUBLE
#endif /* USE_MPI */
#endif /* r16 */

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

#define HDF_FAIL -1

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

/* RH debug */

#define RH_D TRUE

/* RH memory trace */

#define RH_M TRUE

/* RH flow trace */

#define RH_F TRUE

/* Fortran 77, 90, or 95-specific defines */

#ifdef USE_F95
#define _IMAG imag
#else
#define _IMAG aimag
#endif


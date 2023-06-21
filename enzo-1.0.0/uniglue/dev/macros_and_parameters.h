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
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

/* HDF5 Big Endian (Enzo default) */

#define HDF5_R4 H5T_IEEE_F32BE
#define HDF5_R8 H5T_IEEE_F64BE
#define HDF5_B8 H5T_STD_B8BE
#define HDF5_I4 H5T_STD_I32BE
#define HDF5_I8 H5T_STD_I64BE

/* HDF5 Little Endian (UNTESTED!) */

/* #define HDF5_R4 H5T_IEEE_F32LE */
/* #define HDF5_R8 H5T_IEEE_F64LE */
/* #define HDF5_B8 H5T_STD_B8LE */
/* #define HDF5_I4 H5T_STD_I32LE */
/* #define HDF5_I8 H5T_STD_I64LE */

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define nint(A) int((A) + 0.5*sign(A))
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )


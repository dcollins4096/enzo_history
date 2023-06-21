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
/* Macro definitions for HDF4 compatibility */

#if defined(SP2)
typedef void         *VOIDP;
typedef float        float32;
typedef double       float64;
#endif

#if defined(SUN)
typedef void         *VOIDP;
typedef int          int32;
typedef float        float32;
typedef double       float64;
#endif

#if defined(IRIX) || defined(IRIS4) || defined(COMPAQ) || defined(I686) || defined(IA64) || defined(GNU)
typedef void         *VOIDP;
typedef int          int32;
typedef float        float32;
typedef double       float64;
#endif

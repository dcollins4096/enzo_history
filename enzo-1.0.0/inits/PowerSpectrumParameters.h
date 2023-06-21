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
/* Power Spectrum Parameters */

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* Setable parameters. */

EXTERN float sigma8;
EXTERN float PrimordialIndex;
EXTERN float Gamma;
EXTERN float kcutoff;
EXTERN int   RandomSeed;
EXTERN float kmin;
EXTERN float kmax;
EXTERN int   NumberOfkPoints;
EXTERN int   PowerSpectrumType;

EXTERN char *PSFileName[2];

EXTERN float WDMPartMass;  // WDM particle mass

EXTERN float WDMg_x;  // WDM degrees of freedom

/* Space for pre-computed look-up table. */

EXTERN float *PSLookUpTable[MAX_SPECIES];

/* Internal parameters. */

EXTERN float Normalization;
EXTERN float GrowthFactor;
EXTERN float Redshift;
EXTERN float TophatRadius;

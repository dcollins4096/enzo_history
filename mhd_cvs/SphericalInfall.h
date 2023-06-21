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
/  GLOBAL DATA DECLARATIONS FOR THE SPHERICAL INFALL TEST
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define SIEXTERN
#else /* DEFINE_STORAGE */
# define SIEXTERN extern
#endif /* DEFINE_STORAGE */

/* Flag indicating whether to use fixed perturbation mass. */

SIEXTERN int SphericalInfallFixedAcceleration;

/* Perturbation mass. */

SIEXTERN float SphericalInfallFixedMass;

/* Center of peturbation. */

SIEXTERN FLOAT SphericalInfallCenter[MAX_DIMENSION];

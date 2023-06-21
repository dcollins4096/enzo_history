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
/* Cosmology Parameters */

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* The Hubble constant at z=0 (now) in units of 100 km/s/Mpc. */

EXTERN float HubbleConstantNow;

/* The value of Omega due to non-relativistic particles at z=0. */

EXTERN float OmegaMatterNow;

/* The value of Omega due to lamba (the cosmological constant) at z=0. */

EXTERN float OmegaLambdaNow;

/* The value of Omega due to WDM at z=0 */

EXTERN float OmegaWDMNow;

/* The value of Omega due to HDM (neutrinos) at z=0. */

EXTERN float OmegaHDMNow;

/* The value of Omega due to baryons at z=0. */

EXTERN float OmegaBaryonNow;

/* The comoving size of the simulation box (along the x-dir) in h^{-1} Mpc. */

EXTERN float ComovingBoxSize;

/* The initial redshift. */

EXTERN float InitialRedshift;


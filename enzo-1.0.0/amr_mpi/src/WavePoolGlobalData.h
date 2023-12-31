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
/  GLOBAL DATA DECLARATIONS FOR THE WAVE POOL PROBLEM
/
/  written by: Greg Bryan
/  date:       February, 1994
/  modified1:
/
/  PURPOSE:
/    This is global data that pertains only to the Wave Pool test problem.
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define WPEXTERN
#else /* DEFINE_STORAGE */
# define WPEXTERN extern
#endif /* DEFINE_STORAGE */

/* Angle of the vector of propogation with respect to the x-axis.  */

WPEXTERN float WavePoolAngle;

/* Amplitude of the wave, in terms of maximal density perturbation. */

WPEXTERN float WavePoolAmplitude;

/* Wavelength of one complete wave in problem units. */

WPEXTERN float WavePoolWavelength;

/* Number of full waves to enter */

WPEXTERN float WavePoolNumberOfWaves;

/* Density, pressure and initial velocity of the wave pool. */

WPEXTERN float WavePoolDensity;
WPEXTERN float WavePoolPressure;
WPEXTERN float WavePoolVelocity[MAX_DIMENSION];


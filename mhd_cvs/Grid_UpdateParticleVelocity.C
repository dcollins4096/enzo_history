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
/  GRID CLASS (UPDATE PARTICLE VELOCITY FROM ACCELERATIONS)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

#define VELOCITY_METHOD3

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


int grid::UpdateParticleVelocity(float TimeStep)
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0 || SelfGravity == FALSE) return SUCCESS;

  FLOAT a = 1.0, dadt;
#if defined(VELOCITY_METHOD1) || defined(VELOCITY_METHOD2)
  float VelocityMidStep;
#endif 
  int i;

  /* If using comoving coordinates, divide by a(t) first. */

  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time + TimeStep, &a, &dadt) 
	== FAIL) {
      fprintf(stderr, "Error in CsomologyComputeExpansionFactors.\n");
      return FAIL;
    }

  /* Loop over dimensions. */

  for (int dim = 0; dim < GridRank; dim++) {

    /* Error check. */

    if (ParticleAcceleration[dim] == NULL) {
      fprintf(stderr, "No ParticleAccleration present.\n");
      return FAIL;
    }

    /* Update velocities.  */

    if (ComovingCoordinates) {

      FLOAT coef = 0.5*dadt/a*TimeStep;

      /* If using comoving coordinates, subtract the (time-centered) 
	 drag-like term and add the acceleration. The acceleration has
	 already been divided by a(t). */
      
      for (i = 0; i < NumberOfParticles; i++) {

#ifdef VELOCITY_METHOD1

        /* i) partially time-centered. */

	VelocityMidStep = ParticleVelocity[dim][i] + 
	                  ParticleAcceleration[dim][i]*0.5*TimeStep;

	ParticleVelocity[dim][i] += 
	  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;

#endif /* VELOCITY_METHOD1 */

#ifdef VELOCITY_METHOD2

        /* ii) partially backward. */

	VelocityMidStep = ParticleVelocity[dim][i] ;

	ParticleVelocity[dim][i] += 
	  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;

#endif /* VELOCITY_METHOD2 */

#ifdef VELOCITY_METHOD3

        /* iii) Semi-implicit way */

        ParticleVelocity[dim][i] = ((1.0 - coef)*ParticleVelocity[dim][i] +
                                    ParticleAcceleration[dim][i]*TimeStep)/
                                   (1.0 + coef);

#endif /* VELOCITY_METHOD3 */

      }
    }
    else

      /* Otherwise, just add the acceleration. */

      for (i = 0; i < NumberOfParticles; i++)
	ParticleVelocity[dim][i] += ParticleAcceleration[dim][i] * TimeStep;

  }

  return SUCCESS;
}

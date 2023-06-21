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
/  GRID CLASS (INITIALIZE THE GRID FOR A GRAVITY TEST)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


int grid::TestGravityMotionInitializeGrid(float InitialVelocity)
{
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int dim, i, NumberOfNewParticles = 1;

  /* Set particles. */

  if (NumberOfNewParticles > 0) {

    /* Set number of particles for this grid. */

    NumberOfParticles = NumberOfNewParticles;

    /* Allocate space. */

    ParticleMass = new float[NumberOfParticles];
    ParticleNumber = new int[NumberOfParticles];
    for (dim = 0; dim < GridRank; dim++) {
      ParticlePosition[dim] = new FLOAT[NumberOfParticles];
      ParticleVelocity[dim] = new float[NumberOfParticles];
    }

    /* Set up particle positon & velocity. */

    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < NumberOfParticles; i++) {
	ParticlePosition[dim][i] = 0.5*(DomainLeftEdge[dim] +
					DomainRightEdge[dim]);
	ParticleVelocity[dim][i] = 0;
	if (dim == 0)
	  ParticleVelocity[dim][i] = InitialVelocity;
	ParticleMass[i]          = 1.0;
	ParticleNumber[i]        = i;
      }

  }

  return SUCCESS;
}

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
/  GRID CLASS (ALLOCATES AND CLEARS THE ACCELERATION FIELD FOR PARTICLESS)
/
/  written by: Greg Bryan
/  date:       March, 1995
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



int grid::ClearParticleAccelerations()
{

  int i, dim;

  if (NumberOfParticles > 0)

    /* Loop over active dimension */

    for (dim = 0; dim < GridRank+ComputePotential; dim++) {

      /* Error check. */

      if (ParticleAcceleration[dim] != NULL)
	fprintf(stderr, "ClearParticleAccelerations: Field not NULL.\n");

      /* Allocate accleration field. */

      ParticleAcceleration[dim] = new float[NumberOfParticles];

      /* Clear it. */

      for (i = 0; i < NumberOfParticles; i++)
	ParticleAcceleration[dim][i] = 0.0;

    }

  return SUCCESS;
}


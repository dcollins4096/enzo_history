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
/  GRID CLASS (INTERPOLATE PARTICLE POSITIONS FROM THE ACCELERATION FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE: THIS ROUTINE DOES NOT TRANSFER DATA BETWEEN PROCESSORS!
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


int grid::InterpolateParticlePositions(grid *FromGrid, int DifferenceType)
{

//  if (this != FromGrid)
//    return FAIL;

  FLOAT HoldLeftEdge[MAX_DIMENSION];

  /* Loop over all active dimensions */

  if (NumberOfParticles > 0)
    for (int dim = 0; dim < GridRank+ComputePotential; dim++) {

      /* Adjust the grid position if the acceleration is face-centered. */

      if ((DifferenceType == PARTICLES || DifferenceType == ZEUS_GRIDS) &&
	  dim != GridRank) {
	HoldLeftEdge[dim] = FromGrid->CellLeftEdge[dim][0];
	FromGrid->CellLeftEdge[dim][0] -= 0.5*FromGrid->CellWidth[dim][0];
      }

      if (FromGrid->InterpolatePositions(ParticlePosition, dim,
					 ParticleAcceleration[dim],
					 NumberOfParticles) == FAIL) {
	fprintf(stderr, "Error in grid->InterpolatePositions.\n");
	return FAIL;
      }

      /* Adjust back. */

      if ((DifferenceType == PARTICLES || DifferenceType == ZEUS_GRIDS) &&
	  dim != GridRank)
	FromGrid->CellLeftEdge[dim][0] = HoldLeftEdge[dim];
    }

  return SUCCESS;
}

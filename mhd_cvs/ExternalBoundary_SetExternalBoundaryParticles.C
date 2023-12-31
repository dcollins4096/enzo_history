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
/  EXTERNAL BOUNDARY CLASS (HANDLE PARTICLE EXTERNAL BOUNDARY CONDITIONS)
/
/  written by: Greg Bryan
/  date:       June, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
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

//
//
int ExternalBoundary::SetExternalBoundaryParticles(int FieldRank, 
						  int NumberOfParticles,
						  FLOAT *Position[],
						  float *Velocity[])
{

  /* declarations */

  int dim, i;

  /* Error check: grid ranks. */

  if (FieldRank != BoundaryRank) {
    fprintf(stderr, "FieldRank(%d) != BoundaryRank(%d).\n",
            FieldRank, BoundaryRank);
    return FAIL;
  }

  /* Error check: allowed boundary types. */

  if (ParticleBoundaryType != periodic) {
    fprintf(stderr, "only periodic particle boundary conditions supported.\n");
    return FAIL;
  }

  /* PERIODIC BOUNDARY: wrap particles in each dimension. */

  if (ParticleBoundaryType == periodic && NumberOfProcessors == 1)

    for (dim = 0; dim < FieldRank; dim++)

      for (i = 0; i < NumberOfParticles; i++) {

	if (Position[dim][i] < DomainLeftEdge[dim])
	  Position[dim][i] += DomainRightEdge[dim] - DomainLeftEdge[dim];

	if (Position[dim][i] > DomainRightEdge[dim])
	  Position[dim][i] -= DomainRightEdge[dim] - DomainLeftEdge[dim];

      }	
  
  return SUCCESS;

}

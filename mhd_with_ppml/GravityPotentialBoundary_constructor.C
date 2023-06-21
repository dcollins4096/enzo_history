/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Daniel R. Reynolds                                         *
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRAVITY POTENTIAL BOUNDARY CLASS (CONSTRUCTOR)
/
/  written by: Daniel R. Reynolds
/  date:       September, 2005
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "GravityPotentialBoundary.h"

//
// Constructor
//

#ifdef ISO_GRAV
GravityPotentialBoundary::GravityPotentialBoundary()
{

  /* declarations */
  int dim;

  /* Set scalars to zero. */
  BoundaryRank = 0;
  GravityCellSize = 0.0;

  /* Set flag determining preparedness for use */
  ReadyToGo = FALSE;

  /* Initialize BoundaryType to periodic,  clear
     BoundaryValues pointers, set boundary ownership to false,
     initialize subdomain gravity information to zeros. */
  for (dim=0; dim<MAX_DIMENSION; dim++) {
    LocalDimension[dim] = 0;
    GlobalDimension[dim] = 0;
    BoundaryType[dim] = 0;
    BoundaryValues[dim][0] = NULL;
    BoundaryValues[dim][1] = NULL;
    OwnBoundary[dim][0] = false;
    OwnBoundary[dim][1] = false;
    GravityLeftEdge[dim] = 0.0;
    GravityRightEdge[dim] = 0.0;
    LeftEdgeIndex[dim] = 0;
    RightEdgeIndex[dim] = 0;
  }

}
#endif

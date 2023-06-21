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
/  EXTERNAL BOUNDARY CLASS (INITIALIZE BOUNDARY FACE TO A CONSTANT VALUE)
/
/  written by: Greg Bryan
/  date:       November, 1994
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

// Set external boundaries to a face-constant value.
//   (Note: we do not need BoundaryValue if the boundary type is inflow)

int ExternalBoundary::InitializeExternalBoundaryFace(int dim, 
					     boundary_type LeftBoundaryType,
					     boundary_type RightBoundaryType,
					     float LeftBoundaryValue[],
					     float RightBoundaryValue[])
{
  int field, index;

  /* Error check */

  if (dim > BoundaryRank) {
    fprintf(stderr, "Dimension %d > BoundaryRank %d.\n", dim, BoundaryRank);
    return FAIL;
  }

  /* compute size of entire mesh */

  int size = 1;
  for (int i = 0; i < BoundaryRank; i++)
    size = size*BoundaryDimension[i];

  /* set BoundaryType faces to a constant */

  if (BoundaryDimension[dim] != 1)
    for (field = 0; field < NumberOfBaryonFields; field++) {
      BoundaryType[field][dim][0] = 
	new boundary_type[size/BoundaryDimension[dim]];
      BoundaryType[field][dim][1] = 
	new boundary_type[size/BoundaryDimension[dim]];
      for (index = 0; index < size/BoundaryDimension[dim]; index++) {
	BoundaryType[field][dim][0][index] = LeftBoundaryType;
	BoundaryType[field][dim][1][index] = RightBoundaryType;
      }
    }

  /* If required, set BoundaryType faces to a constant (usually inflow) */

  if (BoundaryDimension[dim] != 1)

    for (field = 0; field < NumberOfBaryonFields; field++) {

      if (LeftBoundaryType == inflow) {
	BoundaryValue[field][dim][0] = new float[size/BoundaryDimension[dim]];
	for (index = 0; index < size/BoundaryDimension[dim]; index++)
	  BoundaryValue[field][dim][0][index] = LeftBoundaryValue[field];
      }

      if (RightBoundaryType == inflow) {
	BoundaryValue[field][dim][1] = new float[size/BoundaryDimension[dim]];
	for (index = 0; index < size/BoundaryDimension[dim]; index++)
	  BoundaryValue[field][dim][1][index] = RightBoundaryValue[field];
      }

    } // end of for loop over fields

  return SUCCESS;

}

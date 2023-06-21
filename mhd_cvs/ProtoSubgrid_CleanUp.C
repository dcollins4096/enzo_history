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
/  PROTOSUBGRID CLASS (FREE UNNECESSARY MEMORY ALLOCATIONS)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
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

void WriteListOfInts(FILE *fptr, int N, int nums[]);

int ProtoSubgrid::CleanUp()
{

  delete GridFlaggingField;
  GridFlaggingField = NULL;

  /* Delete signatures unless dim=1 in which case Signature[0] = GFF */

  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    delete Signature[dim];
    Signature[dim] = NULL;
  }

  /* Set the dimension according to the refinement and ghost zones. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    GridDimension[dim] = GridDimension[dim]*RefineBy + 2*DEFAULT_GHOST_ZONES;
  }

  /* Compute efficiency & report. */

  if (debug) {
    printf("ProtoSubgrid: efficiency = %6.1f%% (%d/%d) dims=", 
	   float(NumberFlagged)/float(size)*100.0, NumberFlagged, size);
    for (int i = 0; i < GridRank; i++)
      printf("%d ", (GridDimension[i] - 2*DEFAULT_GHOST_ZONES)/RefineBy);
    printf("\n");
  }


  return SUCCESS;
}

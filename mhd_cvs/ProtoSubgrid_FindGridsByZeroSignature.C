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
/  PROTOSUBGRID CLASS (DIVIDE THIS GRID BY FINDING ZEROS IN THE SIGNATURE)
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


int ProtoSubgrid::FindGridsByZeroSignature(int dim, int &NumberOfNewGrids, 
				     int GridEnds[MAX_NUMBER_OF_SUBGRIDS][2])
{
  /* Error check */

  if (dim >= GridRank) {
    fprintf(stderr, "Passed dim(%d) > GridRank(%d)\n", dim, GridRank);
    return FAIL;
  }

  if (Signature[dim] == NULL) {
    fprintf(stderr, "Signature %d not yet computed.\n", dim);
    return FAIL;
  }

  /* Initialize */

  int i = 0;
  NumberOfNewGrids = 0;

  /* Loop over signature. */

  while (i < GridDimension[dim]) {

    /* Look for the start of a new subgrid. */

    if (Signature[dim][i] != 0) {
      GridEnds[NumberOfNewGrids][0] = StartIndex[dim] + i;

      /* Now find the end of the subgrid. */

      while (i < GridDimension[dim] && Signature[dim][i] != 0)
	i++;
      GridEnds[NumberOfNewGrids++][1] = StartIndex[dim] + i-1;
    }

    /* Next zone in signature. */

    i++;
  }

  return SUCCESS;
}

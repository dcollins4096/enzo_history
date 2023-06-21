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
/  PROTOSUBGRID CLASS (COMPUTES A 1D SIGNATURE (PROJECTION))
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

/* function prototypes */

extern "C" void FORTRAN_NAME(project)(int *rank, int *i1, int *i2, int *i3,
                                      int *iline, int *pdim, int *field,
                                      int *line);

int ProtoSubgrid::ComputeSignature(int dim)
{

  /* Error check. */

  if (dim >= GridRank) {
    fprintf(stderr, "Project: dim = %d > GridRank = %d\n", dim, GridRank);
    return FAIL;
  }

  /* Already done? */

  if (Signature[dim] != NULL)
    return SUCCESS;

  /* Allocate space. */

  Signature[dim] = new int[GridDimension[dim]];

  /* Perform projection (and clear line)
     For 1D grids, no projection necessary, just copy. */

  if (GridRank == 1)
    for (int i = 0; i < GridDimension[dim]; i++)
      Signature[dim][i] = GridFlaggingField[i];
  else
    FORTRAN_NAME(project)(&GridRank, GridDimension, GridDimension+1,
                          GridDimension+2, &GridDimension[dim], &dim,
                          GridFlaggingField, Signature[dim]);
  
  /*  if (debug) {
      printf ("sig[%d]=%d: ", dim, GridDimension[dim]);
      for (int j = 0; j < GridDimension[dim]; j++)
      printf("%d ", Signature[dim][j]);
      printf("\n");
      } */


  return SUCCESS;
}

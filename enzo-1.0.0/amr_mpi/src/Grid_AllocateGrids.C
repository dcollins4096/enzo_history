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
/  GRID CLASS (ALLOCATE SPACE FOR GRIDS -- DIMS, ETC MUST BE SET)
/
/  written by: Greg Bryan
/  date:       July, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//
//  Allocate room for the grids.
//

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void grid::AllocateGrids()
{

  /* Compute grid size. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Allocate room and clear it. */

  for (int field = 0; field < NumberOfBaryonFields; field++) {
    BaryonField[field]    = new float[size];
    for (int i = 0; i < size; i++)
      BaryonField[field][i] = 0.0;
  }

}

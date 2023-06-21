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
/  GRID CLASS (CLEAR THE MASS FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// Allocate and clear the mass flagging field.

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void grid::ClearMassFlaggingField()
{

  /* error check */

  if (MassFlaggingField != NULL) {
    fprintf(stderr, "ClearMassFlaggingField: Warning, field not deleted.\n");
    delete MassFlaggingField;
  }

  /* compute size and allocate */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  MassFlaggingField = new float[size];

  /* Clear it */

  for (int i = 0; i < size; i++)
    *(MassFlaggingField + i) = 0.0;

}

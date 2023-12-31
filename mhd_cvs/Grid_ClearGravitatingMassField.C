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
/  GRID CLASS (ALLOCATE AND CLEAR THE GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       May, 1995
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

/* function prototypes */


int grid::ClearGravitatingMassField()
{

  /* Return is this is not the right processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Error check. */

  if (GravitatingMassFieldCellSize == FLOAT_UNDEFINED) {
    fprintf(stderr, "GravitatingMassField uninitialized.\n");
    return FAIL;
  }

  /* compute size of the gravitating mass field */

  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldDimension[dim];

  /* allocate and clear the field */

  //  if (GravitatingMassField != NULL)
  //    fprintf(stderr, "ClearGravitatingMassField: Warning! Field not NULL.\n");

  if (GravitatingMassField == NULL)
    GravitatingMassField = new float[size];
  if (GravitatingMassField == NULL) {
    fprintf(stderr, "malloc error (out of memory?)\n");
    return FAIL;
  }

  for (i = 0; i < size; i++)
    GravitatingMassField[i] = 0.0;

  return SUCCESS;
}

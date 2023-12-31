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
/  GRID CLASS (ADD FIELD MASS TO MASS FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
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
/* function prototypes */

int FindField(int f, int farray[], int n);


int grid::AddFieldMassToMassFlaggingField()
{

  int i, DensField;
  float CellVolume;

  if (NumberOfBaryonFields > 0) {

    /* Error check */
    
    if (MassFlaggingField == NULL) {
      fprintf(stderr, "MassFlaggingField not present.\n");
      return FAIL;
    }

    /* Find density field */

    if ((DensField = FindField(Density, FieldType, 
                               NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find density.\n");
      return FAIL;
    }
    
    /* calculate grid size */
    
    int dim, size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

    /* Compute cell volume */

    CellVolume = 1.0;
    for (dim = 0; dim < GridRank; dim++)
      CellVolume *= CellWidth[dim][0];
    
    /* Add mass (method 2). */

    for (i = 0; i < size; i++)
      MassFlaggingField[i] += BaryonField[DensField][i]*CellVolume;

  }

  return SUCCESS;
  
}

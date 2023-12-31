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
/  GRID CLASS (ZERO THE SOLUTION VECTOR UNDER THE SUBGRID SPECIFIED)
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:
/
/  PURPOSE:
/
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



int grid::ZeroSolutionUnderSubgrid(grid *Subgrid, int FieldsToZero,
				   float Value)
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* declarations */
    
  int i, j, k, dim, field, index;
  int SubgridStart[MAX_DIMENSION], SubgridEnd[MAX_DIMENSION];

  /* If FieldsToZero is ZERO_UNDER_SUBGRID_FIELD, and Subgrid is NULL,
     initialize the under subgrid fields (which is the next unused 
     BaryonField even though this pretty much a dumb thing to do). */
  
  if (FieldsToZero == ZERO_UNDER_SUBGRID_FIELD && Subgrid == NULL) {

//    printf("ZeroSUS - ZERO_UNDER_SUBGRID_FIELD && Subgrid == NULL\n");

    int size = 1;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      size *= GridDimension[dim];


    delete [] BaryonField[NumberOfBaryonFields];

    if ( (Unigrid == 1) && (ProblemType == 30) ) 
    {
       printf("ZeroSUS - ZERO_UNDER_SUBGRID_FIELD && Subgrid == NULL\n");
       printf("ZeroSUS - Unigrid: %d\n", Unigrid);
       printf("ZeroSUS - ProblemType: %d\n", ProblemType);
       printf("ZeroSUS - NumberOfBaryonFields: %d\n", NumberOfBaryonFields);
       printf("ZeroSUS - Zero field size: %d\n", (int) (size*sizeof(float)));
       printf("ZeroSUS - CAUTION - allocating bogus BaryonField size = 1\n");
       size = 1;
    }

//    else
//    {
//       printf("ZeroSUS - NumberOfBaryonFields: %d\n", NumberOfBaryonFields);
//       printf("ZeroSUS - Zero field size: %d\n", (int) (size*sizeof(float)));
//    }

    BaryonField[NumberOfBaryonFields] = new float[size];
    for (i = 0; i < size; i++)
      BaryonField[NumberOfBaryonFields][i] = 0.0;
    return SUCCESS;
  }


//  printf("ZeroSUS - Subgrid not NULL\n");

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    SubgridStart[dim] = 0;
    SubgridEnd[dim] = 0;
  }

  /* Compute start and stop indices of the active region of the subgrid
     within this grid (and check to make sure subgrid is within this grid). */

  for (dim = 0; dim < GridRank; dim++) {

    if (Subgrid->GridRightEdge[dim] <= GridLeftEdge[dim] ||
	Subgrid->GridLeftEdge[dim]  >= GridRightEdge[dim])
      return SUCCESS;

    SubgridStart[dim] = nint(
        (Subgrid->GridLeftEdge[dim] - GridLeftEdge[dim])/CellWidth[dim][0]
			       ) + GridStartIndex[dim];
    SubgridEnd[dim] = nint(
	(Subgrid->GridRightEdge[dim] - GridLeftEdge[dim])/CellWidth[dim][0]
			       ) + GridStartIndex[dim] - 1;

    SubgridStart[dim] = max(SubgridStart[dim], GridStartIndex[dim]);
    SubgridEnd[dim]   = min(SubgridEnd[dim], GridEndIndex[dim]);

//    printf("  ZeroSUS: %d, %d, %d\n", dim, SubgridStart[dim], SubgridEnd[dim]);
  }

  /* Now that there is overlap, take the appropriate action. */

  if (FieldsToZero == ZERO_ALL_FIELDS){

//    printf("  ZeroSUS == ZERO_ALL_FIELDS\n");

    /* Zero all fields in this region. */

    for (field = 0; field < NumberOfBaryonFields; field++)
      for (k = SubgridStart[2]; k <= SubgridEnd[2]; k++)
	for (j = SubgridStart[1]; j <= SubgridEnd[1]; j++) {
	  index = (k*GridDimension[1] + j)*GridDimension[0];
	  for (i = SubgridStart[0]; i <= SubgridEnd[0]; i++)
	    BaryonField[field][index + i] = 0;
	}
  }
  else if (FieldsToZero == ZERO_UNDER_SUBGRID_FIELD) {

//    printf("  ZeroSUS = ZERO_UNDER_SUBGRID_FIELD\n");

    if (BaryonField[NumberOfBaryonFields] == NULL) {
      fprintf(stderr, "UNDER_SUBGRID_FIELD not allocated.\n");
      return FAIL;
    }
    
    /* Set points under this subgrid to Value (typically 1). */

//    printf("    ZeroSUS Value = %10.4e\n", Value);

    if ( (Unigrid == 1) && (ProblemType == 30) )
    {
      printf("    Ignoring BaryonField Assignment!\n");
    }
    else
    {
    for (k = SubgridStart[2]; k <= SubgridEnd[2]; k++)
      for (j = SubgridStart[1]; j <= SubgridEnd[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0];
	for (i = SubgridStart[0]; i <= SubgridEnd[0]; i++)
	  BaryonField[NumberOfBaryonFields][index + i] = Value;
      }
    }

  }

  else {
    fprintf(stderr, "FieldsToZero = %d not recognized.\n", FieldsToZero);
    return FAIL;
  }

  return SUCCESS;
  
}

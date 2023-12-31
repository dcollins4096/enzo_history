/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY SLOPE)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::FlagCellsToBeRefinedBySlope()
{
  /* declarations */

  int i, j, k, index, dim;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];

  /* error check */

  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  /* Make sure quantities are defined at least to dim 3 */

  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }

  /* loop over all zones */

  for (dim = 0; dim < 3; dim++) {
    Start[dim] = GridStartIndex[dim];
    End[dim]   = GridEndIndex[dim];
  }

  /* compute size */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* allocate a temporary slope field. */

  float *TempBuffer = new float[size];

  /* loop over active dimensions */

  int Offset = 1;
  for (dim = 0; dim < GridRank; dim++)
    if (GridDimension[dim] > 1) {

      for (int field = 0; field < NumberOfBaryonFields; field++) {

	/* zero slope */

	for (i = 0; i < size; i++)
	  *(TempBuffer + i) = 0.0;

	/* compute slope */

	for (k = Start[2]; k <= End[2]; k++)
	  for (j = Start[1]; j <= End[1]; j++)
	    for (i = Start[0]; i <= End[0]; i++) {
	      index = i + j*GridDimension[0] + 
		      k*GridDimension[1]*GridDimension[0];
	      *(TempBuffer + index) = 0.5*fabs(
		     (*(BaryonField[field] + index + Offset) -
		      *(BaryonField[field] + index - Offset)  ) /
		  max(*(BaryonField[field] + index            ), tiny_number));
	    }
	
	/* flag field based on slope */

	for (i = 0; i < size; i++)
	  *(FlaggingField + i) +=
	    (*(TempBuffer + i) > MinimumSlopeForRefinement) ? 1 : 0;

      }  // end loop over field

      Offset *= GridDimension[dim];

    }  // end loop over dimension

  /* delete buffer */

  delete TempBuffer;

  /* Count number of flagged Cells. */

  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    *(FlaggingField + i) = min(*(FlaggingField + i), 1);
    NumberOfFlaggedCells += *(FlaggingField + i);
  }

  return NumberOfFlaggedCells;

}

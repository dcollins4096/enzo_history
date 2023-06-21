/***********************************************************************
/
/  GRID CLASS (FIND MAXIMUM BARYON DENSITY)
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */



int grid::FindMaximumBaryonDensity(float Position[MAX_DIMENSION],
				   float *MaxDensity)
{
  
  int i, j, k, index;

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);

  /* Loop over grid. */

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
	if (BaryonField[DensNum][index] > *MaxDensity) {
	  *MaxDensity = BaryonField[DensNum][index];
	  Position[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  if (GridRank > 0)
	    Position[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  if (GridRank > 1)
	    Position[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	}
    }

  return SUCCESS;
}

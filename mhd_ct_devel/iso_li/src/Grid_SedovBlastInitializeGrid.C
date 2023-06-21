/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR SEDOV BLAST WAVE TEST) 
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, January 2005.
/
/  PURPOSE: Sets the total energy in the initial explosion region.
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

int grid::SedovBlastInitializeGrid(float dr,
				   float SedovBlastInnerTotalEnergy,
				   grid *TopGrid)
{
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* set fields in the initial explosion region: x^2+y^2 < dr^2. */

  int index, i,j,k;
  float zonex, zoney, zonez, dr2 = dr*dr, center[MAX_DIMENSION];

  for(dim=0;dim<3;dim++){
    center[dim]=0.5*(TopGrid->GridRightEdge[dim]-TopGrid->GridLeftEdge[dim]);

  }

  for (k = 0; k < GridDimension[2]; k++) 
    for (j = 0; j < GridDimension[1]; j++) 
      for (i = 0; i < GridDimension[0]; i++) {

	zonex = *(CellLeftEdge[0] + i) + 0.5*(*(CellWidth[0] + i)) - center[0];
	zoney = *(CellLeftEdge[1] + j) + 0.5*(*(CellWidth[1] + j)) - center[1];
	if(GridRank >= 3)
	  zonez=*(CellLeftEdge[2] + k) + 0.5*(*(CellWidth[2] + k)) - center[2];
	  else
	  zonez=0;

	index = i + GridDimension[0]*(j+GridDimension[1]*k);

	if (zonex*zonex + zoney*zoney + zonez*zonez< dr2)
	  BaryonField[1][index] = SedovBlastInnerTotalEnergy;

      }

  return SUCCESS;
}

/***********************************************************************
/
/  GRID CLASS (RETURN FLUX DIMS IN FLUX STRUCTURE)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// Return flux dimensions in flux structure in argument, divided
//   by the refinement factors.
//   It also sets all the flux field pointers to NULL.

#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::ReturnFluxDims(fluxes &Flux, int RefinementFactors[])
{
  for (int i = 0; i < GridRank; i++)
    for (int j = 0; j < GridRank; j++) {
      Flux.LeftFluxStartGlobalIndex[i][j] = 
	BoundaryFluxes.LeftFluxStartGlobalIndex[i][j]/RefinementFactors[j];
      Flux.LeftFluxEndGlobalIndex[i][j] = 
	BoundaryFluxes.LeftFluxEndGlobalIndex[i][j]/RefinementFactors[j];
      Flux.RightFluxStartGlobalIndex[i][j] = 
	BoundaryFluxes.RightFluxStartGlobalIndex[i][j]/RefinementFactors[j];
      Flux.RightFluxEndGlobalIndex[i][j] = 
	BoundaryFluxes.RightFluxEndGlobalIndex[i][j]/RefinementFactors[j];
    }

  /* set all to NULL for easier clean-up. */

  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    for (int field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
      Flux.LeftFluxes[field][dim] = NULL;
      Flux.RightFluxes[field][dim] = NULL;
    }

  return SUCCESS;

}

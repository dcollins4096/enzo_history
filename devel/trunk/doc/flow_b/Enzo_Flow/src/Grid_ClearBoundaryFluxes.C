/***********************************************************************
/
/  GRID CLASS (CLEAR/ALLOCATE BOUNDARY FLUXES FOR THIS GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// Clear the boundary fluxes (allocate first if necesary).
//  Note that we assume the left and right flux boundaries are the same size.
//   

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void grid::ClearBoundaryFluxes()
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return;

  int dim, i, field, size;

  for (dim = 0; dim < GridRank; dim++) {

    /* compute size for this dim's face */

    size = 1;
    for (i = 0; i < GridRank; i++)
      size *= BoundaryFluxes.LeftFluxEndGlobalIndex[dim][i] -
	      BoundaryFluxes.LeftFluxStartGlobalIndex[dim][i] + 1;

    /* allocate if necesary, and clear */

    for (field = 0; field < NumberOfBaryonFields; field++) {
      if (BoundaryFluxes.LeftFluxes[field][dim] == NULL) {
	BoundaryFluxes.LeftFluxes[field][dim] = new float[size];
	BoundaryFluxes.RightFluxes[field][dim] = new float[size];
      }
      for (i = 0; i < size; i++) {
	BoundaryFluxes.LeftFluxes[field][dim][i] = 0.0;
	BoundaryFluxes.RightFluxes[field][dim][i] = 0.0;
      }
    }

    /* set unused pointers to NULL */

    for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	 field++)
      for (i = 0; i < MAX_DIMENSION; i++) {
	BoundaryFluxes.LeftFluxes[field][i]  = NULL;
	BoundaryFluxes.RightFluxes[field][i] = NULL;
      }
    
  } // end loop over dims

  /* set unused pointers to NULL */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++)
    for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
      BoundaryFluxes.LeftFluxes[field][dim]   = NULL;
      BoundaryFluxes.RightFluxes[field][dim] = NULL;
    }

}

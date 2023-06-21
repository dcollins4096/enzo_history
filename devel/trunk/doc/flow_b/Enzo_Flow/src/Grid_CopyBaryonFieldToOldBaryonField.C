/***********************************************************************
/
/  GRID CLASS (COPY BARYON FIELD TO OLD BARYON FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// Copy the current baryon fields to the old baryon fields
//   (allocate old baryon fields if they don't exist).

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::CopyBaryonFieldToOldBaryonField()
{
    
  /* update the old baryon field time */

  OldTime = Time;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* compute the field size */

  int size = 1;
  for (int dim = 0; dim < GridDimension[dim]; dim++)
    size *= GridDimension[dim];

  /* copy fields */
    
  for (int field = 0; field < NumberOfBaryonFields; field++) {

    /* Check to make sure BaryonField exists. */

    if (BaryonField[field] == NULL) {
      fprintf(stderr, "BaryonField missing.\n");
      return FAIL;
    }

    /* Create OldBaryonField if necessary. */

    if ((OldBaryonField[field] == NULL))
      OldBaryonField[field] = new float[size];

    /* Copy. */

    for (int i = 0; i < size; i++)
      OldBaryonField[field][i] = BaryonField[field][i];

  } // end loop over fields

  return SUCCESS;

}

/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (DETACH RANDOM FORCING FIELDS FROM THE LIST OF
/                           BARYON FIELDS)
/
/  written by: Alexei Kritsuk
/  date:       January 2004
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
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
 
int ExternalBoundary::DetachForcingFromBaryonFields()
{
 
  /* BaryonFields. */
  printf("ForcingDetachedFromBoundary, Before: NBF %"ISYM"\n", NumberOfBaryonFields); 
#ifndef PPML
  for (int dim = 0; dim < BoundaryRank; dim++)
    BoundaryFieldType[NumberOfBaryonFields-1-dim] = FieldUndefined;
  NumberOfBaryonFields -= BoundaryRank;
#else //PPML
  for( int dim = 0; dim< RandomForcingNumberOfFields; dim++){
    BoundaryFieldType[NumberOfBaryonFields-1-dim] = FieldUndefined;
  }
  NumberOfBaryonFields -= RandomForcingNumberOfFields;
#endif //PPML

    printf("ForcingDetachedFromBoundary NBF %"ISYM"\n", NumberOfBaryonFields);
 
  return SUCCESS;
 
}

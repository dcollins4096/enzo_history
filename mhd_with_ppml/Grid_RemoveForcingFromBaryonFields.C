/***********************************************************************
/
/  GRID CLASS (REMOVE RANDOM FORCING FIELDS FROM THE LIST OF BARYON FIELDS)
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
 
int grid::RemoveForcingFromBaryonFields()
{
 
  /* check.
 
  if (debug)
    printf("ForcingToBeRemoved[%"ISYM"] NBF %"ISYM"\n", ProcessorNumber,
    NumberOfBaryonFields - GridRank);*/
  if (ProcessorNumber == MyProcessorNumber && debug)
    if (BaryonField[NumberOfBaryonFields-1] == NULL)
      WARNING_MESSAGE;
 
  /* Remove RandomForcingFields from BaryonFields. */

#ifdef PPML
  for (int dim = 0; dim < RandomForcingNumberOfFields; dim++) {
    RandomForcingField[RandomForcingNumberOfFields-dim-1] =
      BaryonField[NumberOfBaryonFields-dim-1];
    BaryonField[NumberOfBaryonFields-dim-1] = NULL;
    FieldType[NumberOfBaryonFields-dim-1] = FieldUndefined;
   }
 
  NumberOfBaryonFields -= RandomForcingNumberOfFields;
#else //PPML
  for (int dim = 0; dim < GridRank; dim++) {
    RandomForcingField[GridRank-dim-1] =
      BaryonField[NumberOfBaryonFields-dim-1];
    BaryonField[NumberOfBaryonFields-dim-1] = NULL;
    FieldType[NumberOfBaryonFields-dim-1] = FieldUndefined;
   }
 
  NumberOfBaryonFields -= GridRank;


#endif //PPML
 
  return SUCCESS;
 
}

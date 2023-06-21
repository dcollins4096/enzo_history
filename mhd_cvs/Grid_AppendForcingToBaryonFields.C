/***********************************************************************
/
/  GRID CLASS (ADD RANDOM FORCING FIELDS TO THE LIST OF BARYON FIELDS)
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

int grid::AppendForcingToBaryonFields()
{

  /* check if this does concern us. */
  
  if (ProcessorNumber == MyProcessorNumber) {
    int size = 1;
    for (int dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
    if (RandomForcingField[0] == NULL) {
      RandomForcingField[0] = new float[size];
      for (int i = 0; i < size; i++)
	RandomForcingField[0][i] = 0.0;
    }
    if (GridRank > 1)
      if (RandomForcingField[1] == NULL) {
	RandomForcingField[1] = new float[size];
	for (int i = 0; i < size; i++)
	  RandomForcingField[1][i] = 0.0;
      }
    if (GridRank > 2)
      if (RandomForcingField[2] == NULL) {
	RandomForcingField[2] = new float[size];
	for (int i = 0; i < size; i++)
	  RandomForcingField[2][i] = 0.0;
      }
  }

  /* Add RandomForcingFields as extra BaryonFields. */
  char *RandomName_0 = "RandomForcing_0";
  char *RandomName_1 = "RandomForcing_1";
  char *RandomName_2 = "RandomForcing_2";
  BaryonField[NumberOfBaryonFields] = RandomForcingField[0];
  DataLabel[NumberOfBaryonFields] = RandomName_0;
  DataUnits[NumberOfBaryonFields] = NULL;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) {
    BaryonField[NumberOfBaryonFields] = RandomForcingField[1];
    DataLabel[NumberOfBaryonFields] = RandomName_1;
    DataUnits[NumberOfBaryonFields] = NULL;
    FieldType[NumberOfBaryonFields++] = Velocity2;
  }
  if (GridRank > 2) {
    BaryonField[NumberOfBaryonFields] = RandomForcingField[2];
    DataLabel[NumberOfBaryonFields] = RandomName_2;
    DataUnits[NumberOfBaryonFields] = NULL;
    FieldType[NumberOfBaryonFields++] = Velocity3;
  }

  /*
  printf("ForcingAppended[Proc %d, MyProc %d] NBF %d\n", ProcessorNumber, MyProcessorNumber,
	 NumberOfBaryonFields);
  */
  return SUCCESS;

}

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
#ifdef PPML
    if( RandomForcingNumberOfFields > GridRank ){
      for(int field=GridRank; field<RandomForcingNumberOfFields; field++){
	if( RandomForcingField[field] == NULL ){
	  RandomForcingField[field] = new float[size];
	  for(int i = 0; i<size; i++)
	    RandomForcingField[field][i] = 0.0;
	}
      }
    }
    //error check that there are enough possible baryon fields
    if( NumberOfBaryonFields + RandomForcingNumberOfFields > MAX_NUMBER_OF_BARYON_FIELDS ){
      fprintf(stderr, "Fatal error: AppendForcingToBaryonFields: MAX_NUMBER_OF_BARYON_FIELDS not large enough.\n");
      fprintf(stderr, "Must be at least as large as NumberOfBaryonFields (%"ISYM") + RandomForcingNumberOfFields (%d)+1\n",
	      NumberOfBaryonFields, RandomForcingNumberOfFields);
      ERROR_MESSAGE;
    }
#endif //PPML    
  } 
  /* Add RandomForcingFields as extra BaryonFields. */
 
  BaryonField[NumberOfBaryonFields] = RandomForcingField[0];
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) {
    BaryonField[NumberOfBaryonFields] = RandomForcingField[1];
    FieldType[NumberOfBaryonFields++] = Velocity2;
  }
  if (GridRank > 2) {
    BaryonField[NumberOfBaryonFields] = RandomForcingField[2];
    FieldType[NumberOfBaryonFields++] = Velocity3;
  }
#ifdef PPML
  if( RandomForcingNumberOfFields > GridRank ){
    int counter=GridRank;
    for( int field=0; field<PPML_NFaces; field++){
      BaryonField[NumberOfBaryonFields] = RandomForcingField[counter++];
      FieldType[NumberOfBaryonFields++] = Velocity1;
    }
    for( int field=0; field<PPML_NFaces; field++){
      BaryonField[NumberOfBaryonFields] = RandomForcingField[counter++];
      FieldType[NumberOfBaryonFields++] = Velocity2;
    }
    for( int field=0; field<PPML_NFaces; field++){
      BaryonField[NumberOfBaryonFields] = RandomForcingField[counter++];
      FieldType[NumberOfBaryonFields++] = Velocity3;
    }
  }
#endif //PPML 
  /*  if (debug)
    printf("ForcingAppended[%"ISYM"] NBF %"ISYM"\n", ProcessorNumber,
	   NumberOfBaryonFields);
  */
 
  return SUCCESS;
 
}

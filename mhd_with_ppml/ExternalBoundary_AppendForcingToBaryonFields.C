/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (ADD RANDOM FORCING FIELDS TO THE LIST OF
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
 
int ExternalBoundary::AppendForcingToBaryonFields()
{
 
  /* BaryonFields. */
 
  BoundaryFieldType[NumberOfBaryonFields++] = Velocity1;
  if (BoundaryRank > 1)
    BoundaryFieldType[NumberOfBaryonFields++] = Velocity2;
  if (BoundaryRank > 2)
  BoundaryFieldType[NumberOfBaryonFields++] = Velocity3;
#ifdef PPML
  if( RandomForcingNumberOfFields > BoundaryRank ){
    int counter=0;
    do{
      BoundaryFieldType[NumberOfBaryonFields++] = Velocity1;
      counter++;
      if (BoundaryRank > 1){
	BoundaryFieldType[NumberOfBaryonFields++] = Velocity2;
	counter++;
      }
      if (BoundaryRank > 2){
	BoundaryFieldType[NumberOfBaryonFields++] = Velocity3;      
	counter++;
      }
    }while(counter < RandomForcingNumberOfFields - BoundaryRank ); //face
  }//more random forcing.
#endif //PPM 

    printf("ForcingAppendedToBoundary NBF %"ISYM"\n", NumberOfBaryonFields);
 
  return SUCCESS;
 
}

/***********************************************************************
/
/  GRID CLASS (ADD RANDOM FORCING FIELDS TO VELOCITIES)
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
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

#ifdef PPML
#include "PPML.h"
#endif //PPML 

int grid::AddRandomForcing(float * norm, float dtTopGrid)
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "GARF: Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
                                                                             
  /* error check. */
  
  if (RandomForcingField[0] == NULL) 
    ERROR_MESSAGE;
  
  if (RandomForcingField[0][0] == 0.0) 
    ERROR_MESSAGE;

  if (dtTopGrid == 0.0)
    ERROR_MESSAGE;


  /* compute the field size */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* update total energy first. */


  float levelNorm = (*norm)*dtFixed/dtTopGrid;

  if (levelNorm <= 0.0)
    WARNING_MESSAGE;

#ifdef PPML
  if( EquationOfState == 0 && RandomForcingNumberOfFields > 3 ){
    fprintf(stderr,"Fatal Error: Face centered forcing not installed in the energy field.\n");
  }
#endif //PPML
  /* do not do the update if using ZEUS */
  if (HydroMethod != Zeus_Hydro )
#ifdef PPML
    if(  EquationOfState == 0 )
#endif //PPML

#ifdef PPML
    if(MHD_Used != TRUE ){
      for (int i = 0; i < size; i++)
	for (int dim = 0; dim < GridRank; dim++)
	  BaryonField[TENum][i] +=
	    BaryonField[Vel1Num+dim][i]*
	    RandomForcingField[dim][i]*levelNorm +
	    0.5*RandomForcingField[dim][i]*levelNorm*
	    RandomForcingField[dim][i]*levelNorm;
    }else{
      //hey-- if you don't tab this over, it will match what's in the trunk.
#endif //PPML
    for (int i = 0; i < size; i++)
      for (int dim = 0; dim < GridRank; dim++)
	BaryonField[TENum][i] +=
	  BaryonField[Vel1Num+dim][i]*
	  RandomForcingField[dim][i]*levelNorm +
	  0.5*RandomForcingField[dim][i]*levelNorm*
	  RandomForcingField[dim][i]*levelNorm;

#ifdef PPML
    }

#endif //PPML
  /* add velocity perturbation to velocity fields. */

  for (int dim = 0; dim < GridRank; dim++)
    for (int i = 0; i < size; i++)
      BaryonField[Vel1Num+dim][i] += RandomForcingField[dim][i]*levelNorm;

#ifdef PPML

  if( HydroMethod == PPM_Local && RandomForcingNumberOfFields > GridRank ){

    PPML_InterfacePointerBundle Face( this ); //filled after InitInterfaceTypesAndLabels
    if( Face.Error == TRUE ){
      fprintf( stderr,"Interface Bundle Error.\n");
      ERROR_MESSAGE;
    }
#ifdef trash    
    if( this->ReturnInterfacePointers( Face ) == FAIL ){
      fprintf(stderr," ReturnInterfacePointers failed.\n"); return FAIL;}
#endif //trash
    int ForcingCounter = GridRank; //picking up where we left off.

    for( int face_part = 0; face_part < PPML_NFaces; face_part++)
      for (int dim = 0; dim < GridRank; dim++){
	for (int i = 0; i < size; i++)
	  Face.All[face_part][Vel1Num+dim][i] += RandomForcingField[ForcingCounter][i]*levelNorm;
	ForcingCounter++;
      }

#ifdef trash
    if( this->DeleteInterfacePointers( Face ) == FAIL ){
      fprintf(stderr,"Grid_PPML_Mono1: error with Delete InterfacePointers");
      return FAIL;
    }
#endif //trash    
  }//extra force
#endif //PPML
  return SUCCESS;
  
}
  

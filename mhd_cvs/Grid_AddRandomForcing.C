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
 
int grid::AddRandomForcing(float * norm, float * bulkMomentum, float dtTopGrid)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum,i,dim;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "GARF: Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
 
  /* error check. */
 
  if (RandomForcingField[0] == NULL) {
    fprintf(stderr,"random forcing field NULL\n");
    ERROR_MESSAGE;
  }
  
  if (RandomForcingField[0][0] == 0.0) {
    fprintf(stderr,"Error with Random forcing field Field[0][0]==0\n");
    ERROR_MESSAGE;
  }
  

  if (dtTopGrid == 0.0) {
    fprintf(stderr,"Random Forcing: dtTopGrid = 0\n");
    ERROR_MESSAGE;
  }
 
  /* compute the field size */
 
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* update total energy first. */
 
 
  float levelNorm = (*norm)*dtFixed/dtTopGrid;
  if (debug & levelNorm < 0.0)
    WARNING_MESSAGE;
  /* Alexei's no-dt^2 fix 
     // do not do the update if using ZEUS 
  float useDensity = 1;
  if (HydroMethod != Zeus_Hydro && EquationOfState == 0)
    for ( i = 0; i < size; i++)
      for (dim = 0; dim < GridRank; dim++){
	if( MHD_Used ) useDensity = BaryonField[DensNum][i];
	BaryonField[TENum][i] +=
	  useDensity*(BaryonField[Vel1Num+dim][i]*
	   //	  RandomForcingField[dim][i])*levelNorm;
  	  (RandomForcingField[dim][i]-bulkMomentum[dim]))*levelNorm;
      }
  */

/* Or using Isothermal EOS. */
  if (HydroMethod != Zeus_Hydro && EquationOfState == 0 )
    if(MHD_Used != TRUE ){
      for (i = 0; i < size; i++)
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[TENum][i] +=
	    BaryonField[Vel1Num+dim][i]*
	    RandomForcingField[dim][i]*levelNorm +
	    0.5*RandomForcingField[dim][i]*levelNorm*
	    RandomForcingField[dim][i]*levelNorm;
    }else{
      for (i = 0; i < size; i++)
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[TENum][i] +=
	    BaryonField[DensNum][i]*(
				     BaryonField[Vel1Num+dim][i]*
				     RandomForcingField[dim][i]*levelNorm +
				     0.5*RandomForcingField[dim][i]*levelNorm*
				     RandomForcingField[dim][i]*levelNorm);
    }
      
  /* add velocity perturbation to the velocity fields;
     keep the center-of-mass velocity zero. */
 
  if (debug)
    fprintf(stderr, "GARF: Bulk Momentum input %g %g %g %g\n", 
	    bulkMomentum[0], bulkMomentum[1], bulkMomentum[2], levelNorm);
      

  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      //      BaryonField[Vel1Num+dim][i] += RandomForcingField[dim][i]*levelNorm;
      BaryonField[Vel1Num+dim][i] += (RandomForcingField[dim][i]-bulkMomentum[dim])*levelNorm;
 
  return SUCCESS;
 
}


#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#include "GlobalStats.h"
#ifdef PPML
//Only for the index pointer map.
#include "PPML.h" 
#else
#include "IndexPointerMap.h"
#endif
int grid::ComputeGlobalStatsGrid(GlobalStats * Stats){

  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;

  IndexPointerMap ind;
  if( this->IdentifyPhysicalQuantities_2(ind) == FAIL ){
    fprintf(stderr," IdentifyPhysicalQuantities_2 failed\n");return FAIL;}



  /* Compute the pressure. */

  int size = GridDimension[0]*GridDimension[1]* GridDimension[2];
  int result;
  float *pressure_field = new float[size];
  if (DualEnergyFormalism)
    result = this->ComputePressureDualEnergyFormalism(Time, pressure_field);
  else
    result = this->ComputePressure(Time, pressure_field);
  
    if (result == FAIL) {
      fprintf(stderr, "Error in grid->ComputeGlobalStats, bad pressure field.\n");
      return FAIL;
    }

  //Temp Variables for computation
  float SoundSpeedSquared, Vsquared, Bsquared;
  Stat_Obj * element;
  //loop over active zones
  int i,j,k, index;
  for( k=GridStartIndex[2]; k<= GridEndIndex[2]; k++)
    for( j=GridStartIndex[1]; j<= GridEndIndex[1]; j++)
      for( i=GridStartIndex[0]; i<= GridEndIndex[0]; i++){
	index = i + GridDimension[0]*(j + GridDimension[1]*k);

	
	Vsquared =  (BaryonField[ ind.VX ][ index] * BaryonField[ ind.VX ][ index] +
		     ( (GridRank < 2 && MHD_Used != TRUE ) ? 0 : 
		       BaryonField[ ind.VY ][ index] * BaryonField[ ind.VY ][ index] ) +
		     ( (GridRank < 3 && MHD_Used != TRUE ) ? 0 :
		       BaryonField[ ind.VZ ][ index] * BaryonField[ ind.VZ ][ index] ) );

	if( element = Stats->elem("V-rms") )
	  element->value += Vsquared;

	if( element = Stats->elem("MaxV^2") )
	  element->value = max( element->value, Vsquared);

	if( element = Stats->elem("AvgVx") )
	  element->value += BaryonField[ ind.VX ][ index];
	if( element = Stats->elem("AvgVy") )
	  element->value += BaryonField[ ind.VY ][ index];
	if( element = Stats->elem("AvgVz") )
	  element->value += BaryonField[ ind.VZ ][ index];
	
	if( MHD_Used == TRUE ){
#ifdef PPML
	  Bsquared = (BaryonField[ ind.BX ][ index] * BaryonField[ ind.BX ][ index] +
		      BaryonField[ ind.BY ][ index] * BaryonField[ ind.BY ][ index] +
		      BaryonField[ ind.BZ ][ index] * BaryonField[ ind.BZ ][ index] );
	  
	  if( element = Stats->elem("B-rms") )
	    element->value += Bsquared;
	  if( element = Stats->elem("AvgBx") )
	    element->value += BaryonField[ind.BX][index];
	  if( element = Stats->elem("AvgBy") )
	    element->value += BaryonField[ind.BY][index];
	  if( element = Stats->elem("AvgBz") )
	    element->value += BaryonField[ind.BZ][index];
#else
	  Bsquared = (CenteredB[ 0 ][index]*CenteredB[ 0 ][index]+
		      CenteredB[ 1 ][index]*CenteredB[ 1 ][index]+
		      CenteredB[ 2 ][index]*CenteredB[ 2 ][index]);

	  if( element = Stats->elem("B-rms") )
	    element->value += Bsquared;
	  if( element = Stats->elem("AvgBx") )
	    element->value += CenteredB[0][index];
	  if( element = Stats->elem("AvgBy") )
	    element->value += CenteredB[1][index];
	  if( element = Stats->elem("AvgBz") )
	    element->value += CenteredB[2][index];
	  
#endif
	}else{
	  Bsquared = 0;
	}
	if( EquationOfState == 0 ){
	  SoundSpeedSquared = Gamma*pressure_field[index]/BaryonField[ ind.D ][index];
	}else{
	  SoundSpeedSquared = IsothermalSoundSpeed*IsothermalSoundSpeed;
	}

	if( element = Stats->elem("MaxB^2") )
	  element->value= max(element->value,Bsquared);

	if( element=Stats->elem("RMS-mach") )
	  element->value += Vsquared / SoundSpeedSquared;

	  //Variance assumes avg density = 1.
	if( element=Stats->elem("DensityVariance") )
	  element->value += (BaryonField[ ind.D ][ index ] - 1.0 )*
	    (BaryonField[ ind.D ][ index ] - 1.0 );

	if( element=Stats->elem("AvgKE") )
	  element->value += 0.5*BaryonField[ ind.D ][index]* Vsquared;

	if( element=Stats->elem("MinD") )
	  element->value = min(element->value, BaryonField[ ind.D ][index] );
	if( element=Stats->elem("MaxD") )
	  element->value = max(element->value, BaryonField[ ind.D ][index] );


	if( fabs( Bsquared ) > tiny_number ){
	  if( element = Stats->elem("RMS-AlvMach") )
	    element->value += Vsquared*BaryonField[ ind.D][index]/Bsquared;
	  if( element = Stats->elem("AvgBeta") )
	    element->value += 2*pressure_field[index]/Bsquared;
	}else{
	  if( element = Stats->elem("RMS-AlvMach") )
	    element->value += 0;
	  if( element = Stats->elem("AvgBeta") )
	    element->value += 0;
	}

	  
	  

	  
	  

      }

  delete pressure_field;
  return SUCCESS;
}

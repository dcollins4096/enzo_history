
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

  float dx = 1, dy = 1, dz = 1;
  float KE;
  if( 0 == 0 ){
    dx = CellWidth[0][0];
    dy = CellWidth[1][0];
    dz = CellWidth[2][0];
  }


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

	if( element = Stats->elem(V_rms) )
	  element->value += Vsquared;

	if( element = Stats->elem(MaxV2) )
	  element->value = max( element->value, Vsquared);

	if( element = Stats->elem(AvgVx) )
	  element->value += BaryonField[ ind.VX ][ index];
	if( element = Stats->elem(AvgVy) )
	  element->value += BaryonField[ ind.VY ][ index];
	if( element = Stats->elem(AvgVz) )
	  element->value += BaryonField[ ind.VZ ][ index];


	if( element = Stats->elem(px) )
	  element->value += BaryonField[ ind.VX ][ index]*BaryonField[ ind.D ][ index ];
	if( element = Stats->elem(py) )
	  element->value += BaryonField[ ind.VY ][ index]*BaryonField[ ind.D ][ index ];
	if( element = Stats->elem(pz) )
	  element->value += BaryonField[ ind.VZ ][ index]*BaryonField[ ind.D ][ index ];
	
	if( MHD_Used == TRUE ){
#ifdef PPML
	  Bsquared = (BaryonField[ ind.BX ][ index] * BaryonField[ ind.BX ][ index] +
		      BaryonField[ ind.BY ][ index] * BaryonField[ ind.BY ][ index] +
		      BaryonField[ ind.BZ ][ index] * BaryonField[ ind.BZ ][ index] );
	  
	  if( element = Stats->elem(B_rms) )
	    element->value += Bsquared;
	  if( element = Stats->elem(AvgBx) )
	    element->value += BaryonField[ind.BX][index];
	  if( element = Stats->elem(AvgBy) )
	    element->value += BaryonField[ind.BY][index];
	  if( element = Stats->elem(AvgBz) )
	    element->value += BaryonField[ind.BZ][index];
#else
	  Bsquared = (CenteredB[ 0 ][index]*CenteredB[ 0 ][index]+
		      CenteredB[ 1 ][index]*CenteredB[ 1 ][index]+
		      CenteredB[ 2 ][index]*CenteredB[ 2 ][index]);

	  if( element = Stats->elem(B_rms) )
	    element->value += Bsquared;
	  if( element = Stats->elem(AvgBx) )
	    element->value += CenteredB[0][index];
	  if( element = Stats->elem(AvgBy) )
	    element->value += CenteredB[1][index];
	  if( element = Stats->elem(AvgBz) )
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

	if( element = Stats->elem(MaxB2) )
	  element->value= max(element->value,Bsquared);

	if( element=Stats->elem(RMS_mach) )
	  element->value += Vsquared / SoundSpeedSquared;

	  //Variance assumes avg density = 1.
	if( element=Stats->elem(DensityVariance) )
	  element->value += (BaryonField[ ind.D ][ index ] - 1.0 )*
	    (BaryonField[ ind.D ][ index ] - 1.0 );

	KE = 0.5*BaryonField[ ind.D ][index]* Vsquared;
	if( element=Stats->elem(AvgKE) )
	  element->value += KE;
	if( element=Stats->elem(AvgTE) )
	  element->value += BaryonField[ind.TE][index];

	if( element=Stats->elem(AvgGE) )
	  element->value += BaryonField[ind.TE][index] - KE - Bsquared*0.5;
	if( element=Stats->elem(AvgBE) )
	  element->value += Bsquared*0.5;



	if( element=Stats->elem(MinD) )
	  element->value = min(element->value, BaryonField[ ind.D ][index] );
	if( element=Stats->elem(MaxD) )
	  element->value = max(element->value, BaryonField[ ind.D ][index] );


	if( fabs( Bsquared ) > tiny_number ){
	  if( element = Stats->elem(RMSAlvMach) )
	    element->value += Vsquared*BaryonField[ ind.D][index]/Bsquared;
	  if( element = Stats->elem(AvgBeta) )
	    element->value += 2*pressure_field[index]/Bsquared;
	  if( element = Stats->elem(AvgAlfv) )
	    element->value += Bsquared/BaryonField[ ind.D][index];

	  if( element = Stats->elem(MaxDivB) )
	    element->value = max( element->value, fabs(
	    (MagneticField[0][indexb1(i+1,j,k)] -MagneticField[0][indexb1(i,j,k)])/ (dx)+
	    ( (GridRank < 2 ) ? 0 : 
	      (MagneticField[1][indexb2(i,j+1,k)] -MagneticField[1][indexb2(i,j,k)])/ (dy )) +
	    ( (GridRank < 3 ) ? 0 : 
	      (MagneticField[2][indexb3(i,j,k+1)] -MagneticField[2][indexb3(i,j,k)])/ (dz )) ) ); 

	}else{
	  if( element = Stats->elem(RMSAlvMach) )
	    element->value += 0;
	  if( element = Stats->elem(AvgBeta) )
	    element->value += 0;
	}

	  
	  
	if( CellLeftEdge[0][i] > 0.5 - 0.1*CellWidth[0][0] 
	    && CellLeftEdge[0][i+1] < 1.0/128.0 + 0.1*CellWidth[0][0])
	if( CellLeftEdge[1][i] > 0.5 - 0.1*CellWidth[1][0] 
	    && CellLeftEdge[1][i+1] < 1.0/128.0 + 0.1*CellWidth[1][0])
	if( CellLeftEdge[2][i] > 0.5 - 0.1*CellWidth[2][0] 
	    && CellLeftEdge[2][i+1] < 1.0/128.0 + 0.1*CellWidth[2][0])
	  {
	    if( element = Stats->elem(C_rho) )
	      element->value += BaryonField[ ind.D ][ index];
	    if( element = Stats->elem(C_ke) )
	      element->value += (0.5*BaryonField[ ind.D ][ index]*
				 (BaryonField[ ind.VX ][ index]*BaryonField[ ind.VX ][ index]+
				  BaryonField[ ind.VY ][ index]*BaryonField[ ind.VY ][ index]+
				  BaryonField[ ind.VZ ][ index]*BaryonField[ ind.VZ ][ index]));
	    if( element = Stats->elem(C_mach) )
	      element->value += ((BaryonField[ ind.VX ][ index]*BaryonField[ ind.VX ][ index]+
				  BaryonField[ ind.VY ][ index]*BaryonField[ ind.VY ][ index]+
				  BaryonField[ ind.VZ ][ index]*BaryonField[ ind.VZ ][ index]));
	      
	    if( element = Stats->elem(C_px) )
	      element->value += BaryonField[ ind.D ][ index]*BaryonField[ ind.VX ][ index];
	    if( element = Stats->elem(C_py) )
	      element->value += BaryonField[ ind.D ][ index]*BaryonField[ ind.VY ][ index];
	    if( element = Stats->elem(C_pz) )
	      element->value += BaryonField[ ind.D ][ index]*BaryonField[ ind.VZ ][ index];
	    if( element = Stats->elem(C_be) )
	      if( element = Stats->elem(C_bx) )
		element->value += CenteredB[0][index];
	    if( element = Stats->elem(C_by) )
		element->value += CenteredB[1][index];
	    if( element = Stats->elem(C_bz) )
		element->value += CenteredB[2][index];
	    if( element = Stats->elem(C_beta) )
		element->value = 2*pressure_field[index]/Bsquared;
	  }

      }

  delete pressure_field;
  return SUCCESS;
}

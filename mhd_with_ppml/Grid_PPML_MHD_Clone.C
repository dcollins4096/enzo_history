
#ifdef PPML
#ifdef MHDF

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "PPML.h"

//Direction = 0: PPML = MagneticField
//Direction = 1: MagneticField = PPML
int grid::PPML_MHD_Clone(int Direction){
  
  int i,j,k, field, indexP, indexB_L, indexB_R;
  PPML_InterfacePointerBundle Face( this );
  IndexPointerMap ind;
  if( this->IdentifyPhysicalQuantities_2(ind) == FAIL ){
    fprintf(stderr," IdentifyPhysicalQuantities_2 failed\n");return FAIL;}

  int Offset[3] = {1, MagneticDims[1][0], MagneticDims[2][0]*MagneticDims[2][1]};
  for( k=0;k<GridDimension[2];k++)
    for( j=0;j<GridDimension[1];j++)
      for( i=0;i<GridDimension[0];i++){
	indexP = i + GridDimension[0] * ( j + GridDimension[1]*k);

	if( 0 ==  Direction ){
	  field = 0;
	  indexB_L = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k);
	  indexB_R = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k) + Offset[field];
	  Face.X_L[ ind.B[field] ][ indexP ] = 
	    MagneticField[ field ][ indexB_L ];
	  Face.X_R[ ind.B[field] ][ indexP ] = 
	    MagneticField[ field ][ indexB_R ];
	  
	  field = 1;
	  indexB_L = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k);
	  indexB_R = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k) + Offset[field];
	  Face.Y_L[ ind.B[field] ][ indexP ] = 
	    MagneticField[ field ][ indexB_L ];
	  Face.Y_R[ ind.B[field] ][ indexP ] = 
	    MagneticField[ field ][ indexB_R ];
	  
	  field = 2;
	  indexB_L = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k);
	  indexB_R = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k) + Offset[field];
	  Face.Z_L[ ind.B[field] ][ indexP ] = 
	    MagneticField[ field ][ indexB_L ];
	  Face.Z_R[ ind.B[field] ][ indexP ] = 
	    MagneticField[ field ][ indexB_R ];
	}else{ //then Direction == 1
	  field = 0;
	  indexB_L = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k);
	  indexB_R = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k) + Offset[field];
	  MagneticField[ field ][ indexB_L ] =
	    Face.X_L[ ind.B[field] ][ indexP ];

	  MagneticField[ field ][ indexB_R ] =
	    Face.X_R[ ind.B[field] ][ indexP ];
	    
	  
	  field = 1;
	  indexB_L = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k);
	  indexB_R = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k) + Offset[field];
	  MagneticField[ field ][ indexB_L ] =
	    Face.Y_L[ ind.B[field] ][ indexP ];
	    
	  MagneticField[ field ][ indexB_R ] =
	    Face.Y_R[ ind.B[field] ][ indexP ];
	  
	  field = 2;
	  indexB_L = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k);
	  indexB_R = i + MagneticDims[field][0]*( j + MagneticDims[field][1]*k) + Offset[field];
	  MagneticField[ field ][ indexB_L ] =
	    Face.Z_L[ ind.B[field] ][ indexP ];
	  MagneticField[ field ][ indexB_R ] =
	    Face.Z_R[ ind.B[field] ][ indexP ];
	    
	}
	
      }//loop
  
  return SUCCESS;
}


#endif //PPML
#endif //MHD

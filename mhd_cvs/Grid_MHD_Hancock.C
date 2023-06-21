

#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "pout.h"
#include "MHD_Athena.h"
#include "MHDFluid.h"
#ifdef ATHENA

//
// Hancock version 1.
//
// U^{n-1/2}_{i-1/2,L}= U^{n}_{i-1} - 0.5*dtdx*( F(U_{i-1,R}) - F(U_{i-1,L} ) + 0.5*dU_{i-1}
//    U_{i-1,R} = U^{n}_{i-1} + 0.5* dU_{i-1}
//    U_{i-1,L} = U^{n}_{i-1} - 0.5* dU_{i-1}
//    dU_{i-1} is some kind of limited slope.
// U^{n-1/2}_{i-1/2,R} = U^{n}_{i} - 0.5*dtdx*( F(U_{i,R} - F(U_{i,L}) ) - 0.5*dU_{i}
//    U_{i,R} = U^{n}_{i} + 0.5* dU_{i}
//    U_{i,L} = U^{n}_{i} - 0.5* dU_{i}
//    dU_{i} is some kind of limited slope.

// Hancock version 2
//
// U^{n-1/2}_{i-1/2,L}= U^{n}_{i-1} - 0.5*dtdx*dFdX(U_{i-1}) + 0.5*dU_{i-1}
// U^{n-1/2}_{i-1/2,R} = U^{n}_{i} - 0.5*dtdx*dFdX(U_{i})  - 0.5*dU_{i}
// Here, dFdX is the analytic flux derivative based on the PRIMITIVE variables.
// This is used for testing a descrepancy with mhd_li.  


// As per my general convention, arrays of face centered variables start at -1/2
//
// map of indicies:   
//     space:         | i-2 | i-1 | i | i+1 | i+2 | 
//     index[];       |  3  |  1  | 0 |  2  |  4  |
//                    Swap i with j or k as needed.
//                    This swap is controlled by the Offset variable.
//                    cell(i  ,j+1,k  ) = i + nx*(j+1) + nx*ny(k)  = i + nx*j + ny*k + nx
//                                      = index[0] + nx
//                    cell(i  ,j  ,k+2) = i + nx*j + nx*ny*(k+2) 
//                                      = index[0] + 2*nx*ny 
//                    etc.
//    This map was chosen for extensibility:  if a new solver with a more complicated stencil is 
//    written, index[] can be added to without altering the structure of existing code.

int grid::MHD_Hancock1(float * Lhs, float * Rhs, int * index, int dim, float dT,
		       float Bx){
  int field; 

  float dtdx = dT * dxI[dim]; //dxI = 1/dx

  //dU_{i-1} and dU{i}
  float LeftDifference[MAX_MHD_WAVES +1], RightDifference[MAX_MHD_WAVES +1];

  //Temporary work space for U_{i-1,R}, U_{i-1,L}, U_{i,R}, U_{i,L}
  float CellLeftTmp[MAX_MHD_WAVES +1], CellRightTmp[MAX_MHD_WAVES +1];

  //The used quantites for U_{i-1,R}, U_{i-1,L}, U_{i,R}, U_{i,L}
  //This object houses and generates the flux.
  Fluid CellLeft, CellRight;

  //Compute dU_{i-1} (left) and dU_{i} (right.)  
  MHD_LinearSlope(LeftDifference,RightDifference,index,dim,TRUE);

  //Compute LHS from U_{i-1} and its differences.
  for( field=0;field<NumberReconstructed; field++){
    CellLeftTmp[field] =State[field][ index[1] ] - 0.5*LeftDifference[field];
    CellRightTmp[field]=State[field][ index[1] ] + 0.5*LeftDifference[field];
  }

  //If we don't have the Longitudinal field in the reconstruction,
  //add it here.
  if( AllCellCentered == 0 ){
    CellLeftTmp[ Sb[2] ] = Bx;
    CellRightTmp[ Sb[2] ] = Bx;
  }


  //Compute the fluxes.  (Note that this does some additional work, 
  //but since everything is in cache, it shouldn't matter.)
  CellLeft.Fill(CellLeftTmp);
  CellRight.Fill(CellRightTmp);

  for( field=0;field<NumberReconstructed;field++){
    Lhs[field] =  
      State[field][ index[1] ] 
      + 0.5*LeftDifference[field]
      - 0.5*dtdx * ( CellRight.Flux[field] - CellLeft.Flux[field] );
  }

  //Compute RHS from U_{i} and its differences.
  for(field=0;field<NumberReconstructed;field++){
    CellLeftTmp[field] = State[field][ index[0] ]-0.5*RightDifference[field];
    CellRightTmp[field]= State[field][ index[0] ]+0.5*RightDifference[field];
  }
  
  if( AllCellCentered == 0 ){
    CellLeftTmp[ Sb[2] ] = Bx;
    CellRightTmp[ Sb[2] ] = Bx;
  }
  
  
  //Compute fluxes from the state.
  CellLeft.Fill(CellLeftTmp);
  CellRight.Fill(CellRightTmp);
  
  for(field=0;field<NumberReconstructed;field++){
    Rhs[field] = 
      State[field][ index[0] ] 
      - 0.5*RightDifference[field]
      - 0.5*dtdx*( CellRight.Flux[field] - CellLeft.Flux[field] );
  }

  return SUCCESS;

}


int grid::MHD_Hancock2(float * Lhs, float * Rhs, int * index, int dim, float dT, float Bx){

  int field; 

  float dtdx = dT * dxI[dim]; //dxI = 1/dx

  //dU_{i-1} and dU{i}
  float LeftDifference[MAX_MHD_WAVES +1], RightDifference[MAX_MHD_WAVES +1];

  //dF/dx for both i-1 and i.
  float dFlux[MAX_MHD_WAVES +1];

  //Temporary work space for U_{i-1}, U_{i}
  float CellCenterTmp[MAX_MHD_WAVES +1];

  //The used quantity for U_{i-1}, U_{i}
  //This object houses and generates the flux. 
  Fluid CellCenter;

  //Compute dU_{i-1} (left) and dU_{i} (right.)  
  MHD_LinearSlope(LeftDifference,RightDifference,index,dim,TRUE);

  //Compute LHS from U_{i-1} and its differences.
  for( field=0;field<NumberReconstructed; field++){
    CellCenterTmp[field] =State[field][ index[1] ];
  }

  //If we don't have the Longitudinal field in the reconstruction,
  //add it here.
  if( AllCellCentered == 0 ){
    CellCenterTmp[ Sb[2] ] = Bx;
  }


  //Compute the fluxes.  (Note that this does some additional work, 
  //but since everything is in cache, it shouldn't matter.)
  CellCenter.Fill(CellCenterTmp);

  //Compute the flux difference.
  CellCenter.PrimitiveDerivative(dFlux, LeftDifference);

  for( field=0;field<NumberReconstructed;field++){
    Lhs[field] =  
      State[field][ index[1] ] 
      + 0.5*LeftDifference[field]
      - 0.5*dtdx * dFlux[field];
  }

  //Compute RHS from U_{i} and its differences.
  for(field=0;field<NumberReconstructed;field++){
    CellCenterTmp[field] = State[field][ index[0] ];
  }
  
  if( AllCellCentered == 0 ){
    CellCenterTmp[ Sb[2] ] = Bx;
  }
  
  
  //Compute fluxes from the state.
  CellCenter.Fill(CellCenterTmp);
  CellCenter.PrimitiveDerivative(dFlux,RightDifference);

  for(field=0;field<NumberReconstructed;field++){
    Rhs[field] =  
      State[field][ index[0] ] 
      - 0.5*RightDifference[field]
      - 0.5*dtdx*dFlux[field];
  }

  return SUCCESS;


}



#endif //ATHENA







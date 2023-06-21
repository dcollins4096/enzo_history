

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

#ifdef ATHENA

//
// Reconstruct.  
//
// Field wise, driven by MHD_PLM_SlopeLocal[field],
// which is derived from the input variable MHD_PLM_Slope
// 1:  piecewise constant (default)
// 2:  min mod
// 3: van leer.

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


//                   min( |dL|, |dR|) * sign(dL) if  sign(dR)=sign(dL)
// minmod( dL, dR) = 0                           if sign(dR) != sign(dL)
//                   So, the one closer to zero, provided they have they same sign.
//                   Using the smallest slope supresses post-shock oscillation.

//
// VanLeerMC( dL, dC, dR) = minmod( 2*dL, dC, 2*dR)
//                        Again, the one closest to zero, provided they have the same sign.
//                        Using the smallest slope supresses post-shock oscillation.
//                        This slope also satisfies the Upwind Range condition, a stronger condition.
//                        See Laney (2001) "Computational Gasdynamics," sections 16.5 and 23.1 for details.
//                        Also Van Leer (1977) "Towards the Ultimate ... III, IV", J.comp. phys.
//          
inline float MinMod(float dL,float dR){ 
  return ( dL*dR < 0 ) ? 0  : ( ( fabs( dL ) < fabs( dR ) ) ? dL : dR ) ;
}

inline float VanLeerMC(float dL, float dC, float dR){
  if( (dL > 0 && dC > 0 && dR > 0 ) ||
      (dL < 0 && dC < 0 && dR < 0 ) ){
    return MinMod( dC, MinMod( 2*dR, 2*dL) );
  }else{
    return 0;
  }

}
int grid::MHD_PLM(float * Lhs, float * Rhs, int * index, int dim){
  int field; 

  for( field=0; field< NumberReconstructed; field ++ ){


    switch( MHD_PLM_SlopeLocal[field] ){
      
    case 1:

      //Piecewise Constant.
      Lhs[ field ] = State[field][ index[1] ];
      Rhs[ field ] = State[field][ index[0] ];

      break;
      
    case 2:
      
      //     space:         | i-2 | i-1 | i | i+1 | i+2 | 
      //     index[];       |  3  |  1  | 0 |  2  |  4  |
      
      //Min Mod limited slope      
      Lhs[ field ] = State[field][ index[1] ] + 0.5 * 
	MinMod( State[field][ index[1] ] - State[field][ index[3] ],
		State[field][ index[0] ] - State[field][ index[1] ] );
      Rhs[ field ] = State[field][ index[0] ] - 0.5 * 
	MinMod( State[field][ index[0] ] - State[field][ index[1] ],
		State[field][ index[2] ] - State[field][ index[0] ] );
      
      
      break;
    case 3:

      //     space:         | i-2 | i-1 | i | i+1 | i+2 | 
      //     index[];       |  3  |  1  | 0 |  2  |  4  |
      
      // Van Leer "monotonized central" limiter.
      Lhs[ field ] = State[field][ index[1] ] + 0.5 * 
	VanLeerMC( (State[field][ index[1] ] - State[field][ index[3] ]),
		   (State[field][ index[0] ] - State[field][ index[3] ])*0.5,
		   (State[field][ index[0] ] - State[field][ index[1] ]) );
      Rhs[ field ] = State[field][ index[0] ] - 0.5 * 
	VanLeerMC( (State[field][ index[0] ] - State[field][ index[1] ]),
		   (State[field][ index[2] ] - State[field][ index[1] ])*0.5,
		   (State[field][ index[2] ] - State[field][ index[0] ]) );
      
      break;
    default:

      fprintf(stderr," PLM Warning: Reconstruction method %d not defined\n",MHD_PLM_SlopeLocal[field]);
      return FAIL;
      break;
    }//slope switch
    
  }
  return SUCCESS;

}

#endif //ATHENA







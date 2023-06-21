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

int grid::MHD_Flatten(float * Lhs, float * Rhs, int * index, int dim){

  int field;    
  float F0, F1; //Max flattening for Rhs, Lhs (respectively)

  F0 = ( ( ShockDirection[ dim + GridRank*( index[0] ) ] == 1 ) ?
	 max(FlatteningField[ dim + GridRank*( index[0] ) ],
	     FlatteningField[ dim + GridRank*( index[2] ) ]) :
	 max(FlatteningField[ dim + GridRank*( index[0] ) ],
	     FlatteningField[ dim + GridRank*( index[1] ) ]) );
  
  F1 = ( ( ShockDirection[ dim + GridRank*( index[1] ) ] == 1 ) ?
	 max(FlatteningField[ dim + GridRank*( index[1] ) ],
	     FlatteningField[ dim + GridRank*( index[0] ) ]) :
	 max(FlatteningField[ dim + GridRank*( index[1] ) ],
	     FlatteningField[ dim + GridRank*( index[3] ) ]) );
  

  //     space:         | i-2 | i-1 | i | i+1 | i+2 | 
  //     index[];       |  3  |  1  | 0 |  2  |  4  |

  for( field=0; field< NumberReconstructed; field ++ ){
			    
    switch( MHD_Flattening ){

    case 1:
    case 3:
      //Rhs and Lhs already have the high order interpolated data.  
      Rhs[ field ] = F0 * State[field][index[0]] +(1-F0)*Rhs[field];
      Lhs[ field ] = F1 * State[field][index[1]] + (1-F1)*Lhs[field];
      
      break;

    case 0:
    default:
      fprintf(stderr," Undefined flattening method %d\n", MHD_Flattening);
      return FAIL;
      
    }//flattening switch

    }//field
  
  return SUCCESS;
}
#endif //ATHENA

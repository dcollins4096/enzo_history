#ifdef PPML
//
// RotateState( dim ) rotates the vector fields along the appropriate axis.
//
// The Athena solver is formulated always in the Longitudinal Direction.
// This saves us from writing reconstructions and riemann solvers for all 3 dimensions,
// AND from coppies from the 3d cube to 2d stripes.  
//
// The thought is that since this is an unsplit scheme, and all data needs to be compued
// at the same time, we'll save memory bandwidth by doing the entire cell at once.
// 
//
// So for the X direction, V = {vx,vy,vz},
//    for the Y direction, V = {vy,vz,vx}, 
//    etc.  
//
//    Sv[0,1,2] = kinetic field (velocity or momentum) that the Solver sees.
//    Ev[0,1,2] = kinetic field that Enzo sees.
//    MapEtoS[dim][0,1,2] = the map (which is dim dependant) from one to the other.
//    Sb[0,1]   = transvers field components for the Solver
//    BNum[dim][0,1] = transverse fields from Enzo.

#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Athena.h"

//inline me.
int Athena::RotateState(int dim){
  int field;

  //Temporary: might not want to be here after all.
  State[ MapEtoS[dim][ E.D ] ] = BaryonField[ E.D ];
  if( EquationOfState == 0 )
      State[ MapEtoS[dim][ E.TE ] ] = BaryonField[ E.TE ];

  for( field=0; field<3; field++){
    State[ MapEtoS[dim][ E.V[field] ] ] = BaryonField[ E.V[field] ];
    State[ MapEtoS[dim][ E.B[field] ] ] = BaryonField[ E.B[field] ];
  }

#ifdef DoThisLater
#ifdef MHD
  for( field=0; field<2; field++)
    State[ S.B[field] ] = CenteredB[ BNum[dim][field] ];
#endif //MHD
  if( HydroMethod == PPM_Local ){
    for( field=0; field<2; field++)
      State[ S.B[field] ] = [ BNum[dim][field] ];
    if( AllCellCentered == 1 ){
      State[ Sb[2] ] = CenteredB[dim];
    }
  }else{
    fprintf(stderr,"AtheanRotate: not functional with non PPML methods yet.\n");
    return FAIL;
  }
#endif //Later
  return SUCCESS;
}

#endif //PPML

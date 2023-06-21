

#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "MHD_Athena.h"
#include "MHDFluid.h"
#ifdef ATHENA

//
// Derivative of Primitive Flux.
// This is used initially to compare with the code used in mhd_li,
// which was crashing for isothermal runs.
// 
// Outcome TBD.

int Fluid::PrimitiveDerivative(float * dFlux, float * dU){

  float tau = 1./rho;
  float cs2 = ( (EquationOfState == 0 ) ? Gamma*pressure/rho : 
		IsothermalSoundSpeed*IsothermalSoundSpeed);
  //Density
  
  dFlux[Sden] = (vx*dU[Sden] + rho*dU[Sv[0]]);

  //Vx, vy, vz
  dFlux[Sv[0]] = (vx*dU[Sv[0]] + tau*(dU[Seng] + by*dU[Sb[0]] + bz*dU[Sb[2]]));
  dFlux[Sv[1]] = (vx*dU[Sv[1]] - bx*tau*dU[Sb[0]]);
  dFlux[Sv[2]] = (vx*dU[Sv[2]] - bx*tau*dU[Sb[2]]);

  //Bx, by, bz
  dFlux[Sb[2]] = 0;
  dFlux[Sb[0]] = (by*dU[Sv[0]] - bx*dU[Sv[1]] + dU[Sb[0]]*vx);
  dFlux[Sb[1]] = (bz*dU[Sv[0]] - bx*dU[Sv[2]] + dU[Sb[2]]*vx);

  //Pressure.  
  dFlux[Seng]  = (cs2*dU[Sv[0]] + vx*dU[Seng]);

  return SUCCESS;

}
#endif //ATHENA

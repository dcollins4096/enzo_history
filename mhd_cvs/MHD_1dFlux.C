//
// Just the MHD Flux.  
// Accounts for the various reconstruction methods available.
//  Different reconstructions on the Kinematic and Energy fields change
//  what variable is stored in Rhs.
// The state vector is called Rhs for historical reasons.
//
// Computed, in addition to the flux, are: 
//  EngR:      TotalEnergy 
//  MagPressR: Total Pressure, p_{gas} + 0.5*B^2
//  PressureR: Gas Pressure:
//  RhoEnthR:  Rho * Enthalpy
// These can be passed out through the last argument.
//
// (I really should Object this bit of the code a little bit. 
//  It's a little too Fortran-ie right now.)
//
// ALSO NOTE!!!!!
//  This routine currently changes the nature of the Kinematic components
//  of Rhs, so that if Momentum is reconstructed, it's converted to velocity.


#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "MHD_Athena.h"

#ifdef ATHENA
int MHD_1dFlux(float * Fluxes, float * Rhs){


  //variable we'll need:  See MHD_Roe1.C for documentation.
  float EngR, MagPressR, RhoEnthR, PressureR;
  int field;
  //
  // Set up derived quantities.
  // based on which reconstruction switch is on.

  // For kinetic-type field (velocity or momentum)
  // Controlled by MHD_ReconField[0]
  //     0: Reconstruct velocity
  //     1: Reconstruct momentum
  // Note that the whole rest of the solver is formulated with velocity, so remove density if
  // reconstructing momentum
  
  if( MHD_ReconField[0] == 1){
    for( field=0; field<3;field++){
      Rhs[Sv[field]]  /= Rhs[Sden];
    }
  }

  // For energy-type field (Total energy, pressure, total pressure, enthalpy):
  // Controled by MHD_ReconField[1].
  //     0: Reconstruct Total Energy
  //     1: Reconstruct gas Pressure
  //     2: Reconstruct Total Pressure
  //     3: Reconstruct Enthalpy

  if( EquationOfState == 0 ){
    switch( MHD_ReconField[1] ) {
      
    case 0:  //reconstruct energy
      //Already have EngL and EngR
      
      EngR = Rhs[Seng];
      
      //Compute Enthalpy, which is used fo Roe averaging.
      // H = (E + P + 1/2 B^2)/rho  &&  P = (E - 1/2V^2 - 1/2 B^2)*(Gamma -1)
      // so H = (Gamma E + 0.5*( (Gamma - 1) (rho v^2) + (Gamma -2 ) (B^2)))/rho
      
      //total pressure = p + 0.5 b^2 = (e - 0.5 rho v^2)*(gamma-1) + 0.5 b^2 * (Gamma -2)
      MagPressR =
	( EngR - 0.5*(Rhs[Sv[0]]*Rhs[Sv[0]] + Rhs[Sv[1]]*Rhs[Sv[1]] + Rhs[Sv[2]]*Rhs[Sv[2]])*Rhs[Sden])*(Gamma-1) - 
	0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]])*(Gamma-2);
      
      RhoEnthR = EngR+ MagPressR;
      break;
    case 1: //reconstruct gas pressure
      
      //Commented out because I don't actually NEED this variable.
      //PressureL = Lhs[Seng];
      //PressureR = Rhs[Seng];
      
      //TEMPORARY.  This isn't the full MagPress that's used in the rest of the solver.
      MagPressR = 0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]]);
    
      EngR = Rhs[Seng]/(Gamma - 1) + MagPressR + 
	0.5*(Rhs[Sv[0]]*Rhs[Sv[0]] + Rhs[Sv[1]]*Rhs[Sv[1]] + Rhs[Sv[2]]*Rhs[Sv[2]])*Rhs[Sden];
      
      MagPressR += Rhs[Seng];
      
      RhoEnthR = EngR+ MagPressR;    
      break;
    case 2: //reconstruct total pressure
      MagPressR = Rhs[Seng];
      
      PressureR = MagPressR -0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]]);
      
      //                            This part here just saves flops.
      EngR = PressureR/(Gamma -1) + 
	(MagPressR - PressureR) + 
	0.5*(Rhs[Sv[0]]*Rhs[Sv[0]] + Rhs[Sv[1]]*Rhs[Sv[1]] + Rhs[Sv[2]]*Rhs[Sv[2]])*Rhs[Sden];
      
      RhoEnthR = EngR+ MagPressR;
      
      break;
    case 3: //reconstruct enthalpy
      //Compute Enthalpy, which is used fo Roe averaging.
      // H = (E + P + 1/2 B^2)/rho  &&  P = (E - 1/2V^2 - 1/2 B^2)*(Gamma -1)
      // so H = (Gamma E + 0.5*( (Gamma - 1) (rho v^2) + (Gamma -2 ) (B^2)))/rho
      // or H = ( Gamma P/(Gamma - 1) + 0.5( rho v^2 + B^2 ) )/rho
      
      RhoEnthR = Rhs[Seng]*Rhs[Sden]; 
      
      //Temporary
      MagPressR = 0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]]);
      
      PressureR = (Gamma - 1)/Gamma * (RhoEnthR  - 2*MagPressR - 
				       0.5*(Rhs[Sv[0]]*Rhs[Sv[0]] + Rhs[Sv[1]]*Rhs[Sv[1]] + Rhs[Sv[2]]*Rhs[Sv[2]])*Rhs[Sden]);
      
      MagPressR += PressureR;
      
      EngR = RhoEnthR - MagPressR;
      
      break;
    }
  }else if( EquationOfState == 1){ //EOS

    MagPressR = 0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]])
      + IsothermalSoundSpeed*IsothermalSoundSpeed*Rhs[Sden];

  }
  
    //density flux = rho*vx

  Fluxes[Sden] = Rhs[Sden] * Rhs[Sv[0]];  

  if( EquationOfState == 0 ){
    //total energy flux  = vx * rho * enthalpy - bx* V \dot B
    //                                enthalpy = E + p
    Fluxes[Seng] = RhoEnthR * Rhs[Sv[0]] - Rhs[Sb[2]]*(Rhs[Sb[2]]*Rhs[Sv[0]]+ Rhs[Sb[0]]*Rhs[Sv[1]]+ Rhs[Sb[1]]*Rhs[Sv[2]]);
  }
  //X momentum = rho u^2 + pmag
  Fluxes[Sv[0]] = Rhs[Sden]*Rhs[Sv[0]]*Rhs[Sv[0]] + MagPressR - Rhs[Sb[2]]*Rhs[Sb[2]];

  //Y momentum = rho*vx*vy - Bx Sb[0]
  Fluxes[Sv[1]] = Rhs[Sden]*Rhs[Sv[0]]*Rhs[Sv[1]] - Rhs[Sb[2]]*Rhs[Sb[0]];

  //Z momentum = rho*vx*vz - Bx Bz
  Fluxes[Sv[2]] = Rhs[Sden]*Rhs[Sv[0]]*Rhs[Sv[2]] - Rhs[Sb[2]]*Rhs[Sb[1]];

  //Sb[0] flux = By Vx - Bx Vy
  Fluxes[Sb[0]] = Rhs[Sb[0]]*Rhs[Sv[0]] - Rhs[Sb[2]]*Rhs[Sv[1]];

  //Sb[1] flux = Bz Vx - Bx Vz
  Fluxes[Sb[1]] = Rhs[Sb[1]]*Rhs[Sv[0]] - Rhs[Sb[2]]*Rhs[Sv[2]];
  
  return SUCCESS;

}
#endif //ATHENA

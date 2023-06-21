


#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "MHD_Athena.h"
#include "MHDFluid.h"
#ifdef ATHENA

int Fluid::Fill(float * Recon ){
  
  int field;

  //The density field.
  rho = Recon[ Sden ];
  
  // For kinetic-type field (velocity or momentum)
  // Controlled by MHD_ReconField[0]
  //     0: Reconstruct velocity
  //     1: Reconstruct momentum
  
  if( MHD_ReconField[0] == 1){
    for( field=0; field<3;field++){
      momentum[field] = Recon[ Sv[field] ] ;
      vel[field] = Recon[ Sv[field] ] / rho;
    }
  }else if (MHD_ReconField[0] == 0){
    for( field=0; field<3;field++){
      momentum[field] = Recon[ Sv[field] ] * rho ;
      vel[field] = Recon[ Sv[field] ];
    }
  }//reconfield[0]

  //Assign magnetic field.  The numbering is due to the Flux not having any Bz.
  bx = Recon[Sb[2]];
  by = Recon[Sb[0]];
  bz = Recon[Sb[1]];

  // For energy-type field (Total energy, pressure, total pressure, enthalpy):
  // Controled by MHD_ReconField[1].
  //     0: Reconstruct Total Energy
  //     1: Reconstruct gas Pressure
  //     2: Reconstruct Total Pressure
  //     3: Reconstruct Enthalpy

  if( EquationOfState == 0 ){
    switch( MHD_ReconField[1] ) {
      
    case 0:  //reconstruct energy
      //Already have EngL and Eng
      
      Eng = Recon[Seng];
      
      //Compute Enthalpy, which is used fo Roe averaging.
      // H = (E + P + 1/2 B^2)/rho  &&  P = (E - 1/2V^2 - 1/2 B^2)*(Gamma -1)
      // so H = (Gamma E + 0.5*( (Gamma - 1) (rho v^2) + (Gamma -2 ) (B^2)))/rho
      
      //total pressure = p + 0.5 b^2 = (e - 0.5 rho v^2)*(gamma-1) + 0.5 b^2 * (Gamma -2)
      press_tot =
	( Eng - 0.5*(vx*vx + vy*vy + vz*vz)*rho)*(Gamma-1) - 
	0.5*(   bx*bx    + by*by + bz*bz)*(Gamma-2);
      
      enth = Eng+ press_tot;
      break;
    case 1: //reconstruct gas pressure
      
      //Commented out because I don't actually NEED this variable.
      //PressureL = Lhs[Seng];
      //pressure = Recon[Seng];
      
      //TEMPORARY.  This isn't the full MagPress that's used in the rest of the solver.
      press_tot = 0.5*(   bx*bx + by*by + bz*bz);
    
      Eng = Recon[Seng]/(Gamma - 1) + press_tot + 
	0.5*(vx*vx + vy*vy + vz*vz)*rho;
      
      press_tot += Recon[Seng];
      
      enth = Eng + press_tot;    
      break;
    case 2: //reconstruct total pressure
      press_tot = Recon[Seng];
      
      pressure = press_tot -0.5*( bx*bx + by*by + bz*bz);
      
      Eng = pressure/(Gamma -1) + 
	(press_tot - pressure) + 
	0.5*(vx*vx + vy*vy + vz*vz)*rho;
      
      enth = Eng+ press_tot;
      
      break;
    case 3: //reconstruct enthalpy
      //Compute Enthalpy, which is used fo Roe averaging.
      // H = (E + P + 1/2 B^2)/rho  &&  P = (E - 1/2V^2 - 1/2 B^2)*(Gamma -1)
      // so H = (Gamma E + 0.5*( (Gamma - 1) (rho v^2) + (Gamma -2 ) (B^2)))/rho
      // or H = ( Gamma P/(Gamma - 1) + 0.5( rho v^2 + B^2 ) )/rho
      
      enth = Recon[Seng]*rho; 
      
      //Temporary
      press_tot = 0.5*( bx*bx + by*by + bz*bz);
      
      pressure = (Gamma - 1)/Gamma * (enth  - 2*press_tot - 
				      0.5*(vx*vx + vy*vy + vz*vz)*rho);
      
      press_tot += pressure;
      
      Eng = enth - press_tot;
      
      break;
    }
  }else if( EquationOfState == 1){ //EOS

    press_tot = 0.5*(   bx*bx    + by*by + bz*bz)
      + IsothermalSoundSpeed*IsothermalSoundSpeed*rho;

  }
  
  //density flux = rho*vx

  Flux[Sden] = rho * vx;  

  if( EquationOfState == 0 ){
    //total energy flux  = vx * rho * enthalpy - bx* V \dot B
    //                                enthalpy = E + p
    Flux[Seng] = enth * vx - bx*(bx*vx+ by*vy+ bz*vz);
  }
  //X momentum = rho u^2 + pmag - bx^2
  Flux[Sv[0]] = rho*vx*vx + press_tot - bx*bx;

  //Y momentum = rho*vx*vy - Bx Sb[0]
  Flux[Sv[1]] = rho*vx*vy - bx*by;

  //Z momentum = rho*vx*vz - Bx Bz
  Flux[Sv[2]] = rho*vx*vz - bx*bz;

  //Sb[0] flux = By Vx - Bx Vy
  Flux[Sb[0]] = by*vx - bx*vy;

  //Sb[1] flux = Bz Vx - Bx Vz
  Flux[Sb[1]] = bz*vx - bx*vz;

  //Sb[2] flux = 0.  No Bx flux.
  Flux[Sb[2]] = 0;
  
  //Fill the vector of conserved quantites.
  cons[ Sden  ] = rho;
  cons[ Sv[0] ] = mx;
  cons[ Sv[1] ] = my;
  cons[ Sv[2] ] = mz;
  if( EquationOfState == 0 ) cons[ Seng ]  = Eng;
  cons[ Sb[0] ] = by;
  cons[ Sb[1] ] = bz;

  //<dbg>
  if( rho < 0.0 ){
    fprintf(stderr,"MHDFluid_Fill: Negative Density!\n");
  }
  //</dbg>
  return SUCCESS;


}
#endif //ATHERNA



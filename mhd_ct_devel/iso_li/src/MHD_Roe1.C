//
// MHD Roe solver, based on the work of Cargo & Gallice, JCP 1997
//

//
// Code flow:
// 1.) Compute averages & differences needed.
// 2.) Compute wave speeds, 
// 3.) Compute right eigen vecotrs
// 4.) Compute Left eigen vectors/ Wave Strength
// 5.) Compute net flux.

#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "MHD_Athena.h"
#include "MHDFluid.h"
void esys_roe_adb_mhd(const float d, const float v1, const float v2, const float v3,
		      const float h, const float b1, const float b2, const float b3, 
		      const float x, const float y,
		      float eigenvalues[],
		      float right_eigenmatrix[][7], float left_eigenmatrix[][7]);

void esys_roe_iso_mhd(const float d, const float v1, const float v2, const float v3,
		      const float Iso_csound2, const float b1, const float b2, const float b3, 
		      const float x, const float y, 
		      float eigenvalues[],
		      float right_eigenmatrix[][7], float left_eigenmatrix[][7]);


#ifdef ATHENA


int MHD_ROE_NEW(float * Fluxes, Fluid * L, Fluid * R ){
  int field, wave;

  //Eigen System.  Note that RhsV and LhsV, and Strength should be removed
  //after I've tested this shit.

  float RhsV[MAX_MHD_WAVES][MAX_MHD_WAVES], LhsV[MAX_MHD_WAVES][MAX_MHD_WAVES];
  float Speeds[MAX_MHD_WAVES], Strength[MAX_MHD_WAVES];
  float MidState[MAX_MHD_WAVES];

  //Failure checks: EjectionSeat checks for nans.  UseOtherSolver checks for density<0.
  int EjectionSeat = 0;
  int UseOtherSolver = 0;
  int ReturnValue = SUCCESS;

  // Components of Left and Right vectors MUST be zeroed.
  // Note the Left Vectors are stored as rows, right are stored as columns.
  for( field = 0; field< NumberOfMHDFluxes; field++){
    for( wave=0; wave<NumberOfMHDFluxes; wave++){
      RhsV[field][wave] = 0;
      LhsV[wave][field] = 0;
    }
  }


  //Roe Averages
  float sqrtRhoL = sqrt( L->rho );
  float sqrtRhoR = sqrt( R->rho );
  float overSqrts = 1.0/(sqrtRhoL + sqrtRhoR);
  float rhoBar = sqrtRhoL*sqrtRhoR;
  float vBar[3];
  float BxBar, ByBar, BzBar;
  float EnthBar;
  float DU[MAX_MHD_WAVES];
  float XX,YY;

  for(field=0;field<3;field++)
    vBar[field] = (sqrtRhoL * L->vel[field]  + sqrtRhoR * R->vel[field] )*overSqrts;

  // Not sure what to do with this code once HLL average of Bx comes in.
  if( AllCellCentered == 1){
    BxBar = (sqrtRhoL * R->bx + sqrtRhoR * L->bx)*overSqrts;
  }else{
    BxBar = L->bx;  //Or R.bx, they're the same for face centered runs.
  }

  ByBar = (sqrtRhoL * R->by + sqrtRhoR * L->by )*overSqrts;
  BzBar = (sqrtRhoL * R->bz + sqrtRhoR * L->bz )*overSqrts;

  if( EquationOfState == 0){
    EnthBar = ( L->enth/sqrtRhoL + R->enth/sqrtRhoR ) *overSqrts;
  }

  for( field=0;field<NumberOfMHDFluxes;field++)
    DU[field] = R->cons[field] - L->cons[field];

  XX   = 0.5*(DU[Sb[0]]*DU[Sb[0]] + DU[Sb[1]]*DU[Sb[1]])*overSqrts*overSqrts;
  YY   = 0.5*(L->rho + R->rho)/rhoBar;

  //
  // Call eigensystem
  //
  if( EquationOfState == 0 ){
    esys_roe_adb_mhd(rhoBar, vBar[0], vBar[1], vBar[2], EnthBar, BxBar, ByBar,BzBar, XX, YY, Speeds, RhsV, LhsV);
  }else if (EquationOfState == 1 ){
    esys_roe_iso_mhd(rhoBar, vBar[0], vBar[1], vBar[2], IsothermalSoundSpeed*IsothermalSoundSpeed, 
		     BxBar, ByBar,BzBar, XX, YY, Speeds, RhsV, LhsV);
  }


  for( wave =0; wave < NumberOfMHDFluxes; wave++){

    Strength[wave] = 0.0;

    for( field=0; field< NumberOfMHDFluxes; field++)
      Strength[wave] += LhsV[wave][field]*DU[field];
  }

  //Find density & pressure in all 6 intermediate states.
  //If any are negative, eject and call HLLE or HLLC.

  for(field = 0; field<NumberOfMHDFluxes; field++) 
    MidState[field] = L->cons[field];

  for( wave=0;wave < NumberOfMHDFluxes; wave++){
    for( field=0; field< NumberOfMHDFluxes; field++){
      MidState[field] += Strength[wave]*RhsV[field][wave];
    }//field
    if( MidState[ Sden ] <= 0 ){
      UseOtherSolver = TRUE;
      //<dbg>
      //fprintf(stderr,"Roe Solver says  d < 0!\n");
      //</dbg>
      break;
    }
 
    //compare pressure to zero
    if( EquationOfState == 0 ){
      if( ( MidState[ Seng ] - 0.5*(MidState[ Sv[0] ]*MidState[ Sv[0] ] +
				    MidState[ Sv[1] ]*MidState[ Sv[1] ] +
				    MidState[ Sv[2] ]*MidState[ Sv[2] ] )/MidState[ Sden ]
	    -0.5*(MidState[ Sb[0] ]*MidState[ Sb[0] ] +
		  MidState[ Sb[0] ]*MidState[ Sb[0] ] +
		  MidState[ Sb[0] ]*MidState[ Sb[0] ] ) ) < 0.0 ){
	UseOtherSolver = TRUE;
	//<dbg>
	//fprintf(stderr,"Hey, man, fuck!  p < 0!\n");
	//</dbg>
	break;
      }
    }//EOS
  }//wave
  
  if( UseOtherSolver == TRUE ){
    ReturnValue = 2;
  }else{
    ReturnValue = SUCCESS;
  }    
  for(field=0;field<NumberOfMHDFluxes; field++){
    Fluxes[field] = 0.5*( L->Flux[field] + R->Flux[field] );
    for( wave =0; wave < NumberOfMHDFluxes; wave++)
      Fluxes[field] -= 0.5*RhsV[field][wave]*fabs(Speeds[wave])*Strength[wave];

    if( Fluxes[field] != Fluxes[field] ) {
      EjectionSeat = 1;
      fprintf(stderr,"Bad flux in MHD_Roe1 solver. Abort. %d\n", field);
    }

  }//field

  if( EjectionSeat == 1 ){
    if( SuggestFailure == FALSE ){
      for( int j=0; j<10;j++)
	WriteInThisA[j] = 30 + j; 
      fprintf(stderr," Roe suggesting failue.\n");
      MHD_WriteDivB = FALSE; //I know that divB != 0, we're rejecting from the middle of the solve.
      fprintf(stderr," n Speeds:  Fluxes: FluxesL: FluxesR\n");
      for( wave=0; wave<NumberOfMHDFluxes;wave++)
	fprintf(stderr," %d %f %f %f %f\n", wave, Speeds[wave], Fluxes[wave], L->Flux[wave],R->Flux[wave]);
      
    }//suggestion
    SuggestFailure = TRUE;
    ReturnValue = FAIL;
  }//eject

  return ReturnValue;//SUCCESS;

}

//
//
//
//  Don't use this one.  It's only here for historical purposes.
//
//
int MHD_Roe1(float * Fluxes, float * Lhs, float * Rhs ){

  //
  // Set up variables
  //

  // Lhs[i] = {rho, vx, vy, vz, by, bz, energy}  
  // declare eigen system: right and left vector sets, eigen values, strengths ( Lhs \dot \Delta U)
  float RhsV[MAX_MHD_WAVES][MAX_MHD_WAVES], LhsV[MAX_MHD_WAVES][MAX_MHD_WAVES];
  float Speeds[MAX_MHD_WAVES], Strength[MAX_MHD_WAVES];
  float FluxL[MAX_MHD_WAVES], FluxR[MAX_MHD_WAVES];

  // Averages: Density, velocities, magnetic fields, enthalpy.  Also Left and Rhs Enthalpy. (actuall, rho*enthalpy)
  float rho, vx, vy, vz, bx, by, bz, enth;

  //Left and right energy-type quantites.  There's some slight redundancy between these quantities
  // and the input states due to the general nature of the switch mechanism installed.
  float RhoEnthL, RhoEnthR, PressureL, PressureR, MagPressL, MagPressR, EngL, EngR;

  // Differences: Density, energy, velocities, magnetic fields
  float DU[MAX_MHD_WAVES];

  // Factors from CG97 (XX,YY), anscillary variables.
  float XX, YY;
  float sqrtRhoL, sqrtRhoR, overSqrts; 
  int field, wave;

  // Components of Left and Right vectors MUST be zeroed.
  // Note the Left Vectors are stored as rows, right are stored as columns.
  for( field = 0; field< NumberOfMHDFluxes; field++){
    for( wave=0; wave<NumberOfMHDFluxes; wave++){
      RhsV[field][wave] = 0;
      LhsV[wave][field] = 0;
    }
  }
  
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
      Lhs[Sv[field]]  /= Lhs[Sden];
    }
  }
  // For energy-type field (Total energy, pressure, total pressure, enthalpy):
  // Controled by MHD_ReconField[1].
  //     0: Reconstruct Total Energy
  //     1: Reconstruct gas Pressure
  //     2: Reconstruct Total Pressure
  //     3: Reconstruct Enthalpy


  if( EquationOfState == 0 ) {
    switch( MHD_ReconField[1] ) {
      
    case 0:  //reconstruct energy
      //Already have EngL and EngR
      
      EngL = Lhs[Seng];
      EngR = Rhs[Seng];
      
      //Compute Enthalpy, which is used fo Roe averaging.
      // H = (E + P + 1/2 B^2)/rho  &&  P = (E - 1/2V^2 - 1/2 B^2)*(Gamma -1)
      // so H = (Gamma E + 0.5*( (Gamma - 1) (rho v^2) + (Gamma -2 ) (B^2)))/rho
      
      
      //total pressure = p + 0.5 b^2 = (e - 0.5 rho v^2)*(gamma-1) + 0.5 b^2 * (Gamma -2)
      MagPressL = ( EngL - 0.5*(Lhs[Sv[0]]*Lhs[Sv[0]] + Lhs[Sv[1]]*Lhs[Sv[1]] + Lhs[Sv[2]]*Lhs[Sv[2]])*Lhs[Sden])*(Gamma-1) - 
	0.5*(   Lhs[Sb[2]]*Lhs[Sb[2]]    + Lhs[Sb[0]]*Lhs[Sb[0]] + Lhs[Sb[1]]*Lhs[Sb[1]]) *(Gamma-2);
      
      MagPressR = ( EngR - 0.5*(Rhs[Sv[0]]*Rhs[Sv[0]] + Rhs[Sv[1]]*Rhs[Sv[1]] + Rhs[Sv[2]]*Rhs[Sv[2]])*Rhs[Sden])*(Gamma-1) - 
	0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]])  *(Gamma-2);

      RhoEnthL = EngL+ MagPressL;
      RhoEnthR = EngR+ MagPressR;
      break;
    case 1: //reconstruct gas pressure
      
      //Commented out because I don't actually NEED this variable.
      //PressureL = Lhs[Seng];
      //PressureR = Rhs[Seng];
      
      //TEMPORARY.  This isn't the full MagPress that's used in the rest of the solver.
      MagPressL = 0.5*(   Lhs[Sb[2]]*Lhs[Sb[2]]    + Lhs[Sb[0]]*Lhs[Sb[0]] + Lhs[Sb[1]]*Lhs[Sb[1]]);
      MagPressR = 0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]]);
      
      EngL = Lhs[Seng]/(Gamma - 1) + MagPressL +  
	0.5*(Lhs[Sv[0]]*Lhs[Sv[0]] + Lhs[Sv[1]]*Lhs[Sv[1]] + Lhs[Sv[2]]*Lhs[Sv[2]])*Lhs[Sden];
      EngR = Rhs[Seng]/(Gamma - 1) + MagPressR + 
	0.5*(Rhs[Sv[0]]*Rhs[Sv[0]] + Rhs[Sv[1]]*Rhs[Sv[1]] + Rhs[Sv[2]]*Rhs[Sv[2]])*Rhs[Sden];
      
      MagPressL += Lhs[Seng];
      MagPressR += Rhs[Seng];
      
      RhoEnthL = EngL+ MagPressL;
      RhoEnthR = EngR+ MagPressR;    
      break;
    case 2: //reconstruct total pressure
      MagPressL = Lhs[Seng];
      MagPressR = Rhs[Seng];
      
      PressureL = MagPressL -0.5*(   Lhs[Sb[2]]*Lhs[Sb[2]]    + Lhs[Sb[0]]*Lhs[Sb[0]] + Lhs[Sb[1]]*Lhs[Sb[1]]);
      PressureR = MagPressR -0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]]);

      //                            This part here just saves flops.
      EngL = PressureL/(Gamma -1 ) +
	(MagPressL - PressureL) + 
	0.5*(Lhs[Sv[0]]*Lhs[Sv[0]] + Lhs[Sv[1]]*Lhs[Sv[1]] + Lhs[Sv[2]]*Lhs[Sv[2]])*Lhs[Sden];
      EngR = PressureR/(Gamma -1) + 
	(MagPressR - PressureR) + 
	0.5*(Rhs[Sv[0]]*Rhs[Sv[0]] + Rhs[Sv[1]]*Rhs[Sv[1]] + Rhs[Sv[2]]*Rhs[Sv[2]])*Rhs[Sden];
      
      RhoEnthL = EngL+ MagPressL;
      RhoEnthR = EngR+ MagPressR;
      
      break;
    case 3: //reconstruct enthalpy
      //Compute Enthalpy, which is used fo Roe averaging.
      // H = (E + P + 1/2 B^2)/rho  &&  P = (E - 1/2V^2 - 1/2 B^2)*(Gamma -1)
      // so H = (Gamma E + 0.5*( (Gamma - 1) (rho v^2) + (Gamma -2 ) (B^2)))/rho
      // or H = ( Gamma P/(Gamma - 1) + 0.5( rho v^2 + B^2 ) )/rho
      
      RhoEnthL = Lhs[Seng]*Lhs[Sden]; 
      RhoEnthR = Rhs[Seng]*Rhs[Sden]; 
      
      //Temporary
      MagPressL = 0.5*(   Lhs[Sb[2]]*Lhs[Sb[2]]    + Lhs[Sb[0]]*Lhs[Sb[0]] + Lhs[Sb[1]]*Lhs[Sb[1]]);
      MagPressR = 0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]]);
      
      PressureL = (Gamma - 1)/Gamma * (RhoEnthL  - 2*MagPressL - 
				       0.5*(Lhs[Sv[0]]*Lhs[Sv[0]] + Lhs[Sv[1]]*Lhs[Sv[1]] + Lhs[Sv[2]]*Lhs[Sv[2]])*Lhs[Sden]);
      PressureR = (Gamma - 1)/Gamma * (RhoEnthR  - 2*MagPressR - 
				       0.5*(Rhs[Sv[0]]*Rhs[Sv[0]] + Rhs[Sv[1]]*Rhs[Sv[1]] + Rhs[Sv[2]]*Rhs[Sv[2]])*Rhs[Sden]);
      
      MagPressL += PressureL;
      MagPressR += PressureR;
      
      EngL = RhoEnthL - MagPressL;
      EngR = RhoEnthR - MagPressR;
      
      break;
    }
  }else if( EquationOfState == 1){ //EOS
    MagPressL = 0.5*(   Lhs[Sb[2]]*Lhs[Sb[2]]    + Lhs[Sb[0]]*Lhs[Sb[0]] + Lhs[Sb[1]]*Lhs[Sb[1]])
      + IsothermalSoundSpeed*IsothermalSoundSpeed*Lhs[Sden];
    
    MagPressR = 0.5*(   Rhs[Sb[2]]*Rhs[Sb[2]]    + Rhs[Sb[0]]*Rhs[Sb[0]] + Rhs[Sb[1]]*Rhs[Sb[1]])
      + IsothermalSoundSpeed*IsothermalSoundSpeed*Rhs[Sden];
  }


  //
  // fill roe averages
  //
  sqrtRhoL = sqrt( Lhs[Sden] );
  sqrtRhoR = sqrt( Rhs[Sden]);
  overSqrts = 1.0/(sqrtRhoL + sqrtRhoR );

  rho = sqrtRhoL * sqrtRhoR;
  vx  = (sqrtRhoL * Lhs[Sv[0]] + sqrtRhoR * Rhs[Sv[0]])*overSqrts;
  vy  = (sqrtRhoL * Lhs[Sv[1]] + sqrtRhoR * Rhs[Sv[1]])*overSqrts;
  vz  = (sqrtRhoL * Lhs[Sv[2]] + sqrtRhoR * Rhs[Sv[2]])*overSqrts;

  if( AllCellCentered == 1 ){
    bx  = (sqrtRhoL * Rhs[Sb[2]] + sqrtRhoR * Lhs[Sb[2]] ) *overSqrts;
  }else{
  bx  = Lhs[Sb[2]];  //no discontinuity here
  }
  by  = (sqrtRhoL * Rhs[Sb[0]] + sqrtRhoR * Lhs[Sb[0]] ) *overSqrts;
  bz  = (sqrtRhoL * Rhs[Sb[1]] + sqrtRhoR * Lhs[Sb[1]] ) *overSqrts;

  if( EquationOfState == 0 ){
    enth= ( RhoEnthL/sqrtRhoL + RhoEnthR/sqrtRhoR )*overSqrts;
  }

  // Differences.  This is a little clumsy, but Rhs[1,2,3] are Velocities,
  // and the method wants Momenta.

  DU[Sden]  = Rhs[Sden] - Lhs[Sden];
  DU[Sv[0]] = Rhs[Sv[0]]*Rhs[Sden] - Lhs[Sv[0]]*Lhs[Sden];
  DU[Sv[1]] = Rhs[Sv[1]]*Rhs[Sden] - Lhs[Sv[1]]*Lhs[Sden];
  DU[Sv[2]] = Rhs[Sv[2]]*Rhs[Sden] - Lhs[Sv[2]]*Lhs[Sden];
  DU[Sb[0]] = Rhs[Sb[0]] - Lhs[Sb[0]];
  DU[Sb[1]] = Rhs[Sb[1]] - Lhs[Sb[1]];
  if( EquationOfState == 0 )
    DU[Seng] = EngR - EngL;



  XX   = 0.5*(DU[Sb[0]]*DU[Sb[0]] + DU[Sb[1]]*DU[Sb[1]])*overSqrts*overSqrts;
  YY   = 0.5*(Lhs[Sden]+Rhs[Sden])/rho;

  //
  // Call eigensystem
  //
  if( EquationOfState == 0 ){
    esys_roe_adb_mhd(rho, vx, vy, vz, enth, bx, by, bz, XX, YY, Speeds, RhsV, LhsV);
  }else if (EquationOfState == 1 ){
    esys_roe_iso_mhd(rho, vx, vy, vz, IsothermalSoundSpeed, bx, by, bz, XX, YY, Speeds, RhsV, LhsV);
  }
  //
  // Compute fluxes.  Note: this is some terrible notation.  It says FluxL[Sv[0]], but it's actually a MOMENTUM flux.
  //
  
  //density flux = rho*vx

  FluxL[Sden] = Lhs[Sden] * Lhs[Sv[0]];   
  FluxR[Sden] = Rhs[Sden] * Rhs[Sv[0]];  

  if( EquationOfState == 0 ){
    //total energy flux  = vx * rho * enthalpy - bx* V \dot B
    //                                enthalpy = E + p
    FluxL[Seng] = RhoEnthL * Lhs[Sv[0]] - Lhs[Sb[2]]*(Lhs[Sb[2]]*Lhs[Sv[0]]+ Lhs[Sb[0]]*Lhs[Sv[1]]+ Lhs[Sb[1]]*Lhs[Sv[2]]);
    FluxR[Seng] = RhoEnthR * Rhs[Sv[0]] - Rhs[Sb[2]]*(Rhs[Sb[2]]*Rhs[Sv[0]]+ Rhs[Sb[0]]*Rhs[Sv[1]]+ Rhs[Sb[1]]*Rhs[Sv[2]]);
  }

  //X momentum = rho u^2 + pmag
  FluxL[Sv[0]] = Lhs[Sden]*Lhs[Sv[0]]*Lhs[Sv[0]] + MagPressL - Lhs[Sb[2]]*Lhs[Sb[2]];
  FluxR[Sv[0]] = Rhs[Sden]*Rhs[Sv[0]]*Rhs[Sv[0]] + MagPressR - Rhs[Sb[2]]*Rhs[Sb[2]];

  //Y momentum = rho*vx*vy - Bx Sb[0]
  FluxL[Sv[1]] = Lhs[Sden]*Lhs[Sv[0]]*Lhs[Sv[1]] - Lhs[Sb[2]]*Lhs[Sb[0]];
  FluxR[Sv[1]] = Rhs[Sden]*Rhs[Sv[0]]*Rhs[Sv[1]] - Rhs[Sb[2]]*Rhs[Sb[0]];

  //Z momentum = rho*vx*vz - Bx Bz
  FluxL[Sv[2]] = Lhs[Sden]*Lhs[Sv[0]]*Lhs[Sv[2]] - Lhs[Sb[2]]*Lhs[Sb[1]];
  FluxR[Sv[2]] = Rhs[Sden]*Rhs[Sv[0]]*Rhs[Sv[2]] - Rhs[Sb[2]]*Rhs[Sb[1]];

  //Sb[0] flux = By Vx - Bx Vy
  FluxL[Sb[0]] = Lhs[Sb[0]]*Lhs[Sv[0]] - Lhs[Sb[2]]*Lhs[Sv[1]];
  FluxR[Sb[0]] = Rhs[Sb[0]]*Rhs[Sv[0]] - Rhs[Sb[2]]*Rhs[Sv[1]];

  //Sb[1] flux = Bz Vx - Bx Vz
  FluxL[Sb[1]] = Lhs[Sb[1]]*Lhs[Sv[0]] - Lhs[Sb[2]]* Lhs[Sv[2]];
  FluxR[Sb[1]] = Rhs[Sb[1]]*Rhs[Sv[0]] - Rhs[Sb[2]]*Rhs[Sv[2]];

  //
  //Compute wave speeds
  //


  for( wave =0; wave < NumberOfMHDFluxes; wave++){

    Strength[wave] = 0.0;

    for( field=0; field< NumberOfMHDFluxes; field++)
      Strength[wave] += LhsV[wave][field]*DU[field];

  }

  int kludge = 0;
  int EjectionSeat = 0;
  if( kludge != 0 )fprintf(stderr,"Roe: ");
  for( field=0; field< NumberOfMHDFluxes; field++){
    
    Fluxes[field] = 0.5*(FluxL[field] + FluxR[field]);

    for( wave =0; wave < NumberOfMHDFluxes; wave++)
      Fluxes[field] -= 0.5*RhsV[field][wave]*fabs(Speeds[wave])*Strength[wave];

    if( Fluxes[field] != Fluxes[field] ) {
      EjectionSeat = 1;
      fprintf(stderr,"Bad flux in MHD_Roe1 solver. Abort. %d\n", field);
    }

    //Some debugging shit
    if( kludge != 0 )Fluxes[field]= FluxR[field];
    if( kludge != 0 )fprintf(stderr," %f ", Fluxes[field]);
  }//field
  
  if( EjectionSeat == 1 ){
    if( SuggestFailure == FALSE ){
      for( int j=0; j<10;j++)
	WriteInThisA[j] = 30 + j; 
      fprintf(stderr," Roe suggesting failue.\n");
      MHD_WriteDivB = FALSE; //I know that divB != 0, we're rejecting from the middle of the solve.
      fprintf(stderr," n Speeds:  Fluxes: FluxesL: FluxesR\n");
      for( wave=0; wave<NumberOfMHDFluxes;wave++)
	fprintf(stderr," %d %f %f %f %f\n", wave, Speeds[wave], Fluxes[wave], FluxL[wave],FluxR[wave]);
      
    }

      
    SuggestFailure = TRUE;
  }

  if( kludge != 0 )fprintf(stderr,"\n");
  return SUCCESS;
}


//Left STOP



#endif


#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "MHD_Athena.h"

#include "MHDFluid.h"
#ifdef ATHENA


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

int MHD_HLLC(float * Fluxes, Fluid * L, Fluid * R ){

  //<dbg>
  //fprintf(stderr," HLLC solver entered  \n");
  //</dbg>

  int field, wave;
  
  //Roe Averages, from Cargo & Galice roe paper.
  float sqrtRhoL = sqrt( L->rho );
  float sqrtRhoR = sqrt( R->rho );
  float overSqrts = 1.0/(sqrtRhoL + sqrtRhoR);
  float rhoBar = sqrtRhoL*sqrtRhoR;
  float vBar[3];
  float BxBar, ByBar, BzBar;
  float EnthBar;
  float DU[MAX_MHD_WAVES];
  float XX,YY;

  
  //Various speeds
  float Evals[MAX_MHD_WAVES];       //Eigen values from the roe matrix 
  float Lcf, Rcf, La2, Ra2, Lt, Rt; //fast mhd speed, sound speed, two temporary variables
  float Ls, Rs, Cs;                 //Outer (Left and Right) and Center (star) signal speed
  float OverSs;                     //Temp variable.

  //Intermediate States
  //HLLE is an approximation to the continuous central variables, B and V.
  //S4 and S5 are so named to relate to the states in Dai and Woodward that they most represent: 
  //the jump across the contact discontinuity.
  Fluid HLLE, S4, S5;
  float OverHLLERho;

  //
  //Compute the Roe Averages
  //


  for(field=0;field<3;field++)
    vBar[field] = (sqrtRhoL * L->vel[field]  + sqrtRhoR * R->vel[field] )*overSqrts;

  if( AllCellCentered == 1){
    BxBar = (sqrtRhoL * R->bx + sqrtRhoR * L->bx)*overSqrts;
  }else{
    BxBar = L->bx;  //Or R->bx, they're the same for face centered runs.
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
  // With NULL for the 2 vector arguments, the vectors aren't computed, only the wave speeds.
  //

  if( EquationOfState == 0 ){
    esys_roe_adb_mhd(rhoBar, vBar[0], vBar[1], vBar[2], EnthBar, BxBar, ByBar,BzBar, XX, YY, Evals, NULL,NULL);
  }else if (EquationOfState == 1 ){
    esys_roe_iso_mhd(rhoBar, vBar[0], vBar[1], vBar[2], IsothermalSoundSpeed,
		     BxBar, ByBar,BzBar, XX, YY, Evals, NULL,NULL);
  }

  //
  //compute Sl, Sr
  //
  
  if( EquationOfState == 0 ){
    Ra2 = Gamma * R->pressure / R->rho;
    La2 = Gamma * L->pressure / L->rho;
  }else if (EquationOfState == 1 ){
    Ra2 = IsothermalSoundSpeed;
    La2 = IsothermalSoundSpeed;
  }

  //Compute left and right Fast Speed
  Rt = Ra2 + (R->bx*R->bx + R->by*R->by + R->bz*R->bz)/R->rho;
  Lt = La2 + (L->bx*L->bx + L->by*L->by + L->bz*L->bz)/L->rho;

  Rcf = sqrt( 0.5*( Rt + sqrt( Rt*Rt - 4*Ra2*R->bx*R->bx/R->rho) ) );
  Lcf = sqrt( 0.5*( Lt + sqrt( Lt*Lt - 4*La2*L->bx*L->bx/L->rho) ) );

  Rs = max( Evals[ NumberOfMHDFluxes - 1] , R->vx + Rcf );
  Ls = min( Evals[0] , L->vx - Lcf );

  OverSs = 1.0/(Rs - Ls);

  //
  // Compute central states
  //


  //Compute u*  (Call it Cs) see Li, eqn (14)
  //note that there's a typo in that equation.  
  Cs =( R->rho*R->vx*(Rs - R->vx) - L->rho*L->vx*(Ls - L->vx) + 
 	(L->press_tot - L->bx*L->bx) - (R->press_tot - R->bx*R->bx) )/
    ( R->rho*(Rs - R->vx) - L->rho*( Ls - L->vx ) );
  
  //Compute HLLE state.  Note that this fluid is only partially initialized.
  //This state is used as an estimate for the underconstrained system of equations
  //See Toro (10.13) or Li (2)

  HLLE.rho = (Rs*R->rho - Ls*L->rho - (R->Flux[ Sden] - L->Flux[ Sden ])) *OverSs;
  OverHLLERho = 1.0/HLLE.rho;
  HLLE.vx =  (Rs*R->mx - Ls*L->mx - (R->Flux[ Sv[0] ] - L->Flux[ Sv[0] ]))*OverSs*OverHLLERho;
  HLLE.vy =  (Rs*R->my - Ls*L->my - (R->Flux[ Sv[1] ] - L->Flux[ Sv[1] ]))*OverSs*OverHLLERho;
  HLLE.vz =  (Rs*R->mz - Ls*L->mz - (R->Flux[ Sv[2] ] - L->Flux[ Sv[2] ]))*OverSs*OverHLLERho;
  HLLE.bx =  (Rs*R->bx - Ls*L->bx )*OverSs;
  HLLE.by =  (Rs*R->by - Ls*L->by - (R->Flux[ Sb[0] ] - L->Flux[ Sb[0] ]))*OverSs;
  HLLE.bz =  (Rs*R->bz - Ls*L->bz - (R->Flux[ Sb[1] ] - L->Flux[ Sb[1] ]))*OverSs;

  //
  //Fill from two internal states
  //

  S4.HLLC_Star(L,&HLLE, Ls,Cs);
  S5.HLLC_Star(R,&HLLE, Rs,Cs);

  //speed switch:
  // F_{hllc} = 
  // Ls > 0 ? L.flux              // very right going
  // Cs > 0 ? L.flux + Ls*(S4 - L)// right going 
  // Cs < 0 ? R.flux + Rs*(S5 - R)// left going
  // Rs < 0 ? R.flux              // very left going

  if( Ls >= 0 ){
    for(field=0;field<NumberOfMHDFluxes; field++)
      Fluxes[field] = L->Flux[field];
  }else if( Cs >= 0 && Ls < 0){
    for(field=0;field<NumberOfMHDFluxes; field++)
      Fluxes[field] = L->Flux[field] + Ls*(S4.cons[field] - L->cons[field] );
  }else if( Cs <= 0 && Rs >= 0){
    for(field=0;field<NumberOfMHDFluxes; field++)
      Fluxes[field] = R->Flux[field] + Rs*(S5.cons[field] - R->cons[field] );
  }else {
    for(field=0;field<NumberOfMHDFluxes; field++)
      Fluxes[field] = R->Flux[field];
  }


  return SUCCESS;

}

// The HLLC_Stared (left and right) states, from Li, J Comp Phys 203 (2005) 344-357
// While the notation used here is for the Left side, it's symmetric for the right side.
// Also note that the EST object isn't completely filled, only with the variables i need.
// L is the Left(or right) state 
// "this" is the Left Star state (see Li (3) )
// To translate between this and Li (2005):
//   Unadorned variables in this code are "stared" variables in Li
//   L-> variables in this code are the unadorned variables in Li.
//   Marginally annoying, but details of the code drove this.
// EST is the estimate for the variables that are continuous across the center wave.
// In Li's case, this is the HLLE approximation.

int Fluid::HLLC_Star( Fluid * L, Fluid * EST, float Ls, float Cs){

  float OverRelVel = 1.0/(Ls - Cs);
  float G = (Ls - L->vx)*OverRelVel;

  bx  = EST->bx;
  by  = EST->by;
  bz  = EST->bz;
  
  rho = L->rho * G;
  mx  = rho * Cs;
  my  = L->my * G - (bx*by - L->bx*L->by)*OverRelVel;
  mz  = L->mz * G - (bx*bz - L->bx*L->bz)*OverRelVel;
  vx  = Cs;
  vy  = my/rho;
  vz  = mz/rho;

  press_tot = L->rho*(Ls - L->vx)*(Cs - L->vx) + L->press_tot - L->bx*L->bx + bx*bx;

  if( EquationOfState == 0 )
    Eng = L->Eng*G + ( (press_tot*Cs - L->press_tot*L->vx) -
		       (bx*(EST->bx*EST->vx+ EST->by*EST->vy + EST->bz * EST->vz )
			- L->bx*(L->bx*L->vx + L->by*L->vy + L->bz*L->vz ) ) )*OverRelVel;

  //Fill the vector of conserved quantites.
  cons[ Sden  ] = rho;
  cons[ Sv[0] ] = mx;
  cons[ Sv[1] ] = my;
  cons[ Sv[2] ] = mz;
  if( EquationOfState == 0 ) cons[ Seng ]  = Eng;
  cons[ Sb[0] ] = by;
  cons[ Sb[1] ] = bz;
  
  return SUCCESS;
}

#endif //ATHENA


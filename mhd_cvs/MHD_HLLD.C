
//
// HLLD Solver, coded by dcollins.  
// from Miyoshit & Kusano, JCP 208, 315-344,  2005.
// Coded seperately from Ustugov in order to find errors.
// States chosen in order to mimic the paper, with 'x' denoting '*'
//

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

int MHD_HLLD(float * Fluxes, Fluid * L, Fluid * R ){

  //<dbg>
  //fprintf(stderr," HLLD solver entered  \n");
  //</dbg>

  int field, wave;
  int dbg =0; 
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
  float SL, SLx, SM, SRx, SR;       //Outer (Left and Right) and Center (star) signal speed
  float PTx;                        //central pressure
  float OverSs;                     //Temp variable: 1/( RhoR(SR-VxR) - RhoL(SL-VxL))

  //Intermediate States
  //Note that this isn't the most well thought out object--Only the needed elements are computed,
  //and there's no check on what's been done.  
  Fluid ULx, URx, ULxx, URxx;

  //
  // Compute the Roe Averages, for sound speed.
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
  //vars I need:
  Rt = Ra2 + (R->bx*R->bx + R->by*R->by + R->bz*R->bz)/R->rho;
  Lt = La2 + (L->bx*L->bx + L->by*L->by + L->bz*L->bz)/L->rho;

  Rcf = sqrt( 0.5*( Rt + sqrt( Rt*Rt - 4*Ra2*R->bx*R->bx/R->rho) ) );
  Lcf = sqrt( 0.5*( Lt + sqrt( Lt*Lt - 4*La2*L->bx*L->bx/L->rho) ) );

  ///note: several options exist for this.  Here's one:
  //SR = max( Evals[ NumberOfMHDFluxes - 1] , R->vx + Rcf );
  //SL = min( Evals[0] , L->vx - Lcf );
  //And here's another.
  SL = min(R->vx,L->vx) - max( Rcf, Lcf);
  SR = max(R->vx,L->vx) + max( Rcf, Lcf);

  if( SL > 0  ){
    for( field=0; field<NumberOfMHDFluxes; field++)
      Fluxes[field] = L->Flux[field];

  } else if( SR < 0 ){
    for( field=0; field<NumberOfMHDFluxes; field++)
      Fluxes[field] = R->Flux[field];

  } else {
    //Flux is inside the fan, so more work must be done.
    
    OverSs = 1./((SR - R->vx )*R->rho - (SL - L->vx)*L->rho );
    SM =(  ( (SR-R->vx)*R->rho*R->vx  - (SL-L->vx)*L->rho*L->vx ) - R->press_tot + L->press_tot )*
      OverSs;
    
    PTx =( (SR - R->vx)*R->rho*L->press_tot - (SL - L->vx)*L->rho*R->press_tot +
	   L->rho*R->rho*(SR - R->vx)*( SL - L->vx)*(R->vx - L->vx) )*OverSs;

    ULx.HLLD_x( L, SL, PTx, SM );
    URx.HLLD_x( R, SR, PTx, SM );
    SLx = SM - fabs( L->bx ) / sqrt( ULx.rho );
    SRx = SM + fabs( R->bx ) / sqrt( URx.rho );
    
    if( SLx > 0 ){
      for( field=0; field<NumberOfMHDFluxes; field++)
	Fluxes[field] =  L->Flux[field] + SL*ULx.cons[field] - SL*L->cons[field];

    }else if ( SRx < 0 ){
      for( field=0; field<NumberOfMHDFluxes; field++)
	Fluxes[field] =  R->Flux[field] + SR*URx.cons[field] - SR*R->cons[field];
    } else { //SLx < 0, SRx > 0
      
      if( SM > 0 ){
	ULxx.HLLD_xx( &ULx, &URx, 0 );
	for( field=0; field<NumberOfMHDFluxes; field++)
	  Fluxes[field] = L->Flux[field] + SLx*ULxx.cons[field] - (SLx - SL )*ULx.cons[field] - SL*L->cons[field];

      }else{ //SM < 0
	URxx.HLLD_xx( &ULx, &URx, 1 );
	for( field=0; field<NumberOfMHDFluxes; field++)
	  Fluxes[field] = R->Flux[field] + SRx*URxx.cons[field] - (SRx - SR )*URx.cons[field] - SR*R->cons[field];
      } //SM < 0

    }//SLx < 0, SRx > 0

  } // SL < 0, SR > 0
  
  if( dbg == 1 ) fprintf(stderr,"\n");
  return SUCCESS;

}

int Fluid::HLLD_xx( Fluid * SLx, Fluid * SRx, int Side){
  //The continuous variables: This isn't completely symmetric.

  float sqrtRhoL = sqrt(SLx->rho), sqrtRhoR = sqrt(SRx->rho);
  float sqrts = 1./(sqrtRhoL + sqrtRhoR);
  int sgn = ( (SLx->bx > 0 ) ? 1 : -1 );  //(note: if bx = 0, this code isn't called.

  //Variables Continuous throughout the fan.
  press_tot = SLx->press_tot;
  vx = SLx->vx;
  bx = SLx->bx;
  
  vy = ( sqrtRhoL * SLx->vy + sqrtRhoR * SRx->vy + (SRx->by - SLx->by) * sgn ) * sqrts;
  vz = ( sqrtRhoL * SLx->vz + sqrtRhoR * SRx->vz + (SRx->bz - SLx->bz) * sgn ) * sqrts;

  by = ( sqrtRhoL * SRx->by + sqrtRhoR * SLx->by + sqrtRhoL*sqrtRhoR*(SRx->vy - SLx->vy) * sgn ) * sqrts;
  bz = ( sqrtRhoL * SRx->bz + sqrtRhoR * SLx->bz + sqrtRhoL*sqrtRhoR*(SRx->vz - SLx->vz) * sgn ) * sqrts;

  if( Side == 0 ){
    rho = SLx->rho;
    Eng = SLx->Eng - sgn*sqrtRhoL*( SLx->vx*SLx->bx + SLx->vy * SLx->by + SLx->vz* SLx->bz
				-(vx*bx + vy*by + vz*bz ) );
  }else{ //Side = 1
    rho = SRx->rho;
    Eng = SRx->Eng + sgn*sqrtRhoR*( SRx->vx*SRx->bx + SRx->vy * SRx->by + SRx->vz* SRx->bz
				-(vx*bx + vy*by + vz*bz ) );
    
  }// side

  mx = vx * rho;
  my = vy * rho;
  mz = vz * rho;
  cons[ Sden  ] = rho;
  cons[ Sv[0] ] = mx;
  cons[ Sv[1] ] = my;
  cons[ Sv[2] ] = mz;
  if( EquationOfState == 0 ) cons[ Seng ]  = Eng;
  cons[ Sb[0] ] = by;
  cons[ Sb[1] ] = bz;
  
  return FAIL;
}

//Fill the first intermediate state, ULx or URx.  These are obtained by
//assuming central longitudinal velocity is constant, and integrating across the jumps.
//Written as the Left state, but the Right state is identical.
int Fluid::HLLD_x( Fluid * L, float SL, float PTx, float SM){

  rho = L->rho* ( SL - L->vx )/( SL - SM );
  float denom = 1./( L->rho*(SL - L->vx)*(SL - SM ) - L->bx * L->bx );
  vx = SM;
  vy = L->vy - L->bx*L->by * (SM - L->vx )*denom;
  vz = L->vz - L->bx*L->bz * (SM - L->vx )*denom;

  bx = L->bx;
  by = L->by * (L->rho*(SL - L->vx)*(SL - L->vx) - L->bx*L->bx )*denom;
  bz = L->bz * (L->rho*(SL - L->vx)*(SL - L->vx) - L->bx*L->bx )*denom;

  Eng =(  (SL - L->vx)*L->Eng - L->press_tot * L->vx + PTx*SM 
	  + L->bx*(L->vx*L->bx + L->vy*L->by + L->vz*L->bz
		 -(vx*bx + vy*by + vz*bz) ) )/(SL - SM);
  press_tot = PTx;
  mx = vx * rho;
  my = vy * rho;
  mz = vz * rho;
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


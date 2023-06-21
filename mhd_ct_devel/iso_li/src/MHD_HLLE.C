
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

int MHD_HLLE(float * Fluxes, Fluid * L, Fluid * R ){

  int field, wave;
  
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

  //various signal speeds: eigen vectors,
  //left & right (fast mhd speed, sound speeds, Temp variable, signal speed for HLLE approx.)
  float Speeds[MAX_MHD_WAVES];
  float Lcf, Rcf, La2, Ra2, Lt, Rt, Ls, Rs;
  float OverSs;

  for(field=0;field<3;field++)
    vBar[field] = (sqrtRhoL * L->vel[field]  + sqrtRhoR * R->vel[field] )*overSqrts;

  // Not sure what to do with this mother fucker once HLL average of Bx comes in.
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
    esys_roe_adb_mhd(rhoBar, vBar[0], vBar[1], vBar[2], EnthBar, BxBar, ByBar,BzBar, XX, YY, Speeds, NULL,NULL);
  }else if (EquationOfState == 1 ){
    esys_roe_iso_mhd(rhoBar, vBar[0], vBar[1], vBar[2], IsothermalSoundSpeed*IsothermalSoundSpeed,
		     BxBar, ByBar,BzBar, XX, YY, Speeds, NULL,NULL);
  }

  //copute Sl, Sr
  //switch based on its sign
  //compute.

  
  if( EquationOfState == 0 ){
    Ra2 = Gamma * R->pressure / R->rho;
    La2 = Gamma * L->pressure / L->rho;
  }else if (EquationOfState == 1 ){
    Ra2 = IsothermalSoundSpeed*IsothermalSoundSpeed;
    La2 = IsothermalSoundSpeed*IsothermalSoundSpeed;
  }

  //Compute left and right Fast Speed
  Rt = Ra2 + (R->bx*R->bx + R->by*R->by + R->bz*R->bz)/R->rho;
  Lt = La2 + (L->bx*L->bx + L->by*L->by + L->bz*L->bz)/L->rho;

  Rcf = sqrt( 0.5*( Rt + sqrt( Rt*Rt - 4*Ra2*R->bx*R->bx/R->rho) ) );
  Lcf = sqrt( 0.5*( Lt + sqrt( Lt*Lt - 4*La2*L->bx*L->bx/L->rho) ) );

  Rs = max( Speeds[ NumberOfMHDFluxes - 1] , R->vx + Rcf );
  Ls = min( Speeds[0] , L->vx - Lcf );

  OverSs = 1.0/(Rs - Ls);


  if( Ls > 0 ){
    for( field=0; field<NumberOfMHDFluxes; field++){
      Fluxes[field] = L->Flux[field];
    }
  }else if( Rs < 0 ){
    for( field=0; field<NumberOfMHDFluxes; field++){
      Fluxes[field] = R->Flux[field];
    }
  }else{
    for( field=0; field<NumberOfMHDFluxes; field++){
      Fluxes[field] = (Rs*L->Flux[field] - Ls*R->Flux[field] + Ls*Rs*DU[field])*OverSs;
    }
  }

  
  int EjectionSeat = 0;
  for(field=0;field<NumberOfMHDFluxes; field++)
    if( Fluxes[field] != Fluxes[field] ) {
      EjectionSeat = 1;
      fprintf(stderr,"Bad flux in MHD_HLLE solver. Abort. %d\n", field);
    }
  if( EjectionSeat == 1 ){
    if( SuggestFailure == FALSE ){
      for( int j=0; j<10;j++)
	WriteInThisA[j] = 30 + j; 
      fprintf(stderr," HLLE suggesting failue.\n");
      MHD_WriteDivB = FALSE; //I know that divB != 0, we're rejecting from the middle of the solve.
      fprintf(stderr," n Speeds:  Fluxes: FluxesL: FluxesR\n");
      for( wave=0; wave<NumberOfMHDFluxes;wave++)
	fprintf(stderr," %d %f %f %f %f\n", wave, Speeds[wave], Fluxes[wave], L->Flux[wave],R->Flux[wave]);
      
    }//suggestion
    SuggestFailure = TRUE;
  }//eject
  

  return SUCCESS;

}



#endif //ATHENA


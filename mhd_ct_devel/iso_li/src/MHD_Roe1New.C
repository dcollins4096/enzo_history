
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

int MHD_Roe1(float * Fluxes, Fluid * L, Fluid * R ){
  int field, wave;

  //Eigen System.  Note that RhsV and LhsV, and Strength should be removed
  //after I've tested this shit.

  float RhsV[MAX_MHD_WAVES][MAX_MHD_WAVES], LhsV[MAX_MHD_WAVES][MAX_MHD_WAVES];
  float Speeds[MAX_MHD_WAVES], Strength[MAX_MHD_WAVES];

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

  if( AllCellCentered == 1){
    BxBar = (sqrtRhoL * R->bx + sqrtRhoR * L->bx)*overSqrts;
  }else{
    BxBar = L->bx;  //Or R.bx, they're the same for face centered runs.
  }

  ByBar = (sqrtRhoL * R->by + sqrtRhoL * L->by )*overSqrts;
  BzBar = (sqrtRhoL * R->bz + sqrtRhoL * L->bz )*overSqrts;

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
    esys_roe_adb_mhd(rhoBar, vBar[0], vBar[1], vBar[2], IsothermalSoundSpeed*IsothermalSoundSpeed, BxBar, ByBar,BzBar, XX, YY, Speeds, RhsV, LhsV);
    //esys_roe_iso_mhd(rho, vx, vy, vz, IsothermalSoundSpeed, bx, by, bz, XX, YY, Speeds, RhsV, LhsV);
  }


  for( wave =0; wave < NumberOfMHDFluxes; wave++){

    Strength[wave] = 0.0;

    for( field=0; field< NumberOfMHDFluxes; field++)
      Strength[wave] += LhsV[wave][field]*DU[field];

  }

  for(field=0;field<NumberOfMHDFluxes; field++){
    Fluxes[field] = 0.5*( L->Flux[field] + R->Flux[field] );
    for( wave =0; wave < NumberOfMHDFluxes; wave++)
      Fluxes[field] -= 0.5*RhsV[field][wave]*fabs(Speeds[wave])*Strength[wave];
  }

  return SUCCESS;

}

#endif //ATHENA



#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "MHD_Athena.h"

#include "MHDFluid.h"
#ifdef ATHENA

extern "C" void FORTRAN_NAME(hllds)(float *ql, float *qr, float *flux, float * gamma, int * test_flag);
extern "C" void FORTRAN_NAME(hlld_iso)(float *ql, float *qr, float *flux, float * Interface,
				       float * sound_speed, int * test_flag);
int MHD_HLLD_Sergey(float * Fluxes, Fluid * L, Fluid * R ){


  if( dccCounter9 == -2 ){
    fprintf(stderr,"car\n");
    IsothermalSoundSpeed = -1;
  }

  //In the isothermal case, the energy term is not filled with anything useful.
  float Right[8] = { R->rho, R->mx, R->my, R->mz, R->bx, R->by, R->bz, R->Eng};
  float Left[8] = { L->rho, L->mx, L->my, L->mz, L->bx, L->by, L->bz, L->Eng};
  float Interface[8] = {0,0,0,0,0,0,0,0};  //Not used in Athena, only in PPML.
  float SolverFluxes[8];

  if( EquationOfState == 0 ){
    FORTRAN_NAME(hllds)(Left, Right, SolverFluxes, &Gamma, &dccCounter9);
  }else{
    FORTRAN_NAME(hlld_iso)(Left,Right,SolverFluxes,Interface,&IsothermalSoundSpeed,&dccCounter9);
  }
  Fluxes[ Sden ]  = SolverFluxes[0]; // density
  Fluxes[ Sv[0] ] = SolverFluxes[1]; // vx
  Fluxes[ Sv[1] ] = SolverFluxes[2]; // vy
  Fluxes[ Sv[2] ] = SolverFluxes[3]; // vz
  //No bx flux.
  Fluxes[ Sb[0] ] = SolverFluxes[5]; // by (transverse 1)
  Fluxes[ Sb[1] ] = SolverFluxes[6]; // bz (transverse 2)
  if( EquationOfState == 0 )
    Fluxes[ Seng  ] = SolverFluxes[7]; // Energy



  if( dccCounter9 == -1 ){
    for(int field = 0; field < NumberOfMHDFluxes; field++)
      Fluxes[ field ] = 9;//L->cons[field];
    fprintf(stderr,"hack flux 9!\n");
  }

  return SUCCESS;
}
#endif //ATHENA


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
int MHD_TestFlux(float * Fluxes, Fluid * L, Fluid * R ){

  //In the isothermal case, the energy term is not filled with anything useful.
  float Right[8] = { R->rho, R->mx, R->my, R->mz, R->bx, R->by, R->bz, R->Eng};
  float Left[8] = { L->rho, L->mx, L->my, L->mz, L->bx, L->by, L->bz, L->Eng};
  float Interface[8] = {0,0,0,0,0,0,0,0};  //Not used in Athena, only in PPML.
  float SolverFluxes[8];

  Fluxes[ Sden ]  = R->rho;
  Fluxes[ Sv[0] ] = R->mx;
  Fluxes[ Sv[1] ] = R->my;
  Fluxes[ Sv[2] ] = R->mz;
  //No bx flux.
  Fluxes[ Sb[0] ] = R->by;
  Fluxes[ Sb[1] ] = R->bz;
  if( EquationOfState == 0 )
    Fluxes[ Seng  ] = R->Eng;


  return SUCCESS;
}
#endif //ATHENA

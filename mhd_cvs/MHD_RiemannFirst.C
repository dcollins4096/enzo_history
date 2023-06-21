
// The dumbest thing you can possibly do in CFD.
// This is, in fact, a bad numerical method.  Used by itself,
// it's dramatically prone to even-odd oscilaltions.
// Equivalent to U^{n+1} = U^n - A( F( U_{i+1} ) - F( U_{i-1} ) )
// It's used in this routine ONLY as a first order predictor step.

#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "MHD_Athena.h"
#include "MHDFluid.h"
#ifdef ATHENA

int MHD_1dFlux(float * Fluxes, float * Rhs);


int MHD_RiemannFirst(float * Fluxes, Fluid * L, Fluid * R ){
  int field;

  for(field=0;field<NumberOfMHDFluxes; field++)
    Fluxes[field] = 0.5*( L->Flux[field] + R->Flux[field] );

  //See?  Pretty dumb.
  return SUCCESS;
}

#endif

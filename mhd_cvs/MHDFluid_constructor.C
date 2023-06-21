
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "MHDFluid.h"

#ifdef ATHENA

Fluid::Fluid(){

  rho = 0;
  mx = 0;
  my = 0;
  mz = 0;
  by = 0;  //The order of the fluid quantities is very important
  bz = 0;  //The *cons pointer is used as an iterator over conserved quantities,
  Eng = 0; //and the operations it's used in expect this ordering.
  bx = 0;  
  vx = 0;
  vy = 0;
  vz = 0;
  pressure = 0;
  enth = 0;
  press_tot = 0;
  
  momentum = &mx;
  vel = & vx;

  //the vectors cons and Flux are filled in 
  //the routine MHDFluid_Fill
  
}

Fluid::~Fluid(){
  
}

#endif //ATHENA

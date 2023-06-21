
//for multiple inclusions
#ifdef ATHENA

#ifndef MHDFLUID
#define MHDFLUID

#include "MHD_Athena.h"

//avast
class Fluid{
 public:

  //The conserved vector
  float rho;
  float mx;
  float my;
  float mz;
  float Eng;
  float by;
  float bz;
  float bx;

  //The primitive quantities
  float vx;
  float vy;
  float vz;
  float pressure;
  float enth;
  float press_tot;
  
  //Iterators over the vectors.
  float cons[MAX_MHD_WAVES];
  float * vel;
  float * momentum;
  float Flux[MAX_MHD_WAVES+1];

  //Data construction routines.
  Fluid();  //constructor
  ~Fluid(); //Destructor

  //Assign values for all these things.
  int Fill(float * ReconstructedState);   

  //Take Primitive Flux Derivative.  
  int PrimitiveDerivative(float * dFlux, float * dU);

  //Other routines
  int HLLC_Star( Fluid * L, Fluid * HLLE, float Ls, float Cs );
  int HLLD_x( Fluid * L, float SL, float PTx, float SM);
  int HLLD_xx( Fluid * SLx, Fluid * SRx, int Side);

};



#endif //MHDFLUID

#endif //ATHENA

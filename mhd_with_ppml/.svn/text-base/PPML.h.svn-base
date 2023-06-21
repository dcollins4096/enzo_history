
#ifndef PPML_HEADER
#define PPML_HEADER

class grid;
//struct PPML_InterfacePointerBundle;
//struct IndexPointerMap;

class PPML_InterfacePointerBundle{

 public:
  PPML_InterfacePointerBundle();
  PPML_InterfacePointerBundle(grid * Grid);
  ~PPML_InterfacePointerBundle();
  int NumberOfFluidQuantities;
  int PPML_NFaces;
  int Error;
  float ** X_L;
  float ** X_R;
  float ** Y_L;
  float ** Y_R;
  float ** Z_L;
  float ** Z_R;
  float *** All;

};

//Both vector V[3] and scalar VX,VY,VZ are given for code readability.
//This is filled in IdentifyPhysicalQuantities.C, the routine IdentifyPhysicalQuantities_2
struct IndexPointerMap{
  int D,TE, VX, VY, VZ, BX,BY,BZ, GE;
  int V[3], B[3];

#ifdef MHDF
  float * CenteredB[3];
#endif //MHDF
};

#endif

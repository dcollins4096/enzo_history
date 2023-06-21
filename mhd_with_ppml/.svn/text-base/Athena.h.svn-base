
#ifndef ATHENA_HEADDER
#define ATHENA_HEADDER

//
// The davethena method.  
//

// Initial refactor of the Athena method, making it more object oriented.
// 
// Nov. 6 2007: 
//   I only want this as a warpper for the viscosity method,
// to use it with PPML.  Trying to keep the rest of the method in mind.
// Might get a refactor when the rest of the code arrives.
//


//I don't need everything that's in this class: We should restructure this headder.
#include "PPML.h"
class grid;
#define MAX_MHD_WAVES 7  
class Athena{

 public:
  Athena();
  Athena(grid * Grid);
  ~Athena();

  int Error; //Error checking on the constructor.  There's got to be a better way....
  
 private:

  //For loop extents.  7 for Adiabatic MHD, 6 for Isothermal.
  int NumberOfMHDFluxes;

  //Loop extents for reconstruction.  
  int NumberReconstructed;  
  //AllCellCentered=1 averages longitudinal magnetic field.
  //               =0 uses the magnetic field at the interface.
  int AllCellCentered;


  // The set of pointers that Athena will loop over.  
  // This map changes depending on the direction of interest.
  float * State[MAX_NUMBER_OF_BARYON_FIELDS + 3];


  //Map from Physics index to Solver Index.  (the hassle of using someone elses code...)
  //The first index is the dimension of the sweep.
  int MapEtoS[3][MAX_MHD_WAVES];
  int BNum[3][2];

  //replaced with the IndexPointerMap
  ///int S.TE, Sden, Sv[3], Sb[3];
  //int E.TE, Eden, Ev[3], Egas; //Note that Egas isn't actually used.
  IndexPointerMap E, S; //Enzo and State.  I do hate one letter variable names...

  // Meta data that's coppied directly from the Grid object.
  int GridRank;                        // number of dimensions
  int GridDimension[MAX_DIMENSION];    // total dimensions of all grids
  int GridStartIndex[MAX_DIMENSION];   // starting index of the active region
                                       //   (zero based)
  int GridEndIndex[MAX_DIMENSION];     // stoping index of the active region
                                       //   (zero based)

  float dtFixed;                       // current (fixed) timestep
  int    NumberOfBaryonFields;         // active baryon fields
  int   NumberOfFluidQuantities;       // Fluid quantities (non-face-centered stuff.)
  float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
  float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
  float *RandomForcingField[MAX_NUMBER_RANDOM_FORCING]; // pointers to arrays //AK
  FLOAT *CellWidth[MAX_DIMENSION];

#ifdef MHDF
  //Not sure if this is the best way to do things-- May get changed in the future.
  //If dcollins graduated, delete this comment.
  float * CenteredB[MAX_FACE_FIELDS];
  float * MagneticField[MAX_FACE_FIELDS];
  float * ElectricField[MAX_FACE_FIELDS];

  int NFaces, NEdges;
  int MagneticDims[MAX_FACE_FIELDS][MAX_DIMENSION];
  int ElectricDims[MAX_EDGE_FIELDS][MAX_DIMENSION];
#endif //MHDF
  // Offset is for dimension-independant derivative formulation.
  //    so instead of Flux[0][ i+1 ] - Flux[0][ i ], Flux[1][ j+1 ] - Flux[1][ j ] 
  //    I say once Flux[dim][ ijk + Offset[dim] ] - Flux[dim][ ijk ].
  //    Only somewhat less readable, much easier to maintian.
  // dxI[dim] = 1/dx.  Initialized to zero to catch errors.  
  
  int Offset[3];
  float dxI[3];
  
  //The athena ghost zone arrays.  
  int AthStart[3], AthEnd[3], GhostZones[2];

  // Solver parameters. Set and re-set depending on dimension.
  // Loosely match the Enzo global parameters of the similar name.
  int ReconLocal, RiemannLocal;
  int MHD_PLM_SlopeLocal[MAX_NUMBER_OF_BARYON_FIELDS + 3];


  float * FlatteningField;
  int * ShockDirection;

  //<dbg> Debugging flags
  // a counter.
  int ath_counter;
  //For noting failures in Roe's method.
  int CheckForSolverProblems;
  float * SolverProblems[3];
  //</dbg>
  

 public:
  int RotateState(int dim);
  int PPMLViscosity(float * FluxX,float * FluxY, float * FluxZ, int SweepDirection);
  float MHD_Viscosity(float * Fluxes, float * Lhs, float * Rhs, int * index, 
		      int dim, int MHD_DiffusionMethodInput, float dT);
#ifdef MHDF
  int ComputeElectricField(float dT, float * Fluxes[]);
#endif //MHDF
  // Here goes the Outer Loop, the Flux Difference, the CT.  Later.


  
};


#endif //ATHENA_HEADDER



//prevent multiple inclusions
#ifndef AthenaHeader
#define AthenaHeader


#ifdef ATHENA


//
// Header file for Athena MHD suite.
//
//
//1) read with    DEFINE_STORAGE defined for the (single) definition (in Grid_MHD_Athena.C)
//2) read without DEFINE_STORAGE defined for external linkage
//

#ifdef DEFINE_STORAGE
# define ATH_EXTERN
#else /* DEFINE_STORAGE */
# define ATH_EXTERN extern
#endif /* DEFINE_STORAGE */

ATH_EXTERN int NumberOfMHDFluxes;
#define MAX_MHD_WAVES 7
//
ATH_EXTERN int AllCellCentered, NumberReconstructed;

//Map from Physics index to Solver Index.  (the hassle of using
//someone elses code...)
//Note that there are 3 coppies, because there are 3 different maps:  this is a 1d solver, after all.

ATH_EXTERN int MapEtoS[3][MAX_MHD_WAVES];
ATH_EXTERN int BNum[3][2];
ATH_EXTERN int Seng, Sden, Sv[3], Sb[3];
//             energy, density, velocity, MagneticField (Sb)
ATH_EXTERN int Eeng, Eden, Ev[3], Egas; //Note that Egas isn't actually used.
ATH_EXTERN int MHD_PLM_SlopeLocal[MAX_NUMBER_OF_BARYON_FIELDS + 3];
//The pointers that will get looped over in the reconstruction step.
ATH_EXTERN float * State[MAX_NUMBER_OF_BARYON_FIELDS + 3];
ATH_EXTERN float * FlatteningField;
ATH_EXTERN int * ShockDirection;
// Offset is for dimension-independant derivative formulation.
// dxI[dim] = 1/dx.  Initialized to zero to catch errors.  
// Both Filled in Grid_MHD_Athena.

ATH_EXTERN int Offset[3];
ATH_EXTERN float dxI[3];

//The athena ghost zone arrays.  
ATH_EXTERN int AthStart[3], AthEnd[3], GhostZones[2];

//<dbg> a counter.
ATH_EXTERN int ath_counter;
//</dbg>

// For solution sensitive switches.  I'm not sure if these definitely need to be here, but it's easier.
ATH_EXTERN int ReconLocal, RiemannLocal;

//<dbg>
//For noting failures in Roe's method.
ATH_EXTERN int CheckForSolverProblems;
ATH_EXTERN float * SolverProblems[3];
//</dbg>


#endif //ATHENA


#endif //AthenaHeader

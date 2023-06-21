/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Single-Group, Multi-species, Gray Flux-Limited Diffusion Implicit 
/  Problem Class
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: This class defines problem-specific functions for an 
/           implicit gray flux-limited diffusion solve.
/
/           The variables are stored in the following order: 
/              0 -> radiation energy density
/              1 -> fluid energy correction
/              2:Nspecies+1 -> chemical species (Nspecies may be 0)
/
************************************************************************/

#ifdef RAD_HYDRO
#ifndef FLD_IMPLICIT_PROBLEM_DEFINED__
#define FLD_IMPLICIT_PROBLEM_DEFINED__

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EnzoVector.h"
#include "ImplicitProblemABC.h"


class gFLDProblem : public virtual ImplicitProblemABC {

 private:
  
  // flag denoted problem preparedness
  bool prepared;
  
  // HYPRE SStruct-specific data
  Eint32 mattype;                // HYPRE matrix type for solve
  Eint32 stSize;                 // stencil size
  HYPRE_SStructGrid grid;        // HYPRE grid object for setup
  HYPRE_SStructStencil stencil;  // stencil object
  HYPRE_SStructGraph graph;      // graph object

  // HYPRE Solver-specific data
  Eint32 sol_zeroguess;          // use a zero initial guess
  Eint32 sol_maxit;              // maximum number of iterations
  Eint32 sol_relch;              // relative change stopping criteria
  Eint32 sol_rlxtype;            // relaxation type
  Eint32 sol_npre;               // num. pre-relaxation sweeps
  Eint32 sol_npost;              // num. post-relaxation sweeps
  Eint32 sol_printl;             // print output level
  Eint32 sol_log;                // amount of logging
  Eint32 SolvIndices[3][2];      // global L/R edge indices of subdomain 

  // HYPRE solver diagnostics
  int totIters;                  // total MG iterations for solves

  // Inexact Newton solver-specific data
  int newt_maxit;                // maximum number of iterations
  int newt_norm;                 // norm for convergence measurement
  float newt_INconst;            // Inexact-Newton constant
  float newt_tol;                // Newton tolerance
  float newt_MinLinesearch;      // minimum allowed line-search length

  // General problem grid information
  bool OnBdry[3][2]; // denotes if proc owns piece of boundary
  int rank;          // Rank of self-gravity problem
  int layout[3];     // number of procs in each dim (1-based)
  int location[3];   // location of this proc in each dim (0-based)
  int NBors[3][2];   // process IDs of L/R neighbors in each dim
  int LocDims[3];    // implicit problem local dims (no ghost or bdry cells)
  int ArrDims[3];    // local array sizes (includes ghost and bdry cells)
  int GhDims[3][2];  // ghost cells at each face
  int GlobDims[3];   // implicit problem global dimensions (active cells only)
  float dx[3];             // mesh size in each dim
  float EdgeVals[3][2];    // L/R edges of this proc's gravity domain
  float *EBdryVals[3][2];  // boundary values for radiation dirichlet BCs

  // time-stepping related data
  float tnew;        // new time
  float told;        // old time
  float dt;          // time step size
  float theta;       // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  int LimImp;        // implicitness of flux limiter:
                     //    0 -> fully lagged to previous time step
                     //    1 -> fully lagged to previous newton iterate
                     //    2 -> lag only temperature dependence
  EnzoVector *U0;    // old time-level state
  EnzoVector *rhs;   // current time-level rhs
  EnzoVector *rhs0;  // old time-level rhs

  // problem defining data
  float ErUnits;     // scaling coefficient for radiation energy density
  FLOAT anew;        // cosmology expansion coefficient (at tnew)
  FLOAT adotnew;     // time-derivative of a (at tnew)
  FLOAT aold;        // cosmology expansion coefficient (at told)
  FLOAT adotold;     // time-derivative of a (at told)
  int Nchem;         // number of chemical species (non-negative integer)
  int Model;         // model choice, 0=>keep HI variable, but omit coupling to ec, Eg
                     //               1=>case B HII recomb., no emissivity
                     //               2=>case A HII recomb., with emissivity
                     //              10=>fully decoupled, constant opacity

  // storage for integrals over radiation spectrum (set during initialization)
  int ESpectrum;            // integer flag determining spectrum choice
  float intSigE;            // int_{nu0}^{inf} sigma_E(nu) d nu
  float intSigESigHI;       // int_{nu0}^{inf} sigma_E(nu)*sigma_HI(nu) d nu
  float intSigESigHeI;      // int_{nu0}^{inf} sigma_E(nu)*sigma_HeI(nu) d nu
  float intSigESigHeII;     // int_{nu0}^{inf} sigma_E(nu)*sigma_HeII(nu) d nu
  float intSigESigHInu;     // int_{nu0}^{inf} sigma_E(nu)*sigma_HI(nu)/nu d nu
  float intSigESigHeInu;    // int_{nu0}^{inf} sigma_E(nu)*sigma_HeI(nu)/nu d nu
  float intSigESigHeIInu;   // int_{nu0}^{inf} sigma_E(nu)*sigma_HeII(nu)/nu d nu

  // linear solver/Jacobian arrays
  EnzoVector **L;    // local Jacobian components 

  // access to Enzo data
  float *vx;         // x0-directional velocity
  float *vy;         // x1-directional velocity
  float *vz;         // x2-directional velocity
  float *rho;        // density
  float *eh;         // total fluid energy
  float *ne;         // electron number density

  // stored arrays for increased efficiency
  float *Temp;
  float *sigmaA;
  float *sigmaS;
  
  // private computation routines
  int LocRHS(EnzoVector *locrhs, float time, 
	     FLOAT a, FLOAT adot, EnzoVector *u);
  int ComputeRHS(EnzoVector *rhsval, float time, EnzoVector *u);
  int ComputeTemperature(float *Temperature, float time, 
			 FLOAT a, EnzoVector *u);
  int MatrixEntries(double *matentries, float *Eg, float *Eg0, 
		    float *Temperature, float *sigmaA, 
		    float *sigmaS, float *adjvec);
  int DiffRHS(float *drhs, float *Eg, float *Eg0, float *Temperature, 
	      float *sigmaA, float *sigmaS, FLOAT *a);
  int LocEcRHS(float *ecrhs, float *ec, float *Eg, float *Temperature,
	       float *nHI, float *nHeI, float *nHeII, FLOAT *a, 
	       FLOAT *adot, float *aUnits, float *DensityUnits,
	       float *TimeUnits, float *LenUnits);
  int LocEgRHS(float *Egrhs, float *Eg, float *Temperature, 
	       float *n_HI, float *Kappa, FLOAT *a, FLOAT *adot);
  int LocNiRHS(float *rhs_HI, float *rhs_HeI, float *rhs_HeII, 
	       float *n_HI, float *n_HeI, float *n_HeII, float *Eg, 
	       float *Temperature, FLOAT *a, FLOAT *adot);
  int BlockSolve(float *Amat, float *xvec, float *bvec, int *N, int *M);
  int Opacity(float *Kappa, float *n_HI, float *n_HeI, float *n_HeII, FLOAT *a);
  float RadiationSpectrum(float nu);
  float CrossSections(float nu, int species);
  int ComputeRadiationIntegrals();

 public:

  // boundary type in each dimension, face (0->periodic, 1->dirichlet, 2->neumann)
  int BdryType[3][2];
  // file name for boundary condition input
  char *BoundaryFName;

  ///////////////////////////////////////
  // FLD-Specific Routines

  // Constructor
  gFLDProblem();
  
  // Destructor
  ~gFLDProblem();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int BdryConst, float *BdryData);

  // Write boundary conditions to file
  int WriteBoundary(FILE *fptr, char *hdfname);

  // Read boundary conditions from file
  int ReadBoundary(FILE *fptr);

  // Enforce boundary conditions onto a vector
  int EnforceBoundary(EnzoVector *vec, int flag);

  // Problem setup
  int Setup(HierarchyEntry &TopGrid, float Dt, float time,
	    float *Er, float *ec, float **ni);

  // Return number of chemical species
  int GetNumChemicalSpecies();

  // Return number of allowed Newton iterations
  int GetMaxNewtonIters();

  // Return choice of norm for Newton convergence
  int GetNewtonNormChoice();

  // Return inexact Newton constant
  float GetInexactNewtonConstant();

  // Return Newton nonlinear tolerance
  float GetNewtonTolerance();

  // Return minimum allowed linesearch length
  float GetNewtonMinLinesearch();

  // Set choice of initial guess
  int SetZeroMGGuess(int zeroguess);

  // Set maximum number of MG iterations
  int SetMaxMGIters(int maxit);

  // Set stopping criteria based on relative residual change
  int SetRelativeChangeStop(int relch);

  // Set MG relaxation type
  int SetMGRelaxationType(int rlxtype);

  // Set number of pre-relaxation sweeps
  int SetNumPreRelaxSweeps(int npre);

  // Set number of post-relaxation sweeps
  int SetNumPostRelaxSweeps(int npost);

  // Set HYPRE diagnostic print level
  int SetPrintLevel(int printl);

  // Set HYPRE logging level
  int SetLogLevel(int log);

  // Set maximum Newton iterations
  int SetNewtonMaxIterations(int maxit);

  // Set Newton norm for convergence
  int SetNewtonNorm(int norm);

  // Set Newton tolerance
  int SetNewtonTolerance(float tol);

  // Set Newton inexactness constant
  int SetInexactNewtonConstant(float INconst);

  // Set minimum Newton Linesearch length
  int SetNewtonMinLinesearch(float MinLength);

  // Problem return
  int Return(EnzoVector *sol, float *Er, float *ec, float **ni);

  // Return example self-gravity vector
  EnzoVector *ExampleVector();
  

  ///////////////////////////////////////
  // Implicit Solver Interface Routines
  
  // Problem-defining nonlinear residual operations (called repeatedly)
  int nlresid(EnzoVector *fu, EnzoVector *u);
  
  // Problem-specific Linear system setup function, sets up the 
  //   linear Newton system matrix J(u) given an updated state u
  //   (called once per Newton iteration)
  int lsetup(EnzoVector *u);
  
  // Problem-specific Linear solver function 
  //   solves J(u)*s = b to tolerance delta
  //   (called once per Newton iteration)
  int lsolve(EnzoVector *s, EnzoVector *b, EnzoVector *u, float delta);
  
};


#endif
#endif

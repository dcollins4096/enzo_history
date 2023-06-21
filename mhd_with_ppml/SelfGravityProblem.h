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
/  Self-Gravity Problem Class
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: This class defines problem-specific functions for an 
/           implicit self-gravity solve.
/
************************************************************************/

#ifdef ISO_GRAV
#ifndef SELF_GRAVITY_PROBLEM_DEFINED__
#define SELF_GRAVITY_PROBLEM_DEFINED__

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "ImplicitProblemABC.h"


class SelfGravityProblem : public virtual ImplicitProblemABC {

 private:

  // flag denoting problem preparedness
  bool prepared;

  // HYPRE SStruct-specific data
  Eint32 mattype;                // HYPRE matrix type for solve
  Eint32 stSize;                 // stencil size
  HYPRE_SStructGrid grid;        // HYPRE grid object for setup
  HYPRE_SStructStencil stencil;  // stencil object
  HYPRE_SStructGraph graph;      // graph object
  HYPRE_SStructMatrix A;         // system matrix

  // HYPRE Solver-specific data
  Eint32 sol_zeroguess;          // use a zero initial guess
  Eint32 sol_maxit;              // maximum number of iterations
  Eint32 sol_relch;              // relative change stopping criteria
  Eint32 sol_rlxtype;            // relaxation type
  Eint32 sol_npre;               // num. pre-relaxation sweeps
  Eint32 sol_npost;              // num. post-relaxation sweeps
  Eint32 sol_printl;             // print output level
  Eint32 sol_log;                // amount of logging
  Eint32 SolvIndices[3][2];      // global L/R edge indices of grav. domain 
                                 // (including boundary cells)

  // HYPRE solver diagnostics
  int AInit;                     // flag denoting initialization of A matrix
  int totIters;                  // total MG iterations for solves


  // General problem information
  bool OnBdry[3][2]; // denotes if proc is on boundary

  int rank;          // Rank of self-gravity problem
  int layout[3];     // number of procs in each dim (1-based)
  int location[3];   // location of this proc in each dim (0-based)
  int NBors[3][2];   // process IDs of L/R neighbors in each dim
  int LocDims[3];    // implicit problem local dims (no ghost or bdry cells)
  int ArrDims[3];    // local array sizes (includes ghost and bdry cells)
  int GlobDims[3];   // implicit problem global dimensions (active cells only)
  int GravGhosts;    // Enzo ghost zones in each dimension
  int BdryType[3];   // boundary type in each dimension (0->per, 1->dir)
  int EdgeIndices[3][2];  // global L/R edge indices of grav. domain
                          // (interior only, no boundary cells)

  float dx[3];            // mesh size in each dim
  float GravEdges[3][2];  // L/R edges of this proc's gravity domain
  float *BdryVals[3][2];  // boundary values for dirichlet BCs

  EnzoVector *rhorhs;     // Self-Gravity problem rhs

 public:

  ///////////////////////////////////////
  // Self-Gravity Specific Routines

  // Constructor (called once per simulation)
  SelfGravityProblem(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Destructor (called once per simulation)
  ~SelfGravityProblem();

  // Problem Boundary Condition setup (called either once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int BdryConst, float *BdryData);

  // Write boundary conditions to file
  int WriteBoundary(FILE *fptr, char *hdfname);

  // Read boundary conditions from file
  int ReadBoundary(FILE *fptr);

  // Enforce boundary conditions onto a vector
  int EnforceBoundary(EnzoVector *vec, int flag);

  // Problem setup (called each time step, SetupBoundary must be called first)
  int Setup(HierarchyEntry &TopGrid);

  // Set choice of initial guess
  int SetZeroGuess(int zeroguess) {
    sol_zeroguess = zeroguess;
    return SUCCESS;
  };

  // Set maximum number of MG iterations
  int SetMaxMGIters(int maxit) {
    sol_maxit = maxit;
    return SUCCESS;
  };

  // Set stopping criteria based on relative residual change
  int SetRelativeChangeStop(int relch) {
    sol_relch = relch;
    return SUCCESS;
  };

  // Set MG relaxation type
  int SetMGRelaxationType(int rlxtype) {
    sol_rlxtype = rlxtype;
    return SUCCESS;
  };

  // Set number of pre-relaxation sweeps
  int SetNumPreRelaxSweeps(int npre) {
    sol_npre = npre;
    return SUCCESS;
  };

  // Set number of post-relaxation sweeps
  int SetNumPostRelaxSweeps(int npost) {
    sol_npost = npost;
    return SUCCESS;
  };

  // Set HYPRE diagnostic print level
  int SetPrintLevel(int printl) {
    sol_printl = printl;
    return SUCCESS;
  };

  // Set HYPRE logging level
  int SetLogLevel(int log) {
    sol_log = log;
    return SUCCESS;
  };

  // Problem return (called following each time step)
  int Return(EnzoVector *sol, float *GravPotential);
  
  // Return example self-gravity vector
  EnzoVector *ExampleVector() {
    return rhorhs;
  };
  

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

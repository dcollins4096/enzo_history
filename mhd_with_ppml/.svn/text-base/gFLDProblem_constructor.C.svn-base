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
/  Gray Flux-Limited Diffusion Implicit Problem Class constructor 
/  routine
/
/  written by: Daniel Reynolds
/  date:       September, 2006
/  modified1:  
/
/  PURPOSE: Initializes all values to illegal numbers, and sets all 
/           arrays to NULL;  Requires call to Initialize to actually 
/           set up these values.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"


gFLDProblem::gFLDProblem()
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::constructor routine\n");

  // initialize prepared flag to false
  prepared = false;

  // initialize HYPRE values to -1/NULL
  mattype = -1;
  stSize = -1;
  grid = NULL;
  stencil = NULL;
  graph = NULL;
  sol_zeroguess = -1;
  sol_maxit = -1;
  sol_relch = -1;
  sol_rlxtype = -1;
  sol_npre = -1;
  sol_npost = -1;
  sol_printl = -1;
  sol_log = -1;
  SolvIndices[0][0] = -1;
  SolvIndices[0][1] = -1;
  SolvIndices[1][0] = -1;
  SolvIndices[1][1] = -1;
  SolvIndices[2][0] = -1;
  SolvIndices[2][1] = -1;
  totIters = -1;

  // initialize Newton solver values to -1/NULL
  newt_maxit = -1;
  newt_norm = -1;
  newt_INconst = -1.0;
  newt_tol = -1.0;
  newt_MinLinesearch = -1.0;

  // initialize problem grid information to -1/NULL
  OnBdry[0][0] = false;
  OnBdry[0][1] = false;
  OnBdry[1][0] = false;
  OnBdry[1][1] = false;
  OnBdry[2][0] = false;
  OnBdry[2][1] = false;
  rank = -1;
  layout[0] = -1;
  layout[1] = -1;
  layout[2] = -1;
  location[0] = -1;
  location[1] = -1;
  location[2] = -1;
  NBors[0][0] = -1;
  NBors[0][1] = -1;
  NBors[1][0] = -1;
  NBors[1][1] = -1;
  NBors[2][0] = -1;
  NBors[2][1] = -1;
  LocDims[0] = -1;
  LocDims[1] = -1;
  LocDims[2] = -1;
  ArrDims[0] = -1;
  ArrDims[1] = -1;
  ArrDims[2] = -1;
  GhDims[0][0] = -1;
  GhDims[0][1] = -1;
  GhDims[1][0] = -1;
  GhDims[1][1] = -1;
  GhDims[2][0] = -1;
  GhDims[2][1] = -1;
  BdryType[0][0] = -1;
  BdryType[0][1] = -1;
  BdryType[1][0] = -1;
  BdryType[1][1] = -1;
  BdryType[2][0] = -1;
  BdryType[2][1] = -1;
  dx[0] = -1.0;
  dx[1] = -1.0;
  dx[2] = -1.0;
  EdgeVals[0][0] = -1.0;
  EdgeVals[0][1] = -1.0;
  EdgeVals[1][0] = -1.0;
  EdgeVals[1][1] = -1.0;
  EdgeVals[2][0] = -1.0;
  EdgeVals[2][1] = -1.0;
  EBdryVals[0][0] = NULL;
  EBdryVals[0][1] = NULL;
  EBdryVals[1][0] = NULL;
  EBdryVals[1][1] = NULL;
  EBdryVals[2][0] = NULL;
  EBdryVals[2][1] = NULL;
  
  // initialize time-stepping related data to -1/NULL
  tnew = -1.0;
  told = -1.0;
  dt = -1.0;
  theta = -1.0;
  LimImp = -1;
  U0 = NULL;
  rhs = NULL;
  rhs0 = NULL;

  // initialize problem defining data 
  ErUnits = 1.0;
  anew = -1.0;
  adotnew = -1.0;
  aold = -1.0;
  adotold = -1.0;
  Nchem = -1;
  Model = -1;
  ESpectrum = -1;
  intSigE = 0.0;
  intSigESigHI = 0.0;
  intSigESigHeI = 0.0;
  intSigESigHeII = 0.0;
  intSigESigHInu = 0.0;
  intSigESigHeInu = 0.0;
  intSigESigHeIInu = 0.0;

  // initialize linear solver/Jacobian arrays to NULL
  L = NULL;

  // initialize access to Enzo arrays to NULL
  vx = NULL;
  vy = NULL;
  vz = NULL;
  rho = NULL;
  eh = NULL;
  ne = NULL;

  // initialize storage arrays to NULL
  Temp = NULL;
  sigmaA = NULL;
  sigmaS = NULL;

  // initialize file name for boundary condition input to NULL
  BoundaryFName = NULL;

}
#endif

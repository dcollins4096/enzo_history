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

#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"


// Return number of chemical species
int gFLDProblem::GetNumChemicalSpecies() {
  return Nchem;
}

// Return number of allowed Newton iterations
int gFLDProblem::GetMaxNewtonIters() {
  return newt_maxit;
}

// Return choice of norm for Newton convergence
int gFLDProblem::GetNewtonNormChoice() {
  return newt_norm;
}

// Return inexact Newton constant
float gFLDProblem::GetInexactNewtonConstant() {
  return newt_INconst;
}

// Return Newton nonlinear tolerance
float gFLDProblem::GetNewtonTolerance() {
  return newt_tol;
}

// Return minimum allowed linesearch length
float gFLDProblem::GetNewtonMinLinesearch() {
  return newt_MinLinesearch;
}

// Set choice of initial guess
int gFLDProblem::SetZeroMGGuess(int zeroguess) {
  sol_zeroguess = zeroguess;
  return SUCCESS;
};

// Set maximum number of MG iterations
int gFLDProblem::SetMaxMGIters(int maxit) {
  sol_maxit = maxit;
  return SUCCESS;
};

// Set stopping criteria based on relative residual change
int gFLDProblem::SetRelativeChangeStop(int relch) {
  sol_relch = relch;
  return SUCCESS;
};

// Set MG relaxation type
int gFLDProblem::SetMGRelaxationType(int rlxtype) {
  sol_rlxtype = rlxtype;
  return SUCCESS;
};

// Set number of pre-relaxation sweeps
int gFLDProblem::SetNumPreRelaxSweeps(int npre) {
  sol_npre = npre;
  return SUCCESS;
};

// Set number of post-relaxation sweeps
int gFLDProblem::SetNumPostRelaxSweeps(int npost) {
  sol_npost = npost;
  return SUCCESS;
};

// Set HYPRE diagnostic print level
int gFLDProblem::SetPrintLevel(int printl) {
  sol_printl = printl;
  return SUCCESS;
};

// Set HYPRE logging level
int gFLDProblem::SetLogLevel(int log) {
  sol_log = log;
  return SUCCESS;
};

// Set maximum Newton iterations
int gFLDProblem::SetNewtonMaxIterations(int maxit) {
  newt_maxit = maxit;
  return SUCCESS;
};

// Set Newton norm for convergence
int gFLDProblem::SetNewtonNorm(int norm) {
  newt_norm = norm;
  return SUCCESS;
};

// Set Newton tolerance
int gFLDProblem::SetNewtonTolerance(float tol) {
  newt_tol = tol;
  return SUCCESS;
};

// Set Newton inexactness constant
int gFLDProblem::SetInexactNewtonConstant(float INconst) {
  newt_INconst = INconst;
  return SUCCESS;
};

// Set minimum Newton Linesearch length
int gFLDProblem::SetNewtonMinLinesearch(float MinLength) {
  newt_MinLinesearch = MinLength;
  return SUCCESS;
};

// Return example self-gravity vector
EnzoVector* gFLDProblem::ExampleVector() {
  return U0;
};
 
#endif

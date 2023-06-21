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
/  Gray Flux-Limited Diffusion Implicit Problem nonlinear residual 
/  function
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Nonlinear residual function that defines the coupled,
/           implicit-time, radiation diffusion/chemistry/fluid energy 
/           system.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"



int gFLDProblem::nlresid(EnzoVector *fu, EnzoVector *u)
{
  
//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::nlresid routine\n");

  // check that the gFLDProblem has been set up
  if (!prepared) {
    fprintf(stderr,"nlresid error: gFLDProblem not yet prepared\n");
    return FAIL;
  }

//   float dtmp = u->l1norm();
//   if (debug)
//     fprintf(stdout,"nlresid chkpnt 0: ||u||_1 = %g\n",dtmp);

  // have u communicate neighbor information
  if (u->exchange() == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::exchange failure\n");
    return FAIL;
  }

  // enforce boundary conditions on state u
  if (this->EnforceBoundary(u,0) == FAIL) {
    fprintf(stderr,"nlresid error: EnforceBoundary failure\n");
    return FAIL;
  }

  // initialize residual to zero
  if (fu->constant(0.0) == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::constant failure\n");
    return FAIL;
  }
  
//   dtmp = fu->l1norm();
//   if (debug)
//     fprintf(stdout,"nlresid chkpnt 1: ||fu||_1 = %g\n",dtmp);

  // compute rhs at current state for updated time
  if (this->ComputeRHS(rhs, tnew, u) == FAIL) {
    fprintf(stderr,"nlresid error: ComputeRHS failure\n");
    return FAIL;
  }

//   dtmp = rhs->l1norm();
//   if (debug)
//     fprintf(stdout,"nlresid chkpnt 1.5: ||rhs||_1 = %g\n",dtmp);

  // combine together via theta-scheme
  //   fu =  u-u0
  if (fu->linearsum(1.0,u,-1.0,U0) == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::linearsum failure\n");
    return FAIL;
  }

//   dtmp = fu->l1norm();
//   if (debug)
//     fprintf(stdout,"nlresid chkpnt 2: ||fu||_1 = %g\n",dtmp);

  //   fu = (u-u0) + dt*(1-theta)*rhs0
  if (fu->axpy(dt*(1.0-theta),rhs0)  == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::axpy failure\n");
    return FAIL;
  }

//   dtmp = fu->l1norm();
//   if (debug)
//     fprintf(stdout,"nlresid chkpnt 3: ||fu||_1 = %g\n",dtmp);

  //   fu = (u-u0) + dt*(1-theta)*rhs0 + dt*theta*rhs
  if (fu->axpy(dt*theta,rhs) == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::axpy failure\n");
    return FAIL;
  }

  //   Enforce boundary conditions on fu (for Dirichlet faces)
  if (this->EnforceBoundary(fu,1) == FAIL) {
    fprintf(stderr,"nlresid error: EnforceBoundary failure\n");
    return FAIL;    
  } 

//   dtmp = fu->l1norm();
//   if (debug)
//     fprintf(stdout,"nlresid chkpnt 4: ||fu||_1 = %g\n",dtmp);

//   dtmp = u->l1norm();
//   if (debug)
//     fprintf(stdout,"nlresid chkpnt 5: ||u||_1 = %g\n",dtmp);

  // return success
  return SUCCESS;

}
#endif

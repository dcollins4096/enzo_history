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
/  Gray Flux-Limited Diffusion Implicit Problem linear Newton system 
/  setup function
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Called by implicit solver to notify the Problem of 
/           updates to the current state (given in the vector u), 
/           so that the linear system matrix J(u) may be updated 
/           if necessary.  For the gray FLD problem, we here compute 
/           only the local Jacobian components over the domain, and 
/           leave the actual matrix setup/solve for the lsolve 
/           routine, since we use a Schur-complement formulation for 
/           the linear system solution.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"



int gFLDProblem::lsetup(EnzoVector *u)
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::lsetup routine\n");

  // check that the gFLDProblem has been set up
  if (!prepared) {
    fprintf(stderr,"lsetup error: gFLDProblem not yet prepared\n");
    return FAIL;
  }

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) {
    fprintf(stderr,"lsetup error: x0 vector dims do not match\n");
    return FAIL;
  }
  if (usz[1] != LocDims[1]) {
    fprintf(stderr,"lsetup error: x1 vector dims do not match\n");
    return FAIL;
  }
  if (usz[2] != LocDims[2]) {
    fprintf(stderr,"lsetup error: x2 vector dims do not match\n");
    return FAIL;
  }
  if (usz[3] != (2+Nchem)) {
    fprintf(stderr,"lsetup error: nspecies do not match\n");
    return FAIL;
  }
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) {
    fprintf(stderr,"lsetup error: x0 vector sizes do not match\n");
    return FAIL;
  }
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) {
    fprintf(stderr,"lsetup error: x1 vector sizes do not match\n");
    return FAIL;
  }
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) {
    fprintf(stderr,"lsetup error: x2 vector sizes do not match\n");
    return FAIL;
  }
  
  // we compute the local components of the Jacobian system via 
  // finite-differencing.  
  //   set up temporary arrays
  EnzoVector *fval = u->clone();
  EnzoVector *utmp = u->clone();
  EnzoVector *ftmp = u->clone();

  //   compute the local rhs at the current state, new time (fval)
  if (this->LocRHS(fval,tnew,anew,adotnew,u) == FAIL) {
    fprintf(stderr,"lsetup error: LocRHS failure\n");
    return FAIL;
  }

  // determine floating-point roundoff
  float epsilon=1.0;
  while ((1.0 + epsilon*0.5) > 1.0)  epsilon*=0.5;
  
  // clear Jacobian data arrays
  int ns;
  for (ns=0; ns<(2+Nchem); ns++)
    L[ns]->constant(0.0);
  
  //   perturb each component of the current state (utmp),
  //   for each perturbation, compute the local rhs (ftmp);
  //   with these, approximate the Jacobian columns wrt that component (L)
  int iz, iy, ix, idx, ns2;
  float sigma, *uarray, *utmparray, *farray, *ftmparray, *Lblock;
  float diagtemp, diagtemp2;
  for (ns=0; ns<(2+Nchem); ns++) {

    // [re]set utmp to the current state
    utmp->copy(u);

    diagtemp = diagtemp2 = 0.0;

    // perturb the appropriate component of utmp
    uarray = u->GetData(ns);
    utmparray = utmp->GetData(ns);
    for (iz=ghZl; iz<ghZl+usz[2]; iz++) {
      for (iy=ghYl; iy<ghYl+usz[1]; iy++) {
	for (ix=ghXl; ix<ghXl+usz[0]; ix++) {
	  idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	  sigma = sqrt(epsilon)*max(fabs(uarray[idx]),1.0);
	  utmparray[idx] += sigma;

	  diagtemp += uarray[idx];
	  diagtemp2 += sigma;
	}
      }
    }

//     if (debug)
//       fprintf(stdout,"gFLDProblem::lsetup, ||u(%"ISYM")||_1 = %g, ||sigma||_1 = %g\n",
// 	      ns,diagtemp,diagtemp2);

    // compute the local rhs due to this perturbation (ftmp)
    if (this->LocRHS(ftmp,tnew,anew,adotnew,utmp) == FAIL) {
      fprintf(stderr,"lsetup error: LocRHS failure\n");
      return FAIL;
    }

    // store the resulting Jacobian approximations
    for (ns2=0; ns2<(2+Nchem); ns2++) {
      farray = fval->GetData(ns2);
      ftmparray = ftmp->GetData(ns2);
      Lblock = (L[ns2])->GetData(ns);
      for (iz=ghZl; iz<ghZl+usz[2]; iz++) {
	for (iy=ghYl; iy<ghYl+usz[1]; iy++) {
	  for (ix=ghXl; ix<ghXl+usz[0]; ix++) {
	    idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	    sigma = sqrt(epsilon)*max(fabs(uarray[idx]),1.0);
	    Lblock[idx] = dt*theta*(ftmparray[idx]-farray[idx])/sigma;
	  }
	}
      }
    }
  }

  // delete temporary vectors
  delete fval;
  delete utmp;
  delete ftmp;

  // return success
  return SUCCESS;
}
#endif

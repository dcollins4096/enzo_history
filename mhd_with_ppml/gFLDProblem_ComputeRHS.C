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
/  Single-Group, Multi-Species Gray Flux-Limited Diffusion Implicit 
/  Problem Class RHS calculation routine.
/
/  written by: Daniel Reynolds
/  date:       September, 2006
/  modified1:  
/
/  PURPOSE: Takes in EnzoVector and returns right-hand side of ODEs for 
/           relevant equations.  This routine will be called repeatedly, 
/           so it should NOT allocate any memory.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"


/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);



int gFLDProblem::ComputeRHS(EnzoVector *rhsvec, float time, EnzoVector *u) 
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::ComputeRHS routine\n");

  // get local mesh description and check input vector sizes
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) {
    fprintf(stderr,"ComputeRHS error: x0 vector dims do not match\n");
    return FAIL;  }
  if (usz[1] != LocDims[1]) {
    fprintf(stderr,"ComputeRHS error: x1 vector dims do not match\n");
    return FAIL;  }
  if (usz[2] != LocDims[2]) {
    fprintf(stderr,"ComputeRHS error: x2 vector dims do not match\n");
    return FAIL;  }
  if (usz[3] != (2+Nchem)) {
    fprintf(stderr,"ComputeRHS error: nspecies dims do not match\n");
    return FAIL;  }
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) {
    fprintf(stderr,"ComputeRHS error: x0 vector sizes do not match\n");
    return FAIL;  }
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) {
    fprintf(stderr,"ComputeRHS error: x1 vector sizes do not match\n");
    return FAIL;  }
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) {
    fprintf(stderr,"ComputeRHS error: x2 vector sizes do not match\n");
    return FAIL;
  }

  // set a, adot to correct time-level values
  FLOAT a = 1.0, adot = 0.0;
  if (ComovingCoordinates) {
    if (time == told) {
      a = aold;
      adot = adotold;
    }
    else if (time == tnew) {
      a = anew;
      adot = adotnew;
    }
    else {
      fprintf(stderr,"LocRHS Warning: time not one of told, tnew\n");
      if (CosmologyComputeExpansionFactor(time, &a, &adot) == FAIL) {
	fprintf(stderr,"LocRHS Error in CosmologyComputeExpansionFactor\n");
	return FAIL;
      }
    }
  }


  // local rhs for all equations
  if (this->LocRHS(rhsvec, time, a, adot, u) == FAIL) {
    fprintf(stderr,"ComputeRHS error: LocRHS failure\n");
    return FAIL;
  }

  // diffusive rhs for radiation energy equation
  EnzoVector *tmp1 = u->clone();
  tmp1->constant(0.0);
  float *tmp1_E = tmp1->GetData(0);
  float *Er = u->GetData(0);
  float *Er0 = U0->GetData(0);
  // note: sigmaA, Temp have already been filled in by LocRHS
  // zero out sigmaS for now (ignore extinction due to scattering)
  for (int i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)    sigmaS[i] = 0.0;
  if (this->DiffRHS(tmp1_E, Er, Er0, Temp, sigmaA, sigmaS, &a) != SUCCESS) {
    fprintf(stderr,"ComputeRHS: Error in DiffRHS routine\n");
    return FAIL;
  }

//   float rnorm = tmp1->rmsnorm();
//   float inorm = tmp1->infnorm();
//   if (debug)
//     fprintf(stdout,"  diffusive ||rhs||_2 = %g, ||rhs||_inf = %g\n",rnorm,inorm);

  //   combine pieces together and delete temporary storage
  rhsvec->axpy_component(1.0,tmp1,0);
  delete tmp1;

//   rnorm = rhsvec->rmsnorm();
//   inorm = rhsvec->infnorm();
//   if (debug)
//     fprintf(stdout,"  current ||rhs||_2 = %g, ||rhs||_inf = %g\n",rnorm,inorm);

  // return success
  return SUCCESS;
}

#endif

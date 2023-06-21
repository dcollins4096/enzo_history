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
/  Gray Flux-Limited Diffusion Implicit Problem Class local rhs
/  calculation routine.
/
/  written by: Daniel Reynolds
/  date:       September, 2006
/  modified1:  
/
/  PURPOSE: Takes in EnzoVector and returns local time-fixed rhs for 
/           relevant equations.  This routine will be called repeatedly, 
/           so it should NOT allocate any memory.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"
#include "CosmologyParameters.h"

/* function prototypes */
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

 
 

int gFLDProblem::LocRHS(EnzoVector *locrhs, float time, 
			FLOAT a, FLOAT adot, EnzoVector *u) 
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::LocRHS routine\n");

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) {
    fprintf(stderr,"LocRHS error: x0 vector dims do not match\n");
    return FAIL;
  }
  if (usz[1] != LocDims[1]) {
    fprintf(stderr,"LocRHS error: x1 vector dims do not match\n");
    return FAIL;
  }
  if (usz[2] != LocDims[2]) {
    fprintf(stderr,"LocRHS error: x2 vector dims do not match\n");
    return FAIL;
  }
  if (usz[3] != (2+Nchem)) {
    fprintf(stderr,"LocRHS error: nspecies dims do not match\n");
    return FAIL;
  }
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) {
    fprintf(stderr,"LocRHS error: x0 vector sizes do not match\n");
    return FAIL;
  }
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) {
    fprintf(stderr,"LocRHS error: x1 vector sizes do not match\n");
    return FAIL;
  }
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) {
    fprintf(stderr,"LocRHS error: x2 vector sizes do not match\n");
    return FAIL;
  }

  // extract fluid energy, radiation energy and chemistry arrays
  float *Er = u->GetData(0);
  float *ec = u->GetData(1);
  float *n_HI, *n_HeI, *n_HeII;
  if (Nchem == 1) {
    n_HI   = u->GetData(2);
    n_HeI  = NULL;
    n_HeII = NULL;
  }
  else if (Nchem == 3) {
    n_HI   = u->GetData(2);
    n_HeI  = u->GetData(3);
    n_HeII = u->GetData(4);
  }
  else {
    fprintf(stderr,"LocRHS ERROR: only valid for Nchem = {1,3}\n");
    return FAIL;
  }

  // initialize output to zero
  locrhs->constant(0.0);

  // If using cosmology, get units
  float TempUnits, DensityUnits, LenUnits, VelUnits, TimeUnits, aUnits;
  TempUnits = DensityUnits = LenUnits = VelUnits = TimeUnits = aUnits = 1;
  if (ComovingCoordinates) {
    if (CosmologyGetUnits(&DensityUnits, &LenUnits, &TempUnits, 
			  &TimeUnits, &VelUnits, time) == FAIL) {
      fprintf(stderr,"Error in CosmologyGetUnits.\n");
      return FAIL;
    }
    aUnits = 1.0/(1.0 + InitialRedshift);
  }

  // compute temperature over domain
  this->ComputeTemperature(Temp,time,a,u);

//   float TempNorm=0.0;
//   for (int idx=0; idx<ArrDims[0]*ArrDims[1]*ArrDims[2]; idx++)  
//     TempNorm += Temp[idx]*Temp[idx];
//   if (debug)
//     fprintf(stdout,"  current ||Temp||_2 = %g\n",TempNorm);

 
  // compute local rhs for fluid energy correction equation
  float *locrhs_e = locrhs->GetData(1);
  if (this->LocEcRHS(locrhs_e, ec, Er, Temp, n_HI, n_HeI, 
		     n_HeII, &a, &adot, &aUnits, &DensityUnits, 
		     &TimeUnits, &LenUnits) != SUCCESS) {
    fprintf(stderr,"LocRHS: Error in LocEcRHS routine\n");
    return FAIL;
  }

//   float rnorm = locrhs->rmsnorm();
//   float inorm = locrhs->infnorm();
//   if (debug)
//     fprintf(stdout,"  current ||rhs||_2 = %g, ||rhs||_inf = %g\n",rnorm,inorm);


  // compute local rhs for chemical species equations
  //    set pointers to the rhs vectors
  float *rhs_HI, *rhs_HeI, *rhs_HeII;
  if (Nchem == 1) {
    rhs_HI = locrhs->GetData(2);            // Hydrogen-only case
    rhs_HeI  = NULL;
    rhs_HeII = NULL;
  }
  else if (Nchem == 3) {
    rhs_HI = locrhs->GetData(2);            // Hydrogen and Helium case
    rhs_HeI  = locrhs->GetData(3);
    rhs_HeII = locrhs->GetData(4);
  }
  if (this->LocNiRHS(rhs_HI, rhs_HeI, rhs_HeII, n_HI, n_HeI, 
		     n_HeII, Er, Temp, &a, &adot) != SUCCESS) {
    fprintf(stderr,"LocRHS: Error in LocNiRHS routine\n");
    return FAIL;
  }

//   rnorm = locrhs->rmsnorm();
//   inorm = locrhs->infnorm();
//   if (debug)
//     fprintf(stdout,"  current ||rhs||_2 = %g, ||rhs||_inf = %g\n",rnorm,inorm);


  // local rhs for radiation energy equation
  //   compute local opacity
  float *Kappa = sigmaA;
  if (Nchem == 0) {
    //     Error message: need chemistry-independent opacity
    fprintf(stderr,"LocRHS ERROR: opacity requires chemistry coupling\n");
    return FAIL;
  }
  if (this->Opacity(Kappa, n_HI, n_HeI, n_HeII, &a) != SUCCESS) {
    fprintf(stderr,"LocRHS: Error in Opacity routine\n");
    return FAIL;
  }

  //   put the pieces together
  float *locrhs_E = locrhs->GetData(0);
  if (this->LocEgRHS(locrhs_E, Er, Temp, n_HI, Kappa, &a, &adot) != SUCCESS) {
    fprintf(stderr,"LocRHS: Error in LocEgRHS routine\n");
    return FAIL;
  }

//   rnorm = locrhs->rmsnorm();
//   inorm = locrhs->infnorm();
//   if (debug)
//     fprintf(stdout,"  current ||rhs||_2 = %g, ||rhs||_inf = %g\n",rnorm,inorm);

  // return success
  return SUCCESS;
}

#endif

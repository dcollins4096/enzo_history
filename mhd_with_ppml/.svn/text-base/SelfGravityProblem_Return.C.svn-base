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
/  Single-Group, Single-Species Flux-Limited-Diffusion Implicit 
/  Problem Class return routine
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: completely pointless here, but provides example that 
/           may be helpful for implicit problems in other contexts.
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"



int SelfGravityProblem::Return(EnzoVector *sol, float *potential) 
{

  // the Enzo arrays are currently attached to the U0 vector, so 
  // unless the updated solution must be copied into U0, there 
  // is nothing to do here.  However, if this is called, just 
  // extract solution and place it in *potential array.
  EnzoVector *PotVec = sol->clone(&potential);
  PotVec->copy(sol);
  delete PotVec;
  
  // return success
  return SUCCESS;
}
#endif

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
/  Self-Gravity Problem linear Newton system solution function
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: Called by implicit solver to notify the Problem of 
/           updates to the current state (given in the vector u), 
/           so that the linear Newton matrix J(u) may be updated 
/           if necessary.  
/
/           Since the self-gravity problem is linear, there is no 
/           Jacobian update required.
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"



int SelfGravityProblem::lsetup(EnzoVector *u)
{

  // check that the SelfGravityProblem has been set up
  if (!prepared) {
    fprintf(stderr,"lsetup error: SelfGravityProblem not yet prepared\n");
    return FAIL;
  }
  
  // since problem is linear, and matrix/BCs set up in 
  // SelfGravityProblem::Setup, there is nothing to do here.

  // return success
  return SUCCESS;
}
#endif

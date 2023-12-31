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
/  Line-search step size
/
/  written by: Daniel R. Reynolds
/  date:       March, 2006
/  modified1:
/
/  PURPOSE:  Performs the 1D line-search to find lambda in the 
/            Newton correction  x = x + lambda*s, 0<lambda<=1
/
************************************************************************/
#include "InexactNewton.h"
#include "global_data.h"


int InexactNewtonSolver::LinesearchStepSize(ImplicitProblemABC *prob, 
					    EnzoVector *x, EnzoVector *s)
{

//   if (debug)
//     fprintf(stdout,"Entering InexactNewtonSolver::LinesearchStepSize routine\n");

  // check whether linesearch needed
  LSsize = 1.0;
  if (tmp1->copy(x) != SUCCESS) {
    fprintf(stderr, "Error in EnzoVector copy routine!\n");
    return FAIL;
  }
  if (tmp1->axpy(-LSsize, s) != SUCCESS) {
    fprintf(stderr, "Error in EnzoVector axpy routine!\n");
    return FAIL;
  }
  if (prob->nlresid(fvec, tmp1) != SUCCESS) {
    fprintf(stderr, "Error in Problem nlresid routine!\n");
    return FAIL;
  }
  fnorm = fvec->rmsnorm();

  if (debug)
    fprintf(stdout, "   lambda = %g  ->  fnorm = %g, fnorm0 = %g\n", 
	    LSsize, fnorm, fnorm0);

  if (fnorm < fnorm0)  return SUCCESS;
  
  // perform backtracking line-search
  while (LSsize > LSmin) {

    // reduce line-search step length
    LSsize *= 0.5;
    
    // check residual with this step length
    if (tmp1->copy(x) != SUCCESS) {
      fprintf(stderr, "Error in EnzoVector copy routine!\n");
      return FAIL;
    }
    if (tmp1->axpy(-LSsize, s) != SUCCESS) {
      fprintf(stderr, "Error in EnzoVector axpy routine!\n");
      return FAIL;
    }
    if (prob->nlresid(fvec, tmp1) != SUCCESS) {
      fprintf(stderr, "Error in Problem nlresid routine!\n");
      return FAIL;
    }
    fnorm = fvec->rmsnorm();

    if (debug)
      fprintf(stdout,"   lambda = %g  ->  fnorm = %g, fnorm0 = %g\n",
	      LSsize, fnorm, fnorm0);

    if (fnorm < fnorm0)  return SUCCESS;
  }
  
  // if we've made it to this point, then LSsize < LSmin, so return error
  fprintf(stderr,"ERROR: Line-search step below minimum step size\n");
  return FAIL;
}

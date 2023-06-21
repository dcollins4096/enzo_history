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
/  Implicit Problem Abstract Base Class, used in conjunction with 
/  nonlinear implicit solver
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: This class defines problem-specific functions required for 
/  any implicit nonlinear solve (see InexactNewton.h for additional 
/  solver information).  For a given problem, the user must create a 
/  derived class that instantiates all of these operations on the 
/  appropriate data structures.
/
************************************************************************/

#ifndef IMPLICIT_PROBLEM_ABSTRACT_BASE_CLASS_DEFINED__
#define IMPLICIT_PROBLEM_ABSTRACT_BASE_CLASS_DEFINED__

/* #include "typedefs.h" */
#include "EnzoVector.h"


class ImplicitProblemABC 
{

 public:

  // Destructor
  //  virtual ~ImplicitProblemABC() = 0;
  
  // Problem-defining nonlinear residual operations
  virtual int nlresid(EnzoVector *fu, EnzoVector *u) = 0;
  
  // Problem-specific Linear system setup function, sets up the 
  //   linear Newton system matrix J(u) given an updated state u.
  virtual int lsetup(EnzoVector *u) = 0;
  
  // Problem-specific Linear solver function 
  //   solves ||J(u)*s - b|| to tolerance delta
  virtual int lsolve(EnzoVector *s, EnzoVector *b, 
		     EnzoVector *u, float delta) = 0;
  
};
  
#endif

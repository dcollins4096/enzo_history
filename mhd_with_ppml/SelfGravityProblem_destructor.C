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
/  Self Gravity Problem destructor routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Frees all memory allocated for the self-gravity problem.
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"


SelfGravityProblem::~SelfGravityProblem()
{
  // delete HYPRE objects
  if (AInit==1) HYPRE_SStructMatrixDestroy(A);
  HYPRE_SStructGraphDestroy(graph);
  HYPRE_SStructStencilDestroy(stencil);
  HYPRE_SStructGridDestroy(grid);

  // delete local EnzoVectors
  delete rhorhs;
}
#endif

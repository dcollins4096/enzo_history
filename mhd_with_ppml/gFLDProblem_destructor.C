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
/  Gray Flux-Limited Diffusion Implicit Problem Class destructor routine
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: Frees all memory allocated for the implicit FLD problem.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"



gFLDProblem::~gFLDProblem()
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::destructor routine\n");

  // delete HYPRE objects
  HYPRE_SStructGraphDestroy(graph);
  HYPRE_SStructStencilDestroy(stencil);
  HYPRE_SStructGridDestroy(grid);

  // delete EnzoVectors
  int i, j;
  delete U0;
  delete rhs0;
  delete rhs;
  delete[] Temp;
  delete[] sigmaA;
  delete[] sigmaS;
  for (i=0; i<(Nchem+2); i++)
    delete L[i];
  delete[] L;

  // delete boundary condition arrays
  for (i=0; i<3; i++)
    for (j=0; j<2; j++)
      if (EBdryVals[i][j] != NULL)  delete[] EBdryVals[i][j];

}
#endif

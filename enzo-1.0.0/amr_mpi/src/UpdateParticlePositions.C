/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  UPDATE PARTICLE POSITIONS
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"

int UpdateParticlePositions(grid *Grid)
{
  float dt = Grid->ReturnTimeStep();

  /* 1) v(n) --> v(n+1/2) with a(n+1/2) */

  Grid->DebugCheck("UpdateParticlePosition step 1");
  if (Grid->UpdateParticleVelocity(0.5*dt) == FAIL) {
    fprintf(stderr, "Error in grid->UpdateParticleVelocity./\n");
    return FAIL;
  }

  /* 2) x(n) --> x(n+1) with v(n+1/2) 
        (problem 23 is TestGravity, so don't move the particles). */

  Grid->DebugCheck("UpdateParticlePosition step 2");
  if (ProblemType != 23)
    if (Grid->UpdateParticlePosition(dt) == FAIL) {
      fprintf(stderr, "Error in grid->UpdateParticlePosition./\n");
      return FAIL;
    }

  /* 3) v(n+1/2) --> v(n+1) with a(n+1/2) */

  Grid->DebugCheck("UpdateParticlePosition step 3");
  if (Grid->UpdateParticleVelocity(0.5*dt) == FAIL) {
    fprintf(stderr, "Error in grid->UpdateParticleVelocity./\n");
    return FAIL;
  }
  

  return SUCCESS;
}

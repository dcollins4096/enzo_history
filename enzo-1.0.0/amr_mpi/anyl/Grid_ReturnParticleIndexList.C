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
/  GRID CLASS (FIND PARTICLES WITHIN SPECIFIED RADIUS OF GIVEN POSITION)
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "StarParticleData.h"

int grid::ReturnParticleIndexList(int *ParticlesFound, int *ParticleNumberList,
				  FILE *fptr)
{
  
  int i;

  /* Loop over particles. */

  if (NumberOfParticleAttributes == 0)
    for (i = 0; i < NumberOfParticles; i++) {
      ParticleNumberList[(*ParticlesFound)++] = ParticleNumber[i];
      fprintf(fptr, "%f %f %f\n", ParticlePosition[0][i], 
	      ParticlePosition[1][i], ParticlePosition[2][i]);
    } 
  else {
    for (i = 0; i < NumberOfParticles; i++) 
      if (ParticleAttribute[0][i] < 0) {
	ParticleNumberList[(*ParticlesFound)++] = ParticleNumber[i];
	fprintf(fptr, "%f %f %f\n", ParticlePosition[0][i], 
		ParticlePosition[1][i], ParticlePosition[2][i]);
      }
  }

  return SUCCESS;
}

/***********************************************************************
/
/  GRID CLASS (FIND PARTICLES WHICH MATCH THE PROVIDED LIST OF NUMBERS
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
 
 
int grid::FindMatchingParticles(int ParticlesFound, int *ParticleNumberList,
			       float LeftEdge[], float RightEdge[], FILE *fptr)
{
 
  int i, dim, LowerIndex, UpperIndex, MidPoint;
 
  /* Loop over particles. */
 
  for (i = 0; i < NumberOfParticles; i++) {
 
    /* Find place in sorted ParticleNumberList by bisection. */
 
    LowerIndex = -1;
    UpperIndex = ParticlesFound;
 
    while (UpperIndex-LowerIndex > 1) {
      MidPoint = (LowerIndex + UpperIndex)/2;
      if (ParticleNumber[i] > ParticleNumberList[MidPoint])
	LowerIndex = MidPoint;
      else
	UpperIndex = MidPoint;
    }
 
    /* If found, set left/right edge and output */
 
    if (ParticleNumber[i] == ParticleNumberList[LowerIndex] ||
	ParticleNumber[i] == ParticleNumberList[UpperIndex]  ) {
 
      for (dim = 0; dim < GridRank; dim++) {
	LeftEdge[dim] = min(LeftEdge[dim], ParticlePosition[dim][i]);
	RightEdge[dim] = max(RightEdge[dim], ParticlePosition[dim][i]);
      }
 
      fprintf(fptr, "%"FSYM" %"FSYM" %"FSYM"\n", ParticlePosition[0][i],
	      ParticlePosition[1][i], ParticlePosition[2][i]);
 
    } // end: if (index != -1)
 
  } // end: loop over particles
 
  return SUCCESS;
}

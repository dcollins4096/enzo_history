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
 
 
int grid::FindParticlesWithinSphere(float SphereCenter[], float SphereRadius,
				  int *ParticlesFound, int *ParticleNumberList)
{
 
  int i, dim;
  float radius2, delx, dely, delz, DomainSize[MAX_DIMENSION];
 
  /* Quick check to see if sphere overlaps this grid. */
 
  for (dim = 0; dim < GridRank; dim++) {
    if (SphereCenter[dim] - SphereRadius > GridRightEdge[dim] ||
	SphereCenter[dim] + SphereRadius < GridLeftEdge[dim]   )
      return SUCCESS;
    DomainSize[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
  }
 
  /* Loop over particles. */
 
  dely = delz = 0;
  for (i = 0; i < NumberOfParticles; i++) {
 
    delx = fabs(ParticlePosition[0][i] - SphereCenter[0]);
    if (delx > 0.5*DomainSize[0]) delx = DomainSize[0] - delx;
    if (GridRank > 1) {
      dely = fabs(ParticlePosition[1][i] - SphereCenter[1]);
      if (dely > 0.5*DomainSize[1]) dely = DomainSize[1] - dely;
    }
    if (GridRank > 2) {
      delz = fabs(ParticlePosition[2][i] - SphereCenter[2]);
      if (delz > 0.5*DomainSize[2]) delz = DomainSize[2] - delz;
    }
 
    radius2 = delx*delx + dely*dely + delz*delz;
 
    if (radius2 <= SphereRadius*SphereRadius)
      ParticleNumberList[(*ParticlesFound)++] = ParticleNumber[i];
 
  }
 
  return SUCCESS;
}

/***********************************************************************
/
/  GIVEN A POSITION, FIND THE ROOT GRID
/
/  written by: John Wise
/  date:       November, 2005
/  modified1:
/
/ PURPOSE: 
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

#ifdef CONFIG_BFLOAT_4
#define ROUNDOFF 1e-6f
#endif
#ifdef CONFIG_BFLOAT_8
#define ROUNDOFF 1e-12
#endif
#ifdef CONFIG_BFLOAT_16
#define ROUNDOFF 1e-16
#endif

int FindRootGrid(int &dummy, grid **Grids0, int nGrids0, 
		 FLOAT rx, FLOAT ry, FLOAT rz, FLOAT ux, FLOAT uy, FLOAT uz)
{

  if (rx <= DomainLeftEdge[0] || rx >= DomainRightEdge[0] ||
      ry <= DomainLeftEdge[1] || ry >= DomainRightEdge[1] ||
      rz <= DomainLeftEdge[2] || rz >= DomainRightEdge[2]) {
    dummy = INT_UNDEFINED;
    return SUCCESS;
  }

  if (NumberOfProcessors == 1) {
    dummy = 0;
    return SUCCESS;
  }

  int Rank, Dims[MAX_DIMENSION], i;
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION], Bump[MAX_DIMENSION];
  FLOAT BumpPos[MAX_DIMENSION];

  Bump[0] = ROUNDOFF * sign(ux);
  Bump[1] = ROUNDOFF * sign(uy);
  Bump[2] = ROUNDOFF * sign(uz);

  // Make sure that bumped position is always within the domain

  BumpPos[0] = max(min(rx+Bump[0],DomainRightEdge[0]-ROUNDOFF), ROUNDOFF);
  BumpPos[1] = max(min(ry+Bump[1],DomainRightEdge[1]-ROUNDOFF), ROUNDOFF);
  BumpPos[2] = max(min(rz+Bump[2],DomainRightEdge[2]-ROUNDOFF), ROUNDOFF);

  for (i = 0; i < nGrids0; i++) {
    Grids0[i]->ReturnGridInfo(&Rank, Dims, Left, Right);

    if (BumpPos[0] > Left[0] && BumpPos[0] <= Right[0] &&
	BumpPos[1] > Left[1] && BumpPos[1] <= Right[1] &&
	BumpPos[2] > Left[2] && BumpPos[2] <= Right[2]) {
      dummy = i;
      return SUCCESS;
    }
    
  }

  fprintf(stderr, "FindRootGrid: Couldn't find root grid.\n"
	  "x y z = %"FSYM" %"FSYM" %"FSYM"\n", 
	  BumpPos[0], BumpPos[1], BumpPos[2]);
  dummy = nGrids0-1;
  ENZO_FAIL("");

}

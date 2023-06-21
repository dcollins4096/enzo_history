/***********************************************************************
/
/  GRID CLASS (COLLECT SUMMARY INFORMATION)
/
/  written by: Greg Bryan
/  date:       September, 1996
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


int grid::CollectGridInformation(int &GridMemory, float &GridVolume,
				 int &NumberOfCells, float &AxialRatio)
{

  GridMemory    = NumberOfBaryonFields*sizeof(float);
  GridVolume    = 1;
  NumberOfCells = 1;
  int MaxDim    = 1, MinDim = GridDimension[0];

  for (int dim = 0; dim < GridRank; dim++) {
    GridMemory    *= GridDimension[dim];
    GridVolume    *= (GridRightEdge[dim]   - GridLeftEdge[dim]  )/
                     (DomainRightEdge[dim] - DomainLeftEdge[dim]);
    NumberOfCells *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    MaxDim = max(MaxDim, (GridEndIndex[dim] - GridStartIndex[dim] + 1));
    MinDim = min(MinDim, (GridEndIndex[dim] - GridStartIndex[dim] + 1));
  }

  GridMemory = GridMemory + NumberOfParticles*
               (sizeof(float)*(GridRank+2+NumberOfParticleAttributes) + 
                sizeof(FLOAT)*GridRank);

  AxialRatio = float(MaxDim)/float(MinDim);

  return SUCCESS;

}

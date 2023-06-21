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
				 int &CellsActive, float &AxialRatio,
				 int &CellsTotal)
{

  GridMemory    = NumberOfBaryonFields*sizeof(float);
  GridVolume    = 1;
  CellsActive   = 1;
  CellsTotal    = 1;
  int MaxDim    = 1, MinDim = GridDimension[0];

  for (int dim = 0; dim < GridRank; dim++) {
    int DimActive = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    GridMemory    *= GridDimension[dim];
    GridVolume    *= (GridRightEdge[dim]   - GridLeftEdge[dim]  )/
                     (DomainRightEdge[dim] - DomainLeftEdge[dim]);
    CellsActive *= DimActive;
    CellsTotal *= GridDimension[dim];
    MaxDim = max(MaxDim, DimActive);
    MinDim = min(MinDim, DimActive);
  }

  GridMemory = GridMemory + NumberOfParticles*
               (sizeof(float)*(GridRank+2+NumberOfParticleAttributes) + 
                sizeof(FLOAT)*GridRank);

  AxialRatio = float(MaxDim)/float(MinDim);

  return SUCCESS;

}

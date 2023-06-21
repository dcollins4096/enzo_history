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
/  GRID CLASS (COMPUTES GRID DERIVED QUANTITES SUCH AS BOUNDARY FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//  Compute derived quantites
//

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void grid::PrepareGridDerivedQuantities()
{

  /*  Baryons only: set up field quantities and allocate fields
      (we assume here the grid is uniform in each dimension) */

  FLOAT AllCellWidth, GridLeftIncludingBoundary;
  int dim, j;

  for (dim = 0; dim < GridRank; dim++) {

    /* create cell position descriptors */

    CellLeftEdge[dim] = new FLOAT[GridDimension[dim]];
    CellWidth[dim]    = new FLOAT[GridDimension[dim]];

    /* assuming uniform grid, compute the cell width and the left edge
       including the ghostzones (for this dimension) */

    AllCellWidth    = (GridRightEdge[dim] - GridLeftEdge[dim])/
                      FLOAT(GridEndIndex[dim] - GridStartIndex[dim] + 1);
    GridLeftIncludingBoundary = GridLeftEdge[dim] - 
                                AllCellWidth*FLOAT(GridStartIndex[dim]);

    /* based on these quantities, set the cell position descriptors */

    for (j = 0; j < GridDimension[dim]; j++) {
      CellLeftEdge[dim][j] = GridLeftIncludingBoundary + AllCellWidth*FLOAT(j);
      CellWidth[dim][j]    = AllCellWidth;
    }

    /* Set up the BoundaryFluxes descriptors (for each dimension [dim], there
       are two faces (left, right) that are each described by two vectors
       with indicies [j].  We set all the starting vectors to the left-most
       corner and the ending vectors to the right-most corner and then
       correct. */

    for (j = 0; j < GridRank; j++) {
      BoundaryFluxes.LeftFluxStartGlobalIndex[j][dim] = 
	nlongint(( GridLeftEdge[dim] - DomainLeftEdge[dim]) / AllCellWidth);
      BoundaryFluxes.LeftFluxEndGlobalIndex[j][dim] = max(
	nlongint((GridRightEdge[dim] - DomainLeftEdge[dim]) / AllCellWidth)-1,
	BoundaryFluxes.LeftFluxStartGlobalIndex[j][dim]);
      BoundaryFluxes.RightFluxStartGlobalIndex[j][dim] = 
	nlongint(( GridLeftEdge[dim] - DomainLeftEdge[dim]) / AllCellWidth);
      BoundaryFluxes.RightFluxEndGlobalIndex[j][dim] = max(
	nlongint((GridRightEdge[dim] - DomainLeftEdge[dim]) / AllCellWidth)-1,
	BoundaryFluxes.RightFluxStartGlobalIndex[j][dim]);
      }

    /* correct so vectors point to correct position */

    BoundaryFluxes.LeftFluxEndGlobalIndex[dim][dim] = 
      BoundaryFluxes.LeftFluxStartGlobalIndex[dim][dim];
    BoundaryFluxes.RightFluxStartGlobalIndex[dim][dim] = 
      BoundaryFluxes.RightFluxEndGlobalIndex[dim][dim];
      
  }

}


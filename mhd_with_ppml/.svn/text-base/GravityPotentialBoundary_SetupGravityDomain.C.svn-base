/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Daniel R. Reynolds                                         *
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRAVITY POTENTIAL BOUNDARY CLASS (SET UP GRAVITY DOMAIN)
/
/  written by: Daniel R. Reynolds
/  date:       December, 2005
/  modified1:
/
/  PURPOSE: This routine is used to set up domain information for the 
/           self-gravity potential field solve.  Information computed 
/           includes boundary types, per-processor domain edges, 
/           domain extents, as well as information as to whether this 
/           grid owns part of any external boundary.
/
/           Isolating BCs require that the gravity domain also include 
/           the external ghost cells, as well as an additional ghost 
/           cell layer at all faces.
/
/     NOTE: This routine must be called *after* the gravity boundary 
/           types have been put into MetaData (during the problem 
/           initialization).
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "GravityPotentialBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"


#ifdef ISO_GRAV
int GravityPotentialBoundary::SetupGravityDomain(HierarchyEntry &TopGrid,
						 TopGridData &MetaData)
{

  /* find the grid corresponding to this process from the Hierarcy */
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("Error: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
    return FAIL;
  }
  
  /* set boundary rank */
  BoundaryRank = ThisGrid->GridData->GetGridRank();

  /*    get processor layout from Grid */
  int layout[3];
  layout[0] = ThisGrid->GridData->GetProcessorLayout(0);
  layout[1] = ThisGrid->GridData->GetProcessorLayout(1);
  layout[2] = ThisGrid->GridData->GetProcessorLayout(2);
  
  /*    get processor location in MPI grid */
  int location[3];
  location[0] = ThisGrid->GridData->GetProcessorLocation(0);
  location[1] = ThisGrid->GridData->GetProcessorLocation(1);
  location[2] = ThisGrid->GridData->GetProcessorLocation(2);

  /* store gravity domain boundary conditions */
  BoundaryType[0] = MetaData.GravityBoundaryFaces[0];
  BoundaryType[1] = MetaData.GravityBoundaryFaces[1];
  BoundaryType[2] = MetaData.GravityBoundaryFaces[2];

  /* check that these give appropriate values, otherwise set dim to periodic */
  if ((BoundaryType[0] < 0) || (BoundaryType[0] > 1))  BoundaryType[0] = 0;
  if ((BoundaryType[1] < 0) || (BoundaryType[1] > 1))  BoundaryType[1] = 0;
  if ((BoundaryType[2] < 0) || (BoundaryType[2] > 1))  BoundaryType[2] = 0;

  /* if dim is periodic, ensure that its BoundaryValues arrays are NULL */
  if (BoundaryType[0] == 0) {
    if (BoundaryValues[0][0] != NULL)  delete[] BoundaryValues[0][0];
    if (BoundaryValues[0][1] != NULL)  delete[] BoundaryValues[0][1];
  }
  if (BoundaryType[1] == 0) {
    if (BoundaryValues[1][0] != NULL)  delete[] BoundaryValues[1][0];
    if (BoundaryValues[1][1] != NULL)  delete[] BoundaryValues[1][1];
  }
  if (BoundaryType[2] == 0) {
    if (BoundaryValues[2][0] != NULL)  delete[] BoundaryValues[2][0];
    if (BoundaryValues[2][1] != NULL)  delete[] BoundaryValues[2][1];
  }

  
  /* set up subdomain information */
  int dim;

  /*   GravityLeftEdge gives the location of the left edge of the 
       gravity domain -- start with Enzo grid size */
  for (dim=0; dim<3; dim++)
    GravityLeftEdge[dim] = ThisGrid->GridData->GetGridLeftEdge(dim);

  /*   GravityRightEdge gives the location of the right edge of the 
       gravity domain -- start with Enzo grid size */
  for (dim=0; dim<3; dim++)
    GravityRightEdge[dim] = ThisGrid->GridData->GetGridRightEdge(dim);

  /*   LocalDimension holds the dimensions of the local gravity 
          domain: we first obtain the Enzo grid size */
  for (dim=0; dim<3; dim++)
    LocalDimension[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;

  /* GravityCellSize gives grid cell size */
  GravityCellSize = (GravityRightEdge[0]-GravityLeftEdge[0])/LocalDimension[0];

  /* set flags denoting if this processor is on the external boundary */
  for (dim=0; dim<3; dim++) {
    OwnBoundary[dim][0] = (location[dim] == 0);
    OwnBoundary[dim][1] = (location[dim] == layout[dim]-1);
  }

  /*   adjust if proc at low boundary: (add 1 for FD grid ownership) 
       for isolating domain, add ghost zones */
  int GravBufferSize = GRAVITY_BUFFER_SIZE+1;
  for (dim=0; dim<3; dim++) {
    if (OwnBoundary[dim][0]) {
      if (BoundaryType[dim] == 1) {
	LocalDimension[dim] += GravBufferSize;
	GravityLeftEdge[dim] = GravityLeftEdge[dim] 
	  - GravBufferSize*GravityCellSize;
      }
    }
  }

  /*   adjust if proc at high boundary:
          for isolating domain, add ghost zones */
  for (dim=0; dim<3; dim++) {
    if ((OwnBoundary[dim][1]) && (BoundaryType[dim] == 1)) {
      LocalDimension[dim] += GravBufferSize;
      GravityRightEdge[dim] = GravityRightEdge[dim] 
	+ GravBufferSize*GravityCellSize;
    }
  }

  /* compute global dimension information */
  for (dim=0; dim<3; dim++) {

    // get hydro. global grid size
    GlobalDimension[dim] = MetaData.TopGridDims[dim];

    // adjust for boundary conditions
    if (BoundaryType[dim] == 1)
      GlobalDimension[dim] += 8;
  }

  /* compute global index information for this gravity subdomain */
  FLOAT fCellsLeft;
  for (dim=0; dim<3; dim++) {

    /* the global indexing is easy if we're at the left edge */
    if (location[dim]==0)  LeftEdgeIndex[dim]=0;

    /* otherwise, we compute the number of intervening cells to left edge */
    else {

      /* get floating point value for number of cells */
      fCellsLeft = (GravityLeftEdge[dim] - DomainLeftEdge[dim])/GravityCellSize;

      /* round floating point value to closest integer */
      LeftEdgeIndex[dim] = lround(fCellsLeft);

      /* add 4 cells padding if isolating (gravity domain larger than hydro) */
      if (BoundaryType[dim]==1)  LeftEdgeIndex[dim] += 4;
    }

    /* add on local size to obtain right edge indices */
    RightEdgeIndex[dim] = LeftEdgeIndex[dim] + LocalDimension[dim] - 1;
  }  

  /* store setup status */
  ReadyToGo = TRUE;

  return SUCCESS;
}
#endif

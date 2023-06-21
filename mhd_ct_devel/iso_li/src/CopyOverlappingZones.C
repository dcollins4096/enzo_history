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
/  COPY OVERLAPPING GRIDS FUNCTION
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
#include "TopGridData.h"
#include "LevelHierarchy.h"

/*
void Pout( char * string, int i1 = -12345, int i2 = -12345,
           int i3 = -12345, int i4 = -12345, int i5 = -12345 );
*/

int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level)
{

  /* Loop over all grids in list (except for self). */

  LevelHierarchyEntry *Temp = LevelArray[level];
  int dccGC = 0;
  //Pout( "    Copy Overlapping Zones: Level ", level);
  while (Temp != NULL) {
    //Pout( "      Grid ", dccGC);
    if (CurrentGrid->CheckForOverlap(Temp->GridData, 
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid)
	== FAIL) {
      fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
    }
    Temp = Temp->NextGridThisLevel;
    dccGC++;
  }

  return SUCCESS;
}

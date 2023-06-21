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
/  CHECK FOR PARTICLE MASS FIELDS OVERLAPPING THIS GRIDS GRAVITATING MASS FIELD
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


int CopyOverlappingParticleMassFields(grid* CurrentGrid, 
				      TopGridData *MetaData, 
				      LevelHierarchyEntry *LevelArray[], 
				      int level)
{

  /* Loop over all grids in list (including self). */

  LevelHierarchyEntry *Temp = LevelArray[level];

  while (Temp != NULL) {
    if (CurrentGrid->CheckForOverlap(Temp->GridData,
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::AddOverlappingParticleMassField)
	== FAIL) {
      fprintf(stderr, "Error in grid->AddOverlappingParticleMassFields.\n");
    }
    Temp = Temp->NextGridThisLevel;
  }

  return SUCCESS;
}

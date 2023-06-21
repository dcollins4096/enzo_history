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
/  READ IN THE DATA HIERARCHY (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// This function reads in the data hierarchy (TopGrid)

#include <string.h>
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

/* function prototypes */

static int ReadDataGridCounter = 0;

int ReadDataHierarchy(FILE *fptr, HierarchyEntry *Grid, int GridID,
		      HierarchyEntry *ParentGrid)
{

  int TestGridID, NextGridThisLevelID, NextGridNextLevelID;

  /* Read header info for this grid */

  if (fscanf(fptr, "\nGrid = %d\n", &TestGridID) != 1) {
    fprintf(stderr, "Error reading Grid # in grid %d.\n", GridID);
    return FAIL;
  }
  if (TestGridID != GridID) {
    fprintf(stderr, "Unexpected GridID = %d while reading grid %d.\n", 
	    TestGridID, GridID);
    return FAIL;
  }

  /* Create new grid and fill out hierarchy entry. */

  Grid->GridData          = new grid;
  Grid->NextGridThisLevel = NULL;
  Grid->NextGridNextLevel = NULL;
  Grid->ParentGrid        = ParentGrid;
  //  if (ParentGrid == NULL) {
  //    if (ParallelRootGridIO == TRUE)
      Grid->GridData->SetProcessorNumber(ReadDataGridCounter++ % 
					 NumberOfProcessors);
      //    else
      //      Grid->GridData->SetProcessorNumber(ROOT_PROCESSOR);
      //  }
      //  else
      //    Grid->GridData->SetProcessorNumber
      //			    (ParentGrid->GridData->ReturnProcessorNumber());

  /* Read grid data for this grid. */

  if (Grid->GridData->ReadGrid(fptr) == FAIL) {
    fprintf(stderr, "Error in grid->ReadGrid (grid %d).\n", GridID);
    return FAIL;
  }

  /* Read pointer information for the next grid this level. */

  if (fscanf(fptr, "Pointer: Grid[%d]->NextGridThisLevel = %d\n",
	     &TestGridID, &NextGridThisLevelID) != 2) {
    fprintf(stderr, "Error reading NextGridThisLevel pointer for grid %d.\n",
	    GridID);
    return FAIL;
  }
  if (TestGridID != GridID) {
    fprintf(stderr, "GridID = %d does not match grid(1) %d.\n", 
	    TestGridID, GridID);
    return FAIL;
  }

  /* If the pointer was non-zero, then read that grid. */

  if (NextGridThisLevelID != 0) {
    Grid->NextGridThisLevel = new HierarchyEntry;
    if (ReadDataHierarchy(fptr, Grid->NextGridThisLevel, NextGridThisLevelID,
			  ParentGrid) == FAIL) {
      fprintf(stderr, "Error in ReadDataHierarchy(1).\n");
      return FAIL;
    }
  }

  /* Read pointer information for the next grid next level. */

  if (fscanf(fptr, "Pointer: Grid[%d]->NextGridNextLevel = %d\n",
	     &TestGridID, &NextGridNextLevelID) != 2) {
    fprintf(stderr, "Error reading NextGridNextLevel pointer for grid %d.\n",
	    GridID);
    return FAIL;
  }
  if (TestGridID != GridID) {
    fprintf(stderr, "GridID = %d does not match grid(2) %d.\n", 
	    TestGridID, GridID);
    return FAIL;
  }

  /* If the pointer was non-zero, then read that grid. */

  if (NextGridNextLevelID != 0) {
    Grid->NextGridNextLevel = new HierarchyEntry;
    if (ReadDataHierarchy(fptr, Grid->NextGridNextLevel, NextGridNextLevelID,
			  Grid)
	== FAIL) {
      fprintf(stderr, "Error in ReadDataHierarchy(2).\n");
      return FAIL;
    }
  }

  return SUCCESS;
}

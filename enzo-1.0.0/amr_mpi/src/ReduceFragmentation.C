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
/  REDUCE FRAGMENTATION FUNCTION
/
/  written by: Greg Bryan
/  date:       February, 1999
/  modified1:
/
/  PURPOSE:
/    This routine deletes the hierarchy and rereads it, hopefully
/    reducing memory fragmentation in the process.
/
************************************************************************/

#include <stdio.h>
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */
#include "performance.h"
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
#include "CosmologyParameters.h"
#include "error.h"

/* function prototypes */

int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
char LastFileNameWritten[MAX_LINE_LENGTH];
int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

/*  function */

int ReduceFragmentation(HierarchyEntry &TopGrid, TopGridData &MetaData,
			ExternalBoundary *Exterior, 
			LevelHierarchyEntry *LevelArray[])
{

  /* Declarations. */

  int level;
  LevelHierarchyEntry *Previous, *Temp;

  /* Delete hierarchy, level array data and grids themselves. */

  fprintf(stderr, "Fragmentation reduction: deleting...");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {

      /* Delete grid itself. */

      delete Temp->GridData;

      /* Delete hierarchy entry except the top one. */

      if (Temp->GridHierarchyEntry != &TopGrid)
	delete Temp->GridHierarchyEntry;

      Previous = Temp;
      Temp = Temp->NextGridThisLevel;

      /* Delete previous level hierarchy entry. */

      delete Previous;

    }
  }

  /* Reload data (unfortunately also reads in ExternalBoundary). */

  fprintf(stderr, "reading %s...", LastFileNameWritten);
#ifdef USE_MPI
  JBPERF_START_MPI_BARRIER("MPI_Barrier");
  CHECK_MPI_ERROR(MPI_Barrier(MPI_COMM_WORLD));
  JBPERF_STOP_MPI_BARRIER("MPI_Barrier");
#endif /* USE_MPI */
  if (ReadAllData(LastFileNameWritten, &TopGrid, 
		  MetaData, Exterior) == FAIL) {
    fprintf(stderr, "Error reloading data: %s\n", LastFileNameWritten);
    return FAIL;
  }
  AddLevel(LevelArray, &TopGrid, 0);
  fprintf(stderr, "done\n");

  /* Set top grid boundary conditions. */

  Temp = LevelArray[0];
  while (Temp != NULL) {
    if (Temp->GridData->SetExternalBoundaryValues(Exterior) == FAIL) {
      fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
      return FAIL;
    }
    if (CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, 0) 
	== FAIL) {
      fprintf(stderr, "Error in CopyOverlappingZones.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }
  
  /* Rebuild the grids from level 0. */
  
  if (RebuildHierarchy(&MetaData, LevelArray, 0) == FAIL) {
    fprintf(stderr, "Error in RebuildHierarchy.\n");
    return FAIL;
  }

  return SUCCESS;

}

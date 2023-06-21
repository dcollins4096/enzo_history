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
/  REBUILD HIERARCHY FUNCTION
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  August, 1995 by GB
/              Rewritten to rebuild an entire level (and below) at a time.
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <string.h>
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
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

/* function prototypes */

void AddLevel(LevelHierarchyEntry *LevelArray[], HierarchyEntry *Grid, 
	      int level);
int FindSubgrids(HierarchyEntry &Grid, int level);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int  ReportMemoryUsage(char *header = NULL);
int CommunicationShareGrids(HierarchyEntry *GridHierarchyPointer[], int grids);
int CommunicationLoadBalanceGrids(HierarchyEntry *GridHierarchyPointer[],
				  int NumberOfGrids);
int CommunicationTransferParticles(grid *GridPointer[], int NumberOfGrids);


/* RebuildHierarchy function */

int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level)
{

  /* declarations */

  if (debug) printf("RebuildHierarchy: level = %d\n", level);
  ReportMemoryUsage("Rebuild pos 1");

  int i, j, k, grids, grids2, subgrids;
  FLOAT ZeroVector[MAX_DIMENSION];
  LevelHierarchyEntry *Temp;
  for (i = 0; i < MAX_DIMENSION; i++)
    ZeroVector[i] = 0;

  ZLAN_STOP(0);

  /* --------------------------------------------------------------------- */
  /* For each grid on this level collect all the particles below it.
     Notice that this must be done even for static hierarchy's.  */

  HierarchyEntry **GridParent = new HierarchyEntry*[MAX_NUMBER_OF_SUBGRIDS];
  grid           **GridPointer = new grid*[MAX_NUMBER_OF_SUBGRIDS];
  grid           **ContigiousGridList = new grid*[MAX_NUMBER_OF_SUBGRIDS];

  for (i = MAX_DEPTH_OF_HIERARCHY-1; i > level; i--) {

    Temp = LevelArray[i];

    /* Find the parents (at level=level) of all the grids (at level=i). */

    grids = 0;
    while (Temp != NULL) {
      GridPointer[grids] = Temp->GridData;
      GridParent[grids] = Temp->GridHierarchyEntry->ParentGrid;
      for (j = i-1; j > level; j--)
	GridParent[grids] = GridParent[grids]->ParentGrid;
      Temp = Temp->NextGridThisLevel;
      grids++;
    }

    /* Collect all the grids with the same parent and pass them all to
       MoveAllParticles (marking which ones have already been passed). */

    for (j = 0; j < grids; j++)
      if (GridPointer[j] != NULL) {
	grids2 = 0;
	for (k = j; k < grids; k++)
	  if (GridParent[k] == GridParent[j]) {
	    ContigiousGridList[grids2++] = GridPointer[k];
	    GridPointer[k] = NULL;
	  }
	if (GridParent[j]->GridData->MoveAllParticles(grids2, 
						      ContigiousGridList) == FAIL) {
	  fprintf(stderr, "Error in grid->MoveAllParticles.\n");
	  return FAIL;
	}
      }

  } // end: loop over levels

    /* --------------------------------------------------------------------- */
    /* if this is level 0 then transfer particles between grids. */

  ReportMemoryUsage("Rebuild pos 2");
  if (level == 0) {
    grids = 0;
    Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->DebugCheck("Before TransferParticles");
      GridPointer[grids++] = Temp->GridData;
      Temp = Temp->NextGridThisLevel;
    }
    if (CommunicationTransferParticles(GridPointer, grids) == FAIL) {
      fprintf(stderr, "Error in CommunicationTransferParticles.\n");
      return FAIL;
    }
  }

  delete [] GridParent;
  delete [] GridPointer;
  delete [] ContigiousGridList;

  /* --------------------------------------------------------------------- */
  /* For dynamic hierarchies, rebuild the grid structure. */

  ReportMemoryUsage("Rebuild pos 3");

  /* Allocate temporary arrays */

  grid **ToGrids = new grid * [MAX_NUMBER_OF_SUBGRIDS];
  HierarchyEntry **GridHierarchyPointer = 
    new HierarchyEntry * [MAX_NUMBER_OF_SUBGRIDS];

  HierarchyEntry **SubgridHierarchyPointer = 
    new HierarchyEntry * [MAX_NUMBER_OF_SUBGRIDS];

  int *NumberOfOverlaps = new int [MAX_NUMBER_OF_SUBGRIDS];


  if (MetaData->StaticHierarchy == FALSE) {
    
//    if (debug) ReportMemoryUsage("Memory usage report: Rebuild 1");

    /* 1) Create a new TempLevelArray in which to keep the old grids. */

    LevelHierarchyEntry* TempLevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (i = level+1; i < MAX_DEPTH_OF_HIERARCHY; i++) {
      TempLevelArray[i] = LevelArray[i];
      LevelArray[i]     = NULL;
    }
    TempLevelArray[level] = LevelArray[level];

    /* 2) Clean up (delete excess baggage) all grids on this level and below.
          And delete the old hierarchy entries at the same time. */

    for (i = level; i < MAX_DEPTH_OF_HIERARCHY; i++) {
      Temp = TempLevelArray[i];

      while (Temp != NULL) {
	Temp->GridData->CleanUp();
	if (i > level)
	  delete Temp->GridHierarchyEntry;
	Temp = Temp->NextGridThisLevel;
      } // end: if (i > level)

    } // end: loop over levels

//    if (debug) ReportMemoryUsage("Memory usage report: Rebuild 3");

    /* 3) Rebuild all grids on this level and below.  Note: All the grids
          in LevelArray[level+] have been deleted. */

    for (i = level; i < MAX_DEPTH_OF_HIERARCHY-1; i++) {

      /* If there are no grids on this level, exit. */

      if (LevelArray[i] == NULL)
	break;

      /* 3a) Generate an array of grids on this level. */

      grids = 0;
      Temp = LevelArray[i];
      while (Temp != NULL) {
	GridHierarchyPointer[grids++] = Temp->GridHierarchyEntry;
	Temp                          = Temp->NextGridThisLevel;
      }

      /* 3b) Loop over grids creating new (but empty!) subgrids
	 (This also properly fills out the GridHierarchy tree). */

      for (j = 0; j < grids; j++)
	if (FindSubgrids(*GridHierarchyPointer[j], i) == FAIL) {
	  fprintf(stderr, "Error in FindSubgrids.\n");
	  return FAIL;
	}

      /* Create a temporary array of the new subgrids (which are on this 
	 processor) for the next step. */

      HierarchyEntry *Temp2;
      subgrids = 0;
      for (j = 0; j < grids; j++) {
	Temp2 = GridHierarchyPointer[j]->NextGridNextLevel;
	while (Temp2 != NULL) {
	  SubgridHierarchyPointer[subgrids++] = Temp2;
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }

      /* 3g) loop over parent, and copy particles to new grids
	     (all local to this processor) . */

      for (j = 0; j < grids; j++)

	if (GridHierarchyPointer[j]->NextGridNextLevel != NULL) {

	  if (debug) printf("grid %d:", j);

	  GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid(
			   NULL, ZERO_UNDER_SUBGRID_FIELD);

	  for (k = 0; k < subgrids; k++) {
	    if (GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid(
		              SubgridHierarchyPointer[k]->GridData, 
		              ZERO_UNDER_SUBGRID_FIELD, float(k+1)) == FAIL) {
	      fprintf(stderr, "Error in grid->ZeroSolutionUnderSubgrid.\n");
	      return FAIL;
	    }
	    ToGrids[k] = SubgridHierarchyPointer[k]->GridData;
	  }

	  if (GridHierarchyPointer[j]->GridData->MoveSubgridParticlesFast(
				 subgrids, ToGrids, TRUE) == FAIL) {
	    fprintf(stderr, "Error in grid->MoveSubgridParticlesFast.\n");
	    return FAIL;
	  }

	}

      /* Share the new grids amoung processors. */

      CommunicationShareGrids(GridHierarchyPointer, grids);

      /* 3c) Combine the many linked-lists of subgrids into the LevelArray
	 linked list. */

      for (j = 0; j < grids; j++)
	if (GridHierarchyPointer[j]->NextGridNextLevel != NULL)
	  AddLevel(LevelArray, GridHierarchyPointer[j]->NextGridNextLevel,i+1);

      /* 3d) Create an array of the new subgrids. */

      subgrids = 0;
      Temp = LevelArray[i+1];
      while (Temp != NULL) {
	SubgridHierarchyPointer[subgrids++] = Temp->GridHierarchyEntry;
	Temp                                = Temp->NextGridThisLevel;
      }

      /* 3e) Loop over the new subgrids and record in the old subgrids how
	 many times they are needed (number of overlaps with new subgrids). */

      int Overlap, oldgrid = 0;
      Temp = TempLevelArray[i+1];
      while (Temp != NULL) {

	NumberOfOverlaps[oldgrid] = 0;
	for (j = 0; j < subgrids; j++) {
	  if (SubgridHierarchyPointer[j]->GridData->CopyZonesFromGridCountOnly(
		                          Temp->GridData, Overlap) == FAIL) {
	    fprintf(stderr, "Error in grid->CopyZonesFromGridCountOnly.\n");
	    return FAIL;
	  }
	  NumberOfOverlaps[oldgrid] += Overlap;
	}

	if (NumberOfOverlaps[oldgrid] == 0) {
	  delete Temp->GridData;
	  Temp->GridData = NULL;
	}

	Temp = Temp->NextGridThisLevel;
	oldgrid++;
      }

      /* 3f) For each new subgrid, interpolate from parent and then copy
	 from old subgrids.  Also, copy particles (if present).  For each 
	 old subgrid, decrement the Overlap counter, deleting the grid 
	 which it reaches zero. */

      for (j = 0; j < subgrids; j++) {
	SubgridHierarchyPointer[j]->ParentGrid->GridData->DebugCheck(
						        "Rebuild parent");
	SubgridHierarchyPointer[j]->GridData->InterpolateFieldValues(
		       SubgridHierarchyPointer[j]->ParentGrid->GridData);
	SubgridHierarchyPointer[j]->GridData->DebugCheck("Rebuild child");
      }

      /* Copy from old grids. */

      for (j = 0; j < subgrids; j++) {

	oldgrid = 0;
	Temp = TempLevelArray[i+1];
	while (Temp != NULL) {

	  if (Temp->GridData != NULL) {

	    /* Copy from old subgrid. */

	    if (SubgridHierarchyPointer[j]->GridData->CopyZonesFromGrid(
                                       Temp->GridData, ZeroVector) == FAIL) {
	      fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
	      return FAIL;
	    }

	    /* Check if we can delete the old subgrid yet. */

	    SubgridHierarchyPointer[j]->GridData->CopyZonesFromGridCountOnly(
		                                  Temp->GridData, Overlap);

	    if (Overlap == TRUE)
	      if (--NumberOfOverlaps[oldgrid] <= 0) {
		delete Temp->GridData;
		Temp->GridData = NULL;
	      }

	  } // end: if (Temp->GridData != NULL)

	  /* Next old subgrid. */

	  Temp = Temp->NextGridThisLevel;
	  oldgrid++;

	} // end: loop over old subgrids

//	if (debug) ReportMemoryUsage("Memory usage report: PostBuild");

	SubgridHierarchyPointer[j]->GridData->ComputeTimeStep();

      } // end: loop over new subgrids

      /* Redistribute grids over processors to Load balance. */

      CommunicationLoadBalanceGrids(SubgridHierarchyPointer, subgrids);

      /* 3h) Clean up the LevelHierarchy entries for the old subgrids.
	     Also, we can check to see if any old subgrids were missed. */

      while (TempLevelArray[i+1] != NULL) {
	Temp = TempLevelArray[i+1]->NextGridThisLevel;

	if (TempLevelArray[i+1]->GridData != NULL) {
	  fprintf(stderr, "An old subgrid was not deleted.  Why?\n");
	  return FAIL;
	}

	/* Remove the LevelHierarchy entry for that grid. */

	delete TempLevelArray[i+1];
	TempLevelArray[i+1] = Temp;
      }

    } // end: loop over levels

    // end: if (StaticHierarchy == FALSE)
  }  else if (MetaData->StaticHierarchy == TRUE) {

  /* --------------------------------------------------------------------- */
  /* Redistribute particles: for each grid, move the particles that belong in
     it's subgrids (this only has to be done if we didn't just do a rebuild
      since the rebuild does this as it goes). */


    for (i = level; i < MAX_DEPTH_OF_HIERARCHY-1; i++) {

      /* If there are no grids on this level, exit. */

      if (LevelArray[i] == NULL)
	break;

      /* 3a) Generate an array of grids on this level. */

      grids = 0;
      Temp = LevelArray[i];
      while (Temp != NULL) {
	GridHierarchyPointer[grids++] = Temp->GridHierarchyEntry;
	Temp                          = Temp->NextGridThisLevel;
      }

      /* 3d) Create an array of the subgrids. */

      subgrids = 0;
      Temp = LevelArray[i+1];
      while (Temp != NULL) {
	SubgridHierarchyPointer[subgrids++] = Temp->GridHierarchyEntry;
	Temp                                = Temp->NextGridThisLevel;
      }

      /* 3g) loop over parent, and copy particles to new grids. */

      for (j = 0; j < grids; j++)

	if (GridHierarchyPointer[j]->NextGridNextLevel != NULL) {

	  GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid(
			   NULL, ZERO_UNDER_SUBGRID_FIELD);

	  for (k = 0; k < subgrids; k++) {
	    if (GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid(
		              SubgridHierarchyPointer[k]->GridData, 
		              ZERO_UNDER_SUBGRID_FIELD, float(k+1)) == FAIL) {
	      fprintf(stderr, "Error in grid->ZeroSolutionUnderSubgrid.\n");
	      return FAIL;
	    }
	    ToGrids[k] = SubgridHierarchyPointer[k]->GridData;
	  }

	  if (GridHierarchyPointer[j]->GridData->MoveSubgridParticlesFast(
				 subgrids, ToGrids, FALSE) == FAIL) {
	    fprintf(stderr, "Error in grid->MoveSubgridParticlesFast.\n");
	    return FAIL;
	  }

	}

      /* Set boundary conditions. */

      LevelHierarchyEntry *Temp = LevelArray[i+1];
      while (Temp != NULL) {

	if (Temp->GridData->InterpolateBoundaryFromParent
	    (Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
	  return FAIL;
	}

	Temp = Temp->NextGridThisLevel;
      }

    } // end: loop over levels

  } // end: if (StaticHierarchy == TRUE)

  /* Deallocate temporary arrays */

  delete [] ToGrids;
  delete [] GridHierarchyPointer;
  delete [] SubgridHierarchyPointer; 
  delete [] NumberOfOverlaps;

  ZLAN_STOP(1);

  /* Done for this level. */

  ReportMemoryUsage("Rebuild pos 4");
  return SUCCESS;

}

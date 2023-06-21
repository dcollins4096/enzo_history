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
/  COMMUNICATION ROUTINE: SHARE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
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
#include "error.h"

/* function prototypes */

float ReturnCPUTime();

#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_PackedGrid;
#endif

int CommunicationShareGrids(HierarchyEntry *GridHierarchyPointer[], 
			    int NumberOfGrids)
{

  if (NumberOfProcessors == 1)
    return SUCCESS;

  /* Declarations. */

  struct PackedGrid {
    int Rank;
    int Dimension[MAX_DIMENSION];
    FLOAT LeftEdge[MAX_DIMENSION];
    FLOAT RightEdge[MAX_DIMENSION];
    int NumberOfParticles;
    int ParentNumber;
  };
  int i;

  /* Count the subgrids on this processor. */

  int GridsToSend = 0;
  HierarchyEntry *Temp;
  for (i = 0; i < NumberOfGrids; i++) {
    Temp = GridHierarchyPointer[i]->NextGridNextLevel;
    while (Temp != NULL) {
      GridsToSend++;
      Temp = Temp->NextGridThisLevel;
    }
  }
//  printf("ShareGrids (%d): NumberOfGrids = %d GridsToSend = %d\n", 
//	 MyProcessorNumber, NumberOfGrids, GridsToSend);

  /* Allocate an array of packed subgrids and fill it out. */

  int Counter = 0;
  PackedGrid *SendList = NULL;
  if (GridsToSend > 0) {
    SendList = new PackedGrid[GridsToSend];
    for (i = 0; i < NumberOfGrids; i++) {
      Temp = GridHierarchyPointer[i]->NextGridNextLevel;
      while (Temp != NULL) {
	Temp->GridData->ReturnGridInfo(&SendList[Counter].Rank,
				     SendList[Counter].Dimension,
				     SendList[Counter].LeftEdge,
				     SendList[Counter].RightEdge);
	SendList[Counter].NumberOfParticles = 
	  Temp->GridData->ReturnNumberOfParticles();
	SendList[Counter++].ParentNumber = i;
	Temp = Temp->NextGridThisLevel;
      }
    }
  }

  /* Allocate the array to receive subgrids. */

  PackedGrid *SharedList = NULL;
  int NumberOfSharedGrids = 0;

#ifdef USE_MPI

  /* Generate a new MPI data type corresponding to the PackedGrid struct. */

  if (FirstTimeCalled) {
    CHECK_MPI_ERROR(MPI_Type_contiguous(sizeof(PackedGrid), MPI_BYTE, 
					&MPI_PackedGrid));
    CHECK_MPI_ERROR(MPI_Type_commit(&MPI_PackedGrid));
    FirstTimeCalled = FALSE;
  }
  int *SharedListCount = new int[NumberOfProcessors],
      *SharedListDisplacements = new int[NumberOfProcessors];

  /* Get counts from each processor to allocate buffers. */

  float time1 = ReturnCPUTime();
  ZLAN_START
  JBPERF_START_MPI_GATHER("MPI_Allgather",1,MPI_INT);
  CHECK_MPI_ERROR(MPI_Allgather(&GridsToSend, 1, MPI_INT, SharedListCount, 1, 
				MPI_INT, MPI_COMM_WORLD));
  JBPERF_STOP_MPI_GATHER("MPI_Allgather",1,MPI_INT);

  /* Allocate buffers and generated displacement list. */

  for (i = 0; i < NumberOfProcessors; i++) {
    SharedListDisplacements[i] = NumberOfSharedGrids;
    NumberOfSharedGrids += SharedListCount[i];
  }
  SharedList = new PackedGrid[NumberOfSharedGrids];

  /* Perform sharing operation. */

  JBPERF_START_MPI_GATHER("MPI_Allgatherv",GridsToSend,MPI_PackedGrid);

  CHECK_MPI_ERROR(MPI_Allgatherv(SendList, GridsToSend, MPI_PackedGrid, 
				 SharedList,SharedListCount, 
				 SharedListDisplacements, MPI_PackedGrid,
				 MPI_COMM_WORLD));
  JBPERF_STOP_MPI_GATHER("MPI_Allgatherv",GridsToSend,MPI_PackedGrid);

  ZLAN_STOP_GLOBAL(13);

  CommunicationTime += ReturnCPUTime() - time1;

  delete [] SharedListCount;
  delete [] SharedListDisplacements;
  
#endif /* USE_MPI */

  /* Unpack the subgrids. */

  HierarchyEntry *PreviousGrid = NULL, *ThisGrid, *SubgridParent;
  for (i = 0; i < NumberOfSharedGrids; i++) {

    /* If this subgrid has a parent not on this processor, then allocate it. */

    SubgridParent = GridHierarchyPointer[SharedList[i].ParentNumber];
    if (SubgridParent->GridData->ReturnProcessorNumber() 
	!= MyProcessorNumber) {

      /* create hierarchy entry */
      
      ThisGrid = new HierarchyEntry;
      
      /* set hierarchy values */

      if (PreviousGrid != NULL)
	if (PreviousGrid->ParentGrid != SubgridParent)
	  PreviousGrid = NULL;

      if (PreviousGrid == NULL)
	SubgridParent->NextGridNextLevel = ThisGrid;
      else
	PreviousGrid->NextGridThisLevel = ThisGrid;
      ThisGrid->NextGridNextLevel = NULL;
      ThisGrid->NextGridThisLevel = NULL;
      ThisGrid->ParentGrid        = SubgridParent;
      PreviousGrid = ThisGrid;
      
      /* create new grid */
      
      ThisGrid->GridData = new grid;
      
      /* set some the new grid's properties (rank, field types, etc.)
	 based on the current grid */
      
      ThisGrid->GridData->InheritProperties(SubgridParent->GridData);
      
      /* Set the new grid's positional parameters.
         (The zero indicates there are no particles (for now). */

      ThisGrid->GridData->PrepareGrid(SharedList[i].Rank,
				      SharedList[i].Dimension,
				      SharedList[i].LeftEdge,
				      SharedList[i].RightEdge,
				      0);

      ThisGrid->GridData->SetNumberOfParticles(SharedList[i].NumberOfParticles);
      if (SubgridParent != NULL)
	SubgridParent->GridData->SetNumberOfParticles(SubgridParent->GridData->ReturnNumberOfParticles() - SharedList[i].NumberOfParticles);

      ThisGrid->GridData->SetProcessorNumber(SubgridParent->GridData->ReturnProcessorNumber());

    }
  } // end: loop over shared grids

  /* CleanUp. */

  delete [] SendList;
  delete [] SharedList;

  return SUCCESS;
}

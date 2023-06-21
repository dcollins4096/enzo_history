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
/  COMMUNICATION ROUTINE: PARTITION GRID
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
#include <math.h>
#ifdef USE_MPI
#include <mpi.h>
#ifdef USE_MPE
#include <mpe.h>
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

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);


int CommunicationPartitionGrid(HierarchyEntry *Grid)
{

  if (NumberOfProcessors == 1)
    return SUCCESS;

  /* Declarations. */

  int Rank, dim, i, j, k, ijk, Dims[MAX_DIMENSION], Layout[] = {0,0,0};
  int TempDims[MAX_DIMENSION] /* , TempStart[MAX_DIMENSION] */ ;
  int *GridDims[MAX_DIMENSION], *StartIndex[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  /* Compute side length of root grid. */

  Grid->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
  for (dim = 0; dim < Rank; dim++)
    Dims[dim] -= 2*DEFAULT_GHOST_ZONES;
  float Edge = POW(float(Dims[0]*Dims[1]*Dims[2])/float(NumberOfProcessors), 
		   1/float(Rank));

  /* If using MPI, use their routine to calculate layout. */

#ifdef USE_MPI

  int LayoutTemp[] = {0,0,0};
  CHECK_MPI_ERROR(MPI_Dims_create(NumberOfProcessors, Rank, LayoutTemp));

  /* Swap layout because we want smallest value to be at Layout[0]. */

  for (dim = 0; dim < Rank; dim++)
    Layout[dim] = LayoutTemp[Rank-1-dim];

  /* Force some distributions if the default is brain-dead. */

  if (Rank == 3 && NumberOfProcessors == 64)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 4;

#endif /* USE_MPI */

  /* Generate arrays of grid dimensions and start positions. */

  int NumberOfNewGrids = 1;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {

    /* Compute number of new grids along this dimension. */

    if (Layout[dim] == 0)
      Layout[dim] = max(nint((float)(Dims[dim])/Edge), 1);

    GridDims[dim] = new int[Layout[dim]];
    StartIndex[dim] = new int[Layout[dim]];

    /* Compute dims and start indexes of the grids along this dim. */

    float ExactDims = float(Dims[dim])/float(Layout[dim]);
    float ExactCount = 0;
    int DisplacementCount = 0;
    for (i = 0; i < Layout[dim]; i++) {
      ExactCount += ExactDims;

      /* Compute the integer number of cells along this dimension
	 (if dim == 0 then make sure it is even as well since the FFT
	 requires this). */

      if (dim == 0)
	GridDims[dim][i] = nint(ExactCount*0.5)*2 - DisplacementCount;
      else
	GridDims[dim][i] = nint(ExactCount) - DisplacementCount;
      StartIndex[dim][i] = DisplacementCount;
      DisplacementCount += GridDims[dim][i];
    }
    NumberOfNewGrids *= Layout[dim];
  }

  if (MyProcessorNumber == ROOT_PROCESSOR)
  {
    printf("PartitionGrid (on all processors): Layout = %d %d %d\n", 
      Layout[0], Layout[1], Layout[2]);
    printf("NumberOfNewGrids = %d\n",NumberOfNewGrids);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
    {
      printf("GridDims[%d]: ",dim);
      for (i = 0; i < Layout[dim]; i++)
      {
        printf(" %d",GridDims[dim][i]);
      }
      printf("\n");
    }
    for (dim = 0; dim < MAX_DIMENSION; dim++) 
    {
      printf("StartIndex[%d]: ",dim);
      for (i = 0; i < Layout[dim]; i++)
      {
        printf(" %d",StartIndex[dim][i]);
      }
      printf("\n");
    }
  }

/*
  if ((ProblemType == 30) && (ParallelRootGridIO == 1) && (ParallelParticleIO == 1))
  {
    printf("Unigrid: %d\n", Unigrid);
    printf("Set Unigrid = 1\n");
    Unigrid = 1;
  }
*/

  /* Initialize the under subgrid field for particle movement. */

  printf("Call ZeroSUS on TopGrid\n");

  Grid->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);

  /* Generate this many grids (on this processor). */

/*
  Unigrid = 0;
  printf("Re-set Unigrid = 0\n");
*/

  printf("Grid structure: %d\n", (int) (sizeof(grid)));
  printf("SubGrids structure: %d\n", (int) ((Layout[0]*Layout[1]*Layout[2])*sizeof(grid)));

  grid *NewGrid, *OldGrid = Grid->GridData;
  grid **SubGrids = new grid*[Layout[0]*Layout[1]*Layout[2]];
  HierarchyEntry *ThisGrid;

  int gridcounter = 0;

  for (k = 0; k < Layout[2]; k++)
    for (j = 0; j < Layout[1]; j++)
      for (i = 0; i < Layout[0]; i++) {

	/* Allocate a new grid hierarchy entry and insert into linked list. */

	if (gridcounter == 0)
	  ThisGrid = Grid;
	else {
	  ThisGrid = new HierarchyEntry;
	  ThisGrid->NextGridThisLevel = Grid->NextGridThisLevel;
	  ThisGrid->NextGridNextLevel = NULL;
	  ThisGrid->ParentGrid        = Grid->ParentGrid;
	  Grid->NextGridThisLevel     = ThisGrid;
	}

	/* Allocate a new grid and prepare it. */

	NewGrid = new grid;
	ThisGrid->GridData = NewGrid;
	NewGrid->InheritProperties(OldGrid);
	NewGrid->SetGravityParameters(OldGrid->ReturnGravityBoundaryType());

	/* Compute grid region. */

      printf("GC K J I: %d %d %d %d\n",gridcounter,k,j,i);

	for (dim = 0; dim < MAX_DIMENSION; dim++) {
	  ijk = (dim == 0) ? i : ((dim == 1) ? j : k);
	  TempDims[dim] = GridDims[dim][ijk];
	  /* TempStart[dim] = StartIndex[dim][ijk]; */
	  LeftEdge[dim] = Left[dim] + (Right[dim] - Left[dim])*
	    FLOAT(StartIndex[dim][ijk])/FLOAT(Dims[dim]);
	  RightEdge[dim] = Left[dim] + (Right[dim] - Left[dim])*
	    FLOAT(StartIndex[dim][ijk]+TempDims[dim])/FLOAT(Dims[dim]);
	  if (dim < Rank)
	    TempDims[dim] += 2*DEFAULT_GHOST_ZONES;

        printf("  LeftEdge[%d] = %8.4f  RightEdge[%d] = %8.4f\n",
                 dim, LeftEdge[dim], dim, RightEdge[dim]);

	}

	NewGrid->PrepareGrid(Rank, TempDims, LeftEdge, RightEdge, 0);

	/* Record this subgrid number in the oldgrid's undersubgrid field. */

        printf("Call ZeroSUS on OldGrid with Value = %10.4e\n", float(gridcounter+1));

	if (OldGrid->ZeroSolutionUnderSubgrid(NewGrid, 
		   ZERO_UNDER_SUBGRID_FIELD, float(gridcounter+1)) == FAIL) {
	  fprintf(stderr, "Error in grid->ZeroSolutionUnderSubgrid.\n");
	  return FAIL;
	}
	SubGrids[gridcounter] = NewGrid;

	gridcounter++;

      }

    Unigrid = 0;
    printf("Re-set Unigrid = 0\n");

  /* Move Particles (while still on same processor). */

  if (OldGrid->MoveSubgridParticlesFast(gridcounter, SubGrids, TRUE) == FAIL) {
    fprintf(stderr, "Error in grid->MoveSubgridParticlesFast.\n");
    return FAIL;
  }
  delete [] SubGrids;

  /* Distribute new grids amoung processors (and copy out fields). */

  gridcounter = 0;
  ThisGrid = Grid;

//printf("Grid distribution\n");

  for (k = 0; k < Layout[2]; k++)
    for (j = 0; j < Layout[1]; j++)
      for (i = 0; i < Layout[0]; i++) {

	grid *NewGrid = ThisGrid->GridData;

	/* Broadcast the number of particles to the other processors
	   (OldGrid is assumed to be on the root processor). */

	if (NumberOfProcessors > 1) {
	  int IntTemp = NewGrid->ReturnNumberOfParticles();
//          printf("NewGrid->ReturnNumberOfParticles: %d\n", IntTemp);
	  CommunicationBroadcastValue(&IntTemp, ROOT_PROCESSOR);
	  NewGrid->SetNumberOfParticles(IntTemp);
	}

	/* Transfer from Old to New (which is still also on root processor) */

	FLOAT Zero[] = {0,0,0};
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  NewGrid->AllocateGrids();
	if (NewGrid->CopyZonesFromGrid(OldGrid, Zero) == FAIL) {
	  fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
	  return FAIL;
	}

	/* Set processor number of new grid.  Cyclic distribution. */

	int NewProc = gridcounter % NumberOfProcessors;

//        printf("K J I: %d %d %d  GC %d  NewProc %d\n",k,j,i,gridcounter,NewProc);
	/* Move Grid from current processor to new Processor. */

	NewGrid->CommunicationMoveGrid(NewProc);

	gridcounter++;
	ThisGrid = ThisGrid->NextGridThisLevel;

      }

  /* Clean up. */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Delete OldGrid\n");

  delete OldGrid;

  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("OldGrid deleted\n");

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    delete [] GridDims[dim];
    delete [] StartIndex[dim];
  }

  printf("Exit CommunicationPartitionGrid\n");
#ifdef USE_MPI
  JBPERF_START_MPI_BARRIER("MPI_Barrier");
  CHECK_MPI_ERROR(MPI_Barrier(MPI_COMM_WORLD));
  JBPERF_STOP_MPI_BARRIER("MPI_Barrier");
#endif

  return SUCCESS;
}

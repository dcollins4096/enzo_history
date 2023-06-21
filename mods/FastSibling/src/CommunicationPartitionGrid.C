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
#include "mpi.h"
#endif /* USE_MPI */
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
  float Edge = pow(float(Dims[0]*Dims[1]*Dims[2])/float(NumberOfProcessors), 
		   1/float(Rank));

  /* If using MPI, use their routine to calculate layout. */

#ifdef USE_MPI

  int LayoutTemp[] = {0,0,0};
  if (MPI_Dims_create(NumberOfProcessors, Rank, LayoutTemp) != MPI_SUCCESS) {
    fprintf(stderr, "Error in MPI_Dims_create.\n");
    return FAIL;
  }

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

      if (dim == 0 && SelfGravity == TRUE)
	GridDims[dim][i] = nint(ExactCount*0.5)*2 - DisplacementCount;
      else
	GridDims[dim][i] = nint(ExactCount) - DisplacementCount;
      StartIndex[dim][i] = DisplacementCount;
      DisplacementCount += GridDims[dim][i];
    }
    NumberOfNewGrids *= Layout[dim];
  }
  if (debug)
    printf("PartitionGrid: Layout = %d %d %d\n", 
	   Layout[0], Layout[1], Layout[2]);

  /* Initialize the under subgrid field for particle movement. */

  Grid->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);

  /* Generate this many grids (on this processor). */

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
	}

	NewGrid->PrepareGrid(Rank, TempDims, LeftEdge, RightEdge, 0);

	/* Record this subgrid number in the oldgrid's undersubgrid field. */

	if (OldGrid->ZeroSolutionUnderSubgrid(NewGrid, 
		   ZERO_UNDER_SUBGRID_FIELD, float(gridcounter+1)) == FAIL) {
	  fprintf(stderr, "Error in grid->ZeroSolutionUnderSubgrid.\n");
	  return FAIL;
	}
	SubGrids[gridcounter] = NewGrid;

	/* Loop over the list of subgrids attached to the original grid and
	   attach any which are within this new grid to the list children. */

	if (gridcounter > 0) {
	  HierarchyEntry *Temp = Grid->NextGridNextLevel;
	  HierarchyEntry **OldTemp = &(Grid->NextGridNextLevel);
	  int TempRank;
	  FLOAT TempLeft[MAX_DIMENSION], TempRight[MAX_DIMENSION];
	  while (Temp != NULL) {

	    /* Check if center of this subgrid lies within the current new
	       patch. */

	    Temp->GridData->ReturnGridInfo(&TempRank, TempDims, 
					   TempLeft, TempRight);
	    int CopyEntry = TRUE;
	    for (dim = 0; dim < TempRank; dim++)
	      if (0.5*(TempLeft[dim]+TempRight[dim]) < LeftEdge[dim] ||
		  0.5*(TempLeft[dim]+TempRight[dim]) >= RightEdge[dim])
		CopyEntry = FALSE;

	    /* If it does, snip this hierarchy entry out of the list of the
	       first patch and move it to the head of the linked list of
	       subgrids for this patch. */

	    if (CopyEntry) {
	      *OldTemp = Temp->NextGridThisLevel;
	      Temp->NextGridThisLevel = ThisGrid->NextGridNextLevel;
	      ThisGrid->NextGridNextLevel = Temp;
	      Temp->ParentGrid = ThisGrid;
	      Temp = *OldTemp;  /* follow the linked list onwards. */
	    } else {

	      /* Otherwise, follow the linked list to the next grid, recording
		 a pointer to the link in case we need to update it. */

	      OldTemp = &(Temp->NextGridThisLevel);
	      Temp = Temp->NextGridThisLevel;
	    }

	  } // end: loop over linked list
	} // end: if (gridcounter > 0)

	gridcounter++;

      } // end: loop over the new patches
  
  /* Move Particles (while still on same processor). */

  if (OldGrid->MoveSubgridParticlesFast(gridcounter, SubGrids, TRUE) == FAIL) {
    fprintf(stderr, "Error in grid->MoveSubgridParticlesFast.\n");
    return FAIL;
  }
  delete [] SubGrids;

  /* Distribute new grids amoung processors (and copy out fields). */

  gridcounter = 0;
  ThisGrid = Grid;
  for (k = 0; k < Layout[2]; k++)
    for (j = 0; j < Layout[1]; j++)
      for (i = 0; i < Layout[0]; i++) {

	NewGrid = ThisGrid->GridData;

	/* Broadcast the number of particles to the other processors
	   (OldGrid is assumed to be on the root processor). */

	if (NumberOfProcessors > 1) {
	  int IntTemp = NewGrid->ReturnNumberOfParticles();
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

	/*
	if (MyProcessorNumber == NewGrid->ReturnProcessorNumber() ||
	    MyProcessorNumber == OldGrid->ReturnProcessorNumber()  )
	  if (OldGrid->CommunicationSendRegion(NewGrid, NewProc, ALL_FIELDS, 
				    NEW_ONLY, TempStart, TempDims) == FAIL) {
	    fprintf(stderr, "Error in grid->CommunicationSendRegion.\n");
	    return FAIL;
	    } */

	/* Set processor number of new grid.  Cyclic distribution. */

	int NewProc = gridcounter % NumberOfProcessors;

	/* Move Grid from current processor to new Processor. */

	NewGrid->CommunicationMoveGrid(NewProc);

	gridcounter++;
	ThisGrid = ThisGrid->NextGridThisLevel;

      }

  /* Clean up. */

  delete OldGrid;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    delete [] GridDims[dim];
    delete [] StartIndex[dim];
  }

  return SUCCESS;
}

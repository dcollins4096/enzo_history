/***********************************************************************
/
/  GRID CLASS (COALESCE NEIGHBOURING SUBGRIDS)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

void grid::CoalesceSubgrids(GridList &SubgridList)
{
  /* declarations */

  int CombinedEndIndex[MAX_DIMENSION], CombinedStartIndex[MAX_DIMENSION];
  int dim, dim2, size, CheckSubgridDim[MAX_DIMENSION], 
      SubgridDim[MAX_DIMENSION], Dims[MAX_DIMENSION];
  float AfterRatio, BeforeRatio, NarrowSubgridFactor;
  GridList *Subgrid, *CheckSubgrid, *PreviousSubgrid;

  /* Traverse subgrid list until no further changes are made. */

  int ChangesMade = TRUE;
  while (ChangesMade) {
    ChangesMade = FALSE;

    /* loop over the list of subgrids */

    Subgrid = SubgridList.NextGrid;
    while (Subgrid != NULL) {

      /* loop over list of remaining subgrids */
	
      CheckSubgrid    = Subgrid->NextGrid;
      PreviousSubgrid = Subgrid;
      while (CheckSubgrid != NULL) {

	/* calculate efficiency of combined subgrid */

	int CombinedSize = 1, SubgridSize = 1, CheckSubgridSize = 1;
	for (dim = 0; dim < GridRank; dim++) {

	  /* Compute start and end indexes of combined grid. */

	  CombinedStartIndex[dim] = min(CheckSubgrid->StartIndex[dim],
			  		     Subgrid->StartIndex[dim]);
	  CombinedEndIndex[dim] = max(CheckSubgrid->EndIndex[dim],
				           Subgrid->EndIndex[dim]);

	  /* Compute active dimension of all three subgrids. */

	  Dims[dim]            = CombinedEndIndex[dim]         - 
                                 CombinedStartIndex[dim]       + 1;
	  SubgridDim[dim]      =      Subgrid->EndIndex[dim]   - 
	                              Subgrid->StartIndex[dim] + 1;
	  CheckSubgridDim[dim] = CheckSubgrid->EndIndex[dim]   - 
	                         CheckSubgrid->StartIndex[dim] + 1;

	  /* Compute active sizes of three subgrids. */

	  CombinedSize     *= Dims[dim];
	  SubgridSize      *= SubgridDim[dim];
	  CheckSubgridSize *= CheckSubgridDim[dim];

	}

	/* Compute longest/shortest length ratios. */

	BeforeRatio = 1.0, AfterRatio = 1.0;
	for (dim = 0; dim < GridRank; dim++)
	  for (dim2 = 0; dim2 < GridRank; dim2++) {
	    BeforeRatio = max(BeforeRatio, float(CheckSubgridDim[dim]) /
			                   float(CheckSubgridDim[dim2]) );
	    BeforeRatio = max(BeforeRatio, float(SubgridDim[dim])      /
			                   float(SubgridDim[dim2])      );
	    AfterRatio  = max(AfterRatio,  float(Dims[dim])            / 
			                   float(Dims[dim2])            );
	  }
			   
	/* Try to eliminate narrow subgrids (< 3 zones) by making the
	   acceptance criterion lower when combining such subgrids. */

	NarrowSubgridFactor = 1.0;
	for (dim = 0; dim < GridRank; dim++) {

	  if ((SubgridDim[dim] < 2 || CheckSubgridDim[dim] < 2) && 
	      Dims[dim] >= 2)
	    NarrowSubgridFactor *= 0.70;
	  else if ((SubgridDim[dim] < 3 || CheckSubgridDim[dim] < 3) && 
		   Dims[dim] >= 3)
	    NarrowSubgridFactor *= 0.80;
//	  else if ((SubgridDim[dim] < 4 || CheckSubgridDim[dim] < 4) && 
//		   Dims[dim] >= 4)
//	    NarrowSubgridFactor *= 0.90;

//	    NarrowSubgridFactor = min(NarrowSubgridFactor, 0.90);
	  
	}

	if (AfterRatio > BeforeRatio)
	  NarrowSubgridFactor = 1.0;

	float CombinedEfficiency = float(CheckSubgrid->NumberFlagged + 
				         Subgrid->NumberFlagged) / 
				   float(CombinedSize);
/*      float SubgridEfficiency  = float(Subgrid->NumberFlagged) / 
	                         float(SubgridSize);
	float CheckSubgridEfficiency = float(CheckSubgrid->NumberFlagged) /
                                     float(CheckSubgridSize); */
	float AverageEfficiency = float(CheckSubgrid->NumberFlagged +
				        Subgrid->NumberFlagged)/
				  float(SubgridSize + CheckSubgridSize);

	/* check to see if the new grid is good enough */

	if (CombinedEfficiency > MinimumEfficiency*NarrowSubgridFactor ||
	    (CombinedEfficiency >= AverageEfficiency)                    ) {

	  /* combine grids */

	  for (dim = 0; dim < GridRank; dim++) {
	    Subgrid->GridDimension[dim] = Dims[dim];
	    Subgrid->StartIndex[dim]    = CombinedStartIndex[dim];
	    Subgrid->EndIndex[dim]      = CombinedEndIndex[dim];
	  }
	  Subgrid->NumberFlagged       += CheckSubgrid->NumberFlagged;
	  
	  /* remove CheckSubgrid */

	  PreviousSubgrid->NextGrid = CheckSubgrid->NextGrid;
	  delete CheckSubgrid;
	  CheckSubgrid = PreviousSubgrid->NextGrid;

	  ChangesMade = TRUE;

	}
	else {
	  PreviousSubgrid = CheckSubgrid;
	  CheckSubgrid = CheckSubgrid->NextGrid;
	}

      } // end of CheckSubgrid loop

      Subgrid = Subgrid->NextGrid;

    } // end of main loop over subgrids

  } // end: while ChangesMade

  /* loop over the list of subgrids and fill out the rest of the information
     that will be needed later (left/right edges and dimensions). */

  Subgrid = SubgridList.NextGrid;
  while (Subgrid != NULL) {

    size = 1;

#define NO_SQUARE_GRIDS

#ifdef SQUARE_GRIDS
    int MaxSize = 0;
    for (dim = 0; dim < GridRank; dim++)
      if (Subgrid->GridDimension[dim] > MaxSize)
	MaxSize = Subgrid->GridDimension[dim];
    for (dim = 0; dim < GridRank; dim++) {
      int DimDiff = int(Subgrid->GridDimension[dim] - MaxSize + 1)/2;
      Subgrid->StartIndex[dim] -= DimDiff;
      Subgrid->GridDimension[dim] = MaxSize;
      Subgrid->EndIndex[dim] = Subgrid->StartIndex[dim] + MaxSize - 1;
    }
#endif /* SQUARE_GRIDS */

    for (dim = 0; dim < GridRank; dim++) {

      size *= Subgrid->GridDimension[dim];

      /* compute Left/Right Edges */
      
      Subgrid->GridLeftEdge[dim] = GridLeftEdge[dim] + 
	FLOAT(Subgrid->StartIndex[dim] - GridStartIndex[dim])*
	  (*CellWidth[dim]);
      Subgrid->GridRightEdge[dim] = GridLeftEdge[dim] + 
	FLOAT(Subgrid->EndIndex[dim] + 1 - GridStartIndex[dim])*
	  (*CellWidth[dim]);

      /* add boundary zones to dimension */

      if (GridDimension[dim] > 1)
	Subgrid->GridDimension[dim] = (Subgrid->GridDimension[dim])*RefineBy +
	                               2*DEFAULT_GHOST_ZONES;

    } // end loop over dimensions

    /* Compute efficiency & report. */

    if (debug)
      printf("CoalesceSubgrids: efficiency = %6.1f%% (%d/%d)\n", 
	     float(Subgrid->NumberFlagged)/float(size)*100.0,
	     Subgrid->NumberFlagged, size);

    /* Next Subgrid. */
    
    Subgrid = Subgrid->NextGrid;

  } // end loop over subgrids
  
}

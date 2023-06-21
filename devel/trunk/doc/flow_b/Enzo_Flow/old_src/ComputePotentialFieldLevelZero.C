/***********************************************************************
/
/  COMPUTE THE POTENTIAL FIELD
/
/  written by: Greg Bryan
/  date:       January, 1998
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

int CommunicationParallelFFT(region *InRegion, int NumberOfInRegions,
			     region **OutRegion, int *NumberOfOutRegions,
			     int DomainDim[], int Rank,
			     int direction, int TransposeOnCompletion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids)
{

  /* Static declarations (for Green's function). */

  static int FirstCall = TRUE, NumberOfGreensRegions;
  static region *GreensRegion;

  /* Declarations. */

  region *OutRegion = NULL;
  int NumberOfOutRegions, DomainDim[MAX_DIMENSION];
  int i, j, n, grid, grid2;

  /* Allocate space for grid info. */

  int NumberOfRegions = NumberOfGrids;
  region *InitialRegion = new region[NumberOfRegions];

  /* Compute adot/a at time = t+1/2dt (time-centered). */

  FLOAT a = 1, dadt, MidTime = Grids[0]->GridData->ReturnTime() + 
                           0.5*Grids[0]->GridData->ReturnTimeStep();
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(MidTime, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }

  /* ------------------------------------------------------------------- */
  /* If this is the first time this routine has been called, then generate
     the Green's function. */

  if (FirstCall) {

    if (MetaData->GravityBoundary == TopGridPeriodic) {

      /* Periodic -- Prepare in k-space. */

      NumberOfGreensRegions = NumberOfGrids;
      GreensRegion = new region[NumberOfGreensRegions];
      for (grid = 0; grid < NumberOfGrids; grid++)
	if (Grids[grid]->GridData->PreparePeriodicGreensFunction(
					     &(GreensRegion[grid])) == FAIL) {
	  fprintf(stderr, "Error in grid->PreparePeriodicGreensFunction.\n");
	  return FAIL;
	}

    } else {

      fprintf(stderr, "Isolated BC's not yet implemented.\n");
      return FAIL;

#ifdef UNUSED

      for (grid = 0; grid < NumberOfGrids; grid++) {

	/* Generate Green's function in real space (doesn't work!). */
	
	if (Grids[grid]->GridData->PrepareGreensFunction() == FAIL) {
	  fprintf(stderr, "Error in grid->PrepareGreensFunction.\n");
	  return FAIL;
	}
      
	/* Turn it into regions. */
	
	if (Grids[grid]->GridData->PrepareFFT(&InitialRegion[grid], 
					      POTENTIAL_FIELD, DomainDim) 
	    == FAIL) {
	  fprintf(stderr, "Error in grid->PrepareFFT.\n");
	  return FAIL;
	}

      } // end loop over grids

      /* Forward FFT Green's function. */

      if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
				   &GreensRegion, &NumberOfGreensRegions,
				   DomainDim, MetaData->TopGridRank,
				   FFT_FORWARD, FALSE) == FAIL) {
	fprintf(stderr, "Error in CommunicationParallelFFT.\n");
	return FAIL;
      }

#endif /* UNUSED */

    } // end: if (Periodic)

    FirstCall = FALSE;

  } // end: if (FirstCall)

  /* ------------------------------------------------------------------- */
  /* Generate FFT regions for density field. */

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->PrepareFFT(&InitialRegion[grid], 
					  GRAVITATING_MASS_FIELD, DomainDim) 
	== FAIL) {
      fprintf(stderr, "Error in grid->PrepareFFT.\n");
      return FAIL;
    }

  /* Forward FFT density field. */

  if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
			       &OutRegion, &NumberOfOutRegions,
			       DomainDim, MetaData->TopGridRank,
			       FFT_FORWARD, TRUE) == FAIL) {
    fprintf(stderr, "Error in CommunicationParallelFFT.\n");
    return FAIL;
  }

  /* Quick error check. */

  if (NumberOfOutRegions != NumberOfGreensRegions) {
    fprintf(stderr, "OutRegion(%d) != GreensRegion(%d)\n", NumberOfOutRegions,
	    NumberOfGreensRegions);
    return FAIL;
  }

  /* Compute coefficient for Greens function. */

  float coef = GravitationalConstant/a;
  //  for (int dim = 0; dim < MetaData->TopGridRank; dim++)
  //    coef *= (DomainRightEdge[dim] - DomainLeftEdge[dim])/float(DomainDim[dim]);
			       
  /* Multiply density by Green's function to get potential. */

  for (i = 0; i < NumberOfGreensRegions; i++)
    if (OutRegion[i].Data != NULL) {
      int size = OutRegion[i].RegionDim[0]*OutRegion[i].RegionDim[1]*
	         OutRegion[i].RegionDim[2];
      for (n = 0, j = 0; j < size; j += 2, n++) {
	OutRegion[i].Data[j  ] *= coef*GreensRegion[i].Data[n];
	OutRegion[i].Data[j+1] *= coef*GreensRegion[i].Data[n];
      }
    }

  /* Inverse FFT potential field. */

  if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
			       &OutRegion, &NumberOfOutRegions,
			       DomainDim, MetaData->TopGridRank,
			       FFT_INVERSE, TRUE) == FAIL) {
    fprintf(stderr, "Error in CommunicationParallelFFT.\n");
    return FAIL;
  }

  /* Copy Potential in active region into while grid. */

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->FinishFFT(&InitialRegion[grid], POTENTIAL_FIELD,
			       DomainDim) == FAIL) {
      fprintf(stderr, "Error in grid->FinishFFT.\n");
      return FAIL;
    }

  /* Update boundary regions of potential
     (first set BCTempL/R which are fluid BC's because that's the format
      that CheckForOverlap takes). */

  boundary_type BCTempLeft[MAX_DIMENSION], BCTempRight[MAX_DIMENSION];
  if (Grids[0]->GridData->ReturnGravityBoundaryType() == TopGridPeriodic) {
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      BCTempLeft[dim] = BCTempRight[dim] = periodic;
  } else {
    fprintf(stderr, "Only periodic gravity BC's allowed.\n");
    return FAIL;
  }

  for (grid = 0; grid < NumberOfGrids; grid++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				      BCTempLeft, BCTempRight,
     	                              &grid::CopyPotentialField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	return FAIL;
      }

  /* Clean up. */

  delete [] InitialRegion;
  if (OutRegion != InitialRegion)
    delete [] OutRegion;

  return SUCCESS;
}

/***********************************************************************
/
/  EVOLVE LEVEL ROUTINES (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modified1:
/
/  PURPOSE:  This is a collection of routines called by EvolveLevel.
/            These have been optimized for enhanced message passing 
/            performance by performing two passes -- one which generates
/            sends and the second which receives them.
/
************************************************************************/

#include <stdio.h>
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

/* function prototypes */

int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
				LevelHierarchyEntry *LevelArray[], int level);
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);



extern int CopyPotentialFieldAverage;

/* ======================================================================= */
/* This routine sets all the boundary conditions for Grids by either
   interpolating from their parents or copying from sibling grids. */

int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior)
{
  int grid, grid2;

  /* -------------- FIRST PASS ----------------- */

  CommunicationDirection = COMMUNICATION_SEND;

  for (grid = 0; grid < NumberOfGrids; grid++) {

    /* a) Interpolate boundaries from the parent grid or set external
       boundary conditions. */

    if (level == 0) {
      if (Grids[grid]->GridData->SetExternalBoundaryValues(Exterior) 
	  == FAIL) {
	fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
	return FAIL;
      }
    }
    else {
      if ((Grids[grid]->GridData->InterpolateBoundaryFromParent
	   (Grids[grid]->ParentGrid->GridData)) == FAIL) {
	fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
	return FAIL;
      }
    }
    
    /* b) Copy any overlapping zones for sibling grids.  */

    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid)
	== FAIL) {
      fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
    }

  } // end loop over grids

  /* -------------- SECOND PASS ----------------- */

  CommunicationDirection = COMMUNICATION_RECEIVE;

  for (grid = 0; grid < NumberOfGrids; grid++) {

    /* a) Interpolate boundaries from the parent grid or set external
       boundary conditions. */

    if (level > 0)
      if ((Grids[grid]->GridData->InterpolateBoundaryFromParent
	   (Grids[grid]->ParentGrid->GridData)) == FAIL) {
	fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
	return FAIL;
      }
    
    /* b) Copy any overlapping zones for sibling grids.  */

    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid)
	== FAIL) {
      fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
    }

  } // end loop over grids

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  return SUCCESS;
}



/* ======================================================================= */
/* This routine prepares the density field for all the grids on this level,
   both particle and baryonic densities.  It also calculates the potential
   field if this is level 0 (since this involves communication). */

int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			int level, TopGridData *MetaData)
{


  //=================================================================

  JBPERF_INIT;
  JBPERF_START_LEVEL("PDF00"); // PrepareDensityField

  int grid, grid2;

  /* Set the time for evaluation of the fields, etc. */

  FLOAT EvaluateTime = LevelArray[level]->GridData->ReturnTime() +
                   0.5*LevelArray[level]->GridData->ReturnTimeStep();

  /* If level is above MaximumGravityRefinementLevel, then just update the
     gravity at the MaximumGravityRefinementLevel. */

  int reallevel = level;
  level = min(level, MaximumGravityRefinementLevel);

  /* Create an array (Grids) of all the grids. */

  JBPERF_START_LEVEL("PDF01"); // create Grids
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  JBPERF_STOP_LEVEL("PDF01");

  /* Grids: Deposit particles in their GravitatingMassFieldParticles. */

  JBPERF_START_LEVEL("PDF02"); // deposit particles (send)
  CommunicationDirection = COMMUNICATION_SEND;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (DepositParticleMassField(Grids[grid], EvaluateTime) == FAIL) {
      fprintf(stderr, "Error in DepositParticleMassField.\n");
      return FAIL;
    }
  JBPERF_STOP_LEVEL("PDF02");

  JBPERF_START_LEVEL("PDF03"); // deposit particles (recv)
  CommunicationDirection = COMMUNICATION_RECEIVE;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (DepositParticleMassField(Grids[grid], EvaluateTime) == FAIL) {
      fprintf(stderr, "Error in DepositParticleMassField.\n");
      return FAIL;
    }
  JBPERF_STOP_LEVEL("PDF03");

  /* Grids: compute the GravitatingMassField (baryons & particles). */

  JBPERF_START_LEVEL("PDF04"); //GravitatingMassField (send)
  CommunicationDirection = COMMUNICATION_SEND;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (PrepareGravitatingMassField(Grids[grid], MetaData, LevelArray,
				    level) == FAIL) {
      fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
      return FAIL;
    }
  JBPERF_STOP_LEVEL("PDF04");

  JBPERF_START_LEVEL("PDF05"); // GravitatingMassField (recv)
  CommunicationDirection = COMMUNICATION_RECEIVE;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (PrepareGravitatingMassField(Grids[grid], MetaData, LevelArray,
				    level) == FAIL) {
      fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
      return FAIL;
    }
  JBPERF_STOP_LEVEL("PDF05");

  /* Copy overlapping mass fields to ensure consistency and B.C.'s. */

  JBPERF_START_LEVEL("PDF06"); // copy overlap (send)
  CommunicationDirection = COMMUNICATION_SEND;
  for (grid = 0; grid < NumberOfGrids; grid++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyOverlappingMassField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
	return FAIL;
      }
  JBPERF_STOP_LEVEL("PDF06");

  JBPERF_START_LEVEL("PDF07"); // copy overlap (recv)
  CommunicationDirection = COMMUNICATION_RECEIVE;
  for (grid = 0; grid < NumberOfGrids; grid++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyOverlappingMassField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
	return FAIL;
      }
  JBPERF_STOP_LEVEL("PDF07");

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  /* Compute the potential for the top grid. */

  JBPERF_START_LEVEL("PDF08"); // top grid potential
  if (level == 0)
    if (ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids)
	== FAIL) {
      fprintf(stderr, "Error in ComputePotentialFieldLevelZero.\n");
      return FAIL;
    }
  JBPERF_STOP_LEVEL("PDF08");

  /* Compute a first iteration of the potential and share BV's. */

#define ITERATE_POTENTIAL
#ifdef ITERATE_POTENTIAL
      if (level > 0) {
	CopyPotentialFieldAverage = 1;
	for (int iterate = 0; iterate < 1; iterate++) {

	  if (iterate > 0)
	    CopyPotentialFieldAverage = 2;

  JBPERF_START_LEVEL("PDF09"); // mid-grid potential
	  int Dummy, grid2;
	  for (grid = 0; grid < NumberOfGrids; grid++)
	    if (Grids[grid]->GridData->SolveForPotential(Dummy, level,
							 EvaluateTime) 
		== FAIL) {
	      fprintf(stderr, "Error in grid->SolveForPotential.\n");
	      return FAIL;
	    }
  JBPERF_STOP_LEVEL("PDF09");

  JBPERF_START_LEVEL("PDF10"); // overlap check (send)
	  CommunicationDirection = COMMUNICATION_SEND;
	  for (grid = 0; grid < NumberOfGrids; grid++)
	    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	     if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyPotentialField) == FAIL) {
	       fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	       return FAIL;
	     }
  JBPERF_STOP_LEVEL("PDF10");

  JBPERF_START_LEVEL("PDF11"); // overlap check (recv)

	  CommunicationDirection = COMMUNICATION_RECEIVE;
	  for (grid = 0; grid < NumberOfGrids; grid++)
	    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	     if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyPotentialField) == FAIL) {
	       fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	       return FAIL;
	     }

  JBPERF_STOP_LEVEL("PDF11");

	  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;


	}
	CopyPotentialFieldAverage = 0;
      }
#endif /* ITERATE_POTENTIAL */

  /* if level > MaximumGravityRefinementLevel, then do final potential
     solve (and acceleration interpolation) here rather than in the main
     EvolveLevel since it involves communications. */

  if (reallevel > MaximumGravityRefinementLevel) {

    /* compute potential and acceleration on coarser level [LOCAL]
       (but only if there is at least a subgrid -- it should be only
        if there is a subgrrid on reallevel, but this is ok). */

  JBPERF_START_LEVEL("PDF12"); // coarse-grid potential
    for (grid = 0; grid < NumberOfGrids; grid++) 
      if (Grids[grid]->NextGridNextLevel != NULL) {
	Grids[grid]->GridData->SolveForPotential(level,
					       MaximumGravityRefinementLevel);
	Grids[grid]->GridData->ComputeAccelerationField(
                           (HydroMethod == Zeus_Hydro) ? ZEUS_GRIDS : GRIDS,
					       MaximumGravityRefinementLevel);
      }
  JBPERF_STOP_LEVEL("PDF12");

    /* Interpolate potential for reallevel grids from coarser grids. */

    int Dummy;
    LevelHierarchyEntry *Temp = LevelArray[reallevel];

  JBPERF_START_LEVEL("PDF13"); // interpolate potential (send)
    CommunicationDirection = COMMUNICATION_SEND;
    while (Temp != NULL) {
      HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
      for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	Temp3 = Temp3->ParentGrid;
      if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->InterpolateAccelerations.\n");
	return FAIL;
      }
      Temp = Temp->NextGridThisLevel;
    }
  JBPERF_STOP_LEVEL("PDF13");

  JBPERF_START_LEVEL("PDF14"); // interpolate potential (recv)
    CommunicationDirection = COMMUNICATION_RECEIVE;
    Temp = LevelArray[reallevel];
    while (Temp != NULL) {
      HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
      for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	Temp3 = Temp3->ParentGrid;
      if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->InterpolateAccelerations.\n");
	return FAIL;
      }
      Temp = Temp->NextGridThisLevel;
    }
  JBPERF_STOP_LEVEL("PDF14");

    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  } // end: if (reallevel > MaximumGravityRefinementLevel)

  JBPERF_STOP_LEVEL("PDF00")

  return SUCCESS;
}


/* ======================================================================= */
/* This routines does the flux correction and project for all grids on this
   level from the list of subgrids. */

int UpdateFromFinerGrids(HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[], 
			 fluxes **SubgridFluxesEstimate[])

{

  int grid, subgrid;
  HierarchyEntry *NextGrid;

  /* Define a temporary flux holder for the refined fluxes. */

  fluxes SubgridFluxesRefined;

  /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */

  /* -------------- FIRST PASS ----------------- */

  CommunicationDirection = COMMUNICATION_SEND;

  for (grid = 0; grid < NumberOfGrids; grid++) {

    /* Loop over subgrids for this grid. */

    NextGrid = Grids[grid]->NextGridNextLevel;
    subgrid = 0;
    while (NextGrid != NULL && FluxCorrection) {

      /* Project subgrid's refined fluxes to the level of this grid. */

      if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		      Grids[grid]->GridData, SubgridFluxesRefined) == FAIL) {
	fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	return FAIL;
      }

      NextGrid = NextGrid->NextGridThisLevel;
      subgrid++;
    }

    /* Loop over subgrids for this grid: replace solution. */

    NextGrid = Grids[grid]->NextGridNextLevel;
    while (NextGrid != NULL) {

      /* Project the subgrid solution into this grid. */

      if (NextGrid->GridData->ProjectSolutionToParentGrid
	                                   (*Grids[grid]->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	return FAIL;
      }

      NextGrid = NextGrid->NextGridThisLevel;
    }

  } // end of loop over subgrids

  /* -------------- SECOND PASS ----------------- */

  CommunicationDirection = COMMUNICATION_RECEIVE;

  for (grid = 0; grid < NumberOfGrids; grid++) {

    /* Loop over subgrids for this grid. */

    NextGrid = Grids[grid]->NextGridNextLevel;
    subgrid = 0;
    while (NextGrid != NULL && FluxCorrection) {

      /* Project subgrid's refined fluxes to the level of this grid. */

      if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		      Grids[grid]->GridData, SubgridFluxesRefined) == FAIL) {
	fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	return FAIL;
      }
	
      /* Correct this grid for the refined fluxes (step #19)
	 (this also deletes the fields in SubgridFluxesRefined). */
      
      if (Grids[grid]->GridData->CorrectForRefinedFluxes
	  (SubgridFluxesEstimate[grid][subgrid], &SubgridFluxesRefined, 
	   SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1]     )
	  == FAIL) {
	fprintf(stderr, "Error in grid->CorrectForRefinedFluxes.\n");
	return FAIL;
      }

      NextGrid = NextGrid->NextGridThisLevel;
      subgrid++;
    }

    /* Loop over subgrids for this grid: replace solution. */

    NextGrid = Grids[grid]->NextGridNextLevel;
    while (NextGrid != NULL) {

      /* Project the subgrid solution into this grid. */

      if (NextGrid->GridData->ProjectSolutionToParentGrid
	                                   (*Grids[grid]->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	return FAIL;
      }

      NextGrid = NextGrid->NextGridThisLevel;
    }

  } // end of loop over subgrids

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  return SUCCESS;
}



/* ======================================================================= */
/* This routine simply converts a linked list of grids into an array of
   pointers. */

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[])
{

  /* Count the number of grids on this level. */

  int NumberOfGrids = 0, counter = 0;
  LevelHierarchyEntry *Temp = LevelArray[level];
  while (Temp != NULL) {
    NumberOfGrids++;
    Temp             = Temp->NextGridThisLevel;
  }

  /* Create a list of pointers and number of subgrids (and fill it out). */

  *Grids = new HierarchyEntry *[NumberOfGrids];
  Temp = LevelArray[level];
  while (Temp != NULL) {
    (*Grids)[counter++] = Temp->GridHierarchyEntry;
    Temp              = Temp->NextGridThisLevel;
  }

  return NumberOfGrids;
}
  

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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "communication.h"

/* function prototypes */

int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
int PrepareGravitatingMassField1(HierarchyEntry *Grid);
int PrepareGravitatingMassField2(HierarchyEntry *Grid, int grid1, 
			      SiblingGridList SiblingList[],
		              TopGridData *MetaData, int level);
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
			       HierarchyEntry *Grids[], int NumberOfGrids);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL);
double ReturnWallTime();


extern int CopyPotentialFieldAverage;

#define GRIDS_PER_LOOP 20000

/* ======================================================================= */
/* This routine sets all the boundary conditions for Grids by either
   interpolating from their parents or copying from sibling grids. */

int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior)
{
  int grid1, grid2, StartGrid, EndGrid;

  /* Do a batch of grids at a time; this is a loop over the batches. */

  double start_time = ReturnWallTime();
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);
    
    /* -------------- FIRST PASS ----------------- */
    /* Here, we just generate the calls to generate the receive buffers,
       without actually doing anything. */

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* a) Interpolate boundaries from the parent grid or set external
	 boundary conditions. */

      CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      if (level == 0) {
	if (Grids[grid1]->GridData->SetExternalBoundaryValues(Exterior) 
	    == FAIL) {
	  fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
	  return FAIL;
	}
      }
      else {
	if ((Grids[grid1]->GridData->InterpolateBoundaryFromParent
	     (Grids[grid1]->ParentGrid->GridData)) == FAIL) {
	  fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
	  return FAIL;
	}
      }

    } // end loop over grids

    /* -------------- SECOND PASS ----------------- */
    /* Now we generate all the sends, and do all the computation for grids
       which are on the same processor as well. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* a) Interpolate boundaries from the parent grid or set external
	 boundary conditions. */

      if (level > 0)
	if ((Grids[grid1]->GridData->InterpolateBoundaryFromParent
	     (Grids[grid1]->ParentGrid->GridData)) == FAIL) {
	  fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
	  return FAIL;
	}

    } // end loop over this batch of grids

    /* -------------- THIRD PASS ----------------- */
    /* In this final step, we get the messages as they come in and then
       match them to the methods which generate the receive handle. */

    if (CommunicationReceiveHandler() == FAIL)
      return FAIL;

  } // end loop over batchs of grids

  PerformanceTimers[4] += ReturnWallTime() - start_time;
  start_time = ReturnWallTime();
  
  /* Do a batch of grids at a time; this is a loop over the batches. */

  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);
    
    /* -------------- FIRST PASS ----------------- */
    /* b) Copy any overlapping zones for sibling grids.  */

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	if (Grids[grid1]->GridData->CheckForOverlap(
				     SiblingList[grid1].GridList[grid2],
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid)
	    == FAIL) {
	  fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
	}
    } // end loop over grids

    /* -------------- SECOND PASS ----------------- */
    /* b) Copy any overlapping zones for sibling grids.  */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	if (Grids[grid1]->GridData->CheckForOverlap(
				     SiblingList[grid1].GridList[grid2],
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid)
	    == FAIL) {
	  fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
	}
    } // end loop over this batch of grids

    /* -------------- THIRD PASS ----------------- */

    if (CommunicationReceiveHandler() == FAIL)
      return FAIL;

  } // end loop over batchs of grids

  PerformanceTimers[17] += ReturnWallTime() - start_time;
  return SUCCESS;
}



/* ======================================================================= */
/* This routine prepares the density field for all the grids on this level,
   both particle and baryonic densities.  It also calculates the potential
   field if this is level 0 (since this involves communication). */

int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData)
{

  int grid1, grid2, StartGrid, EndGrid;
  double start_time = ReturnWallTime();

  /* Set the time for evaluation of the fields, etc. */

  FLOAT EvaluateTime = LevelArray[level]->GridData->ReturnTime() +
                   0.5*LevelArray[level]->GridData->ReturnTimeStep();

  /* If level is above MaximumGravityRefinementLevel, then just update the
     gravity at the MaximumGravityRefinementLevel. */

  int reallevel = level;
  level = min(level, MaximumGravityRefinementLevel);

  /* Create an array (Grids) of all the grids. */

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

  /* Grids: Deposit particles in their GravitatingMassFieldParticles.
     (Do a batch of grids at a time; this is a loop over the batches) */

  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    //    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      if (DepositParticleMassField(Grids[grid1], EvaluateTime) == FAIL) {
	fprintf(stderr, "Error in DepositParticleMassField.\n");
	return FAIL;
      }

    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      if (DepositParticleMassField(Grids[grid1], EvaluateTime) == FAIL) {
	fprintf(stderr, "Error in DepositParticleMassField.\n");
	return FAIL;
      }

    /* Finally, receive the data and process it. */
    
    if (CommunicationReceiveHandler() == FAIL)
      return FAIL;

  } // end loop over batches of grids

  PerformanceTimers[18] += ReturnWallTime() - start_time;
  double time1 = ReturnWallTime();

  /* Grids: compute the GravitatingMassField (baryons & particles). */
  /*   This is now split into two section. */

  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 1 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    //    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      if (PrepareGravitatingMassField1(Grids[grid1]) == FAIL) {
	fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
	return FAIL;
      }

    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      if (PrepareGravitatingMassField1(Grids[grid1]) == FAIL) {
	fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
	return FAIL;
      }

    /* Finally, receive the data and process it. */
    
    if (CommunicationReceiveHandler() == FAIL)
      return FAIL;

  } // end loop over batches of grids

  PerformanceTimers[19] += ReturnWallTime() - time1;
  time1 = ReturnWallTime();

  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 2 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    //    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      if (PrepareGravitatingMassField2(Grids[grid1], grid1, SiblingList,
				       MetaData, level) == FAIL) {
	fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
	return FAIL;
      }

    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      if (PrepareGravitatingMassField2(Grids[grid1], grid1, SiblingList,
				       MetaData, level) == FAIL) {
	fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
	return FAIL;
      }

    /* Finally, receive the data and process it. */
    
    if (CommunicationReceiveHandler() == FAIL)
      return FAIL;

  } // end loop over batches of grids

  PerformanceTimers[20] += ReturnWallTime() - time1;
  time1 = ReturnWallTime();

  /* Copy overlapping mass fields to ensure consistency and B.C.'s. */

  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    //    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	if (Grids[grid1]->GridData->CheckForOverlap(
				   SiblingList[grid1].GridList[grid2],
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyOverlappingMassField) == FAIL) {
	  fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
	  return FAIL;
	}

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	if (Grids[grid1]->GridData->CheckForOverlap(
				   SiblingList[grid1].GridList[grid2],
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyOverlappingMassField) == FAIL) {
	  fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
	  return FAIL;
	}

    if (CommunicationReceiveHandler() == FAIL)
      return FAIL;

  } // end loop over batches of grids

  PerformanceTimers[21] += ReturnWallTime() - time1;
  time1 = ReturnWallTime();

  /* Compute the potential for the top grid. */

  if (level == 0)
    if (ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids)
	== FAIL) {
      fprintf(stderr, "Error in ComputePotentialFieldLevelZero.\n");
      return FAIL;
    }

  /* Compute a first iteration of the potential and share BV's. */

#define ITERATE_POTENTIAL
#ifdef ITERATE_POTENTIAL
      if (level > 0) {
	CopyPotentialFieldAverage = 1;
	for (int iterate = 0; iterate < 1; iterate++) {

	  if (iterate > 0)
	    CopyPotentialFieldAverage = 2;

	  int Dummy;
	  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
	    if (Grids[grid1]->GridData->SolveForPotential(Dummy, level,
							 EvaluateTime) 
		== FAIL) {
	      fprintf(stderr, "Error in grid->SolveForPotential.\n");
	      return FAIL;
	    }

	  for (StartGrid = 0; StartGrid < NumberOfGrids; 
	       StartGrid += GRIDS_PER_LOOP) {
	    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

	    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
	    //	    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
	    CommunicationReceiveIndex = 0;
	    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
	    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
		if (Grids[grid1]->GridData->CheckForOverlap(
				   SiblingList[grid1].GridList[grid2],
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyPotentialField) == FAIL) {
		  fprintf(stderr, "Error in grid->CopyPotentialField.\n");
		  return FAIL;
		}

	    CommunicationDirection = COMMUNICATION_SEND;
	    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
		if (Grids[grid1]->GridData->CheckForOverlap(
				   SiblingList[grid1].GridList[grid2],
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyPotentialField) == FAIL) {
		  fprintf(stderr, "Error in grid->CopyPotentialField.\n");
		  return FAIL;
		}

	    if (CommunicationReceiveHandler() == FAIL)
	      return FAIL;

	  } // end loop over batches of grids

	} /* loop over iterations */
	CopyPotentialFieldAverage = 0;
      } /* end: if (level > 0) */

#endif /* ITERATE_POTENTIAL */

  PerformanceTimers[22] += ReturnWallTime() - time1;
  time1 = ReturnWallTime();

  /* if level > MaximumGravityRefinementLevel, then do final potential
     solve (and acceleration interpolation) here rather than in the main
     EvolveLevel since it involves communications. */

  if (reallevel > MaximumGravityRefinementLevel) {

    /* compute potential and acceleration on coarser level [LOCAL]
       (but only if there is at least a subgrid -- it should be only
        if there is a subgrrid on reallevel, but this is ok). */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) 
      if (Grids[grid1]->NextGridNextLevel != NULL) {
	Grids[grid1]->GridData->SolveForPotential(level,
					       MaximumGravityRefinementLevel);
	Grids[grid1]->GridData->ComputeAccelerationField(
                           (HydroMethod == Zeus_Hydro) ? ZEUS_GRIDS : GRIDS,
					       MaximumGravityRefinementLevel);
      }

    /* Interpolate potential for reallevel grids from coarser grids. */

    int Dummy, GridCount;
    LevelHierarchyEntry *Temp, *LastTemp, *FirstTemp = LevelArray[reallevel];

    do {

      GridCount = 0;
      CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      CommunicationReceiveIndex = 0;
      CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      Temp = FirstTemp;
      while (Temp != NULL && GridCount++ < GRIDS_PER_LOOP) {
	HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
	for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	  Temp3 = Temp3->ParentGrid;
	if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) 
	    == FAIL) {
	  fprintf(stderr, "Error in grid->InterpolateAccelerations.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
      LastTemp = Temp;

      CommunicationDirection = COMMUNICATION_SEND;
      Temp = FirstTemp;
      while (Temp != LastTemp) {
	HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
	for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	  Temp3 = Temp3->ParentGrid;
	if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) 
                                                                   == FAIL) {
	  fprintf(stderr, "Error in grid->InterpolateAccelerations.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
      FirstTemp = LastTemp;

      if (CommunicationReceiveHandler() == FAIL)
	return FAIL;

    } while (LastTemp != NULL);

  } // end: if (reallevel > MaximumGravityRefinementLevel)

  PerformanceTimers[5] += ReturnWallTime() - start_time;
  return SUCCESS;
}


/* ======================================================================= */
/* This routines does the flux correction and project for all grids on this
   level from the list of subgrids. */

int UpdateFromFinerGrids(HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[], 
			 fluxes **SubgridFluxesEstimate[])

{

  int grid1, subgrid, StartGrid, EndGrid;
  HierarchyEntry *NextGrid;
  double start_time = ReturnWallTime();

  /* Define a temporary flux holder for the refined fluxes. */

  fluxes SubgridFluxesRefined;

  /* For each grid,
     (a) project the subgrid's solution into this grid (step #18) */

  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* -------------- FIRST PASS ----------------- */

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* Loop over subgrids for this grid. */

      NextGrid = Grids[grid1]->NextGridNextLevel;
      subgrid = 0;
      CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      while (NextGrid != NULL && FluxCorrection) {

	/* Project subgrid's refined fluxes to the level of this grid. */

#ifdef USE_MPI
	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = grid1;
	CommunicationReceiveArgumentInt[1][CommunicationReceiveIndex] = subgrid;
#endif /* USE_MPI */
	    
	if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		       Grids[grid1]->GridData, SubgridFluxesRefined) == FAIL) {
	  fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	  return FAIL;
	}

	NextGrid = NextGrid->NextGridThisLevel;
	subgrid++;
      }

    } // end of loop over subgrids

    /* -------------- SECOND PASS ----------------- */

    CommunicationDirection = COMMUNICATION_SEND;

    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* Loop over subgrids for this grid. */

      NextGrid = Grids[grid1]->NextGridNextLevel;
      subgrid = 0;
      while (NextGrid != NULL && FluxCorrection) {

	/* Project subgrid's refined fluxes to the level of this grid. */

	if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		      Grids[grid1]->GridData, SubgridFluxesRefined) == FAIL) {
	  fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	  return FAIL;
	}

        /* Correct this grid for the refined fluxes (step #19)
           (this also deletes the fields in SubgridFluxesRefined). 
	   (only call it if the grid and sub-grid are on the same
	    processor, otherwise handled in CommunicationReceiveHandler.) */

	if (NextGrid->GridData->ReturnProcessorNumber() ==
	    Grids[grid1]->GridData->ReturnProcessorNumber())
	  if (Grids[grid1]->GridData->CorrectForRefinedFluxes
	      (SubgridFluxesEstimate[grid1][subgrid], &SubgridFluxesRefined,
	       SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1]     )
	      == FAIL) {
	    fprintf(stderr, "Error in grid->CorrectForRefinedFluxes.\n");
	    return FAIL;
	  }

	NextGrid = NextGrid->NextGridThisLevel;
	subgrid++;
      }

    } // end of loop over subgrids

    /* -------------- THIRD PASS ----------------- */

    if (CommunicationReceiveHandler(SubgridFluxesEstimate,
				    NumberOfSubgrids) == FAIL)
      return FAIL;

  } // end of loop over batches of grids

  /* For each grid,
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */

  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* -------------- FIRST PASS ----------------- */

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* Loop over subgrids for this grid: replace solution. */

      CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      NextGrid = Grids[grid1]->NextGridNextLevel;
      while (NextGrid != NULL) {

	/* Project the subgrid solution into this grid. */

	if (NextGrid->GridData->ProjectSolutionToParentGrid
	                                   (*Grids[grid1]->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	  return FAIL;
	}

	NextGrid = NextGrid->NextGridThisLevel;
      }

    } // end of loop over subgrids

    /* -------------- SECOND PASS ----------------- */

    CommunicationDirection = COMMUNICATION_SEND;

    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* Loop over subgrids for this grid: replace solution. */

      NextGrid = Grids[grid1]->NextGridNextLevel;
      while (NextGrid != NULL) {

	/* Project the subgrid solution into this grid. */

	if (NextGrid->GridData->ProjectSolutionToParentGrid
	                                   (*Grids[grid1]->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	  return FAIL;
	}

	NextGrid = NextGrid->NextGridThisLevel;
      }

    } // end of loop over subgrids

    /* -------------- THIRD PASS ----------------- */

    if (CommunicationReceiveHandler() == FAIL)
      return FAIL;

  } // end of loop over batches of grids

  PerformanceTimers[6] += ReturnWallTime() - start_time;

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

  typedef HierarchyEntry* HierarchyEntryPointer;
  *Grids = new HierarchyEntryPointer[NumberOfGrids];
  Temp = LevelArray[level];
  while (Temp != NULL) {
    (*Grids)[counter++] = Temp->GridHierarchyEntry;
    Temp              = Temp->NextGridThisLevel;
  }

  return NumberOfGrids;
}

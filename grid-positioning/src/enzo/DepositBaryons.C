/***********************************************************************
/
/  DEPOSIT BARYONS IN THIS GRID AND ALL SUBGRIDS INTO
/      GRAVITATINGMASSFIELD IN THIS GRID
/
/  written by: Greg Bryan
/  date:       March, 1999
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */ 

#include <stdio.h>
#include "ErrorExceptions.h"
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
 
int DepositBaryonsChildren(HierarchyEntry *DepositGrid,
			   HierarchyEntry *Grid, FLOAT Time);
 
int DepositBaryons(HierarchyEntry *Grid, FLOAT When)
{
 
  /* Get the time and dt for this grid.  Compute time+1/2 dt. */
 
  FLOAT TimeMidStep =     Grid->GridData->ReturnTime() +
                      When*Grid->GridData->ReturnTimeStep();
 
  if (CommunicationDirection == COMMUNICATION_SEND ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
 
    /* Set the under_subgrid field (indicating if a cell is refined or not)
       on this grid. */
 
    HierarchyEntry *Temp = Grid->NextGridNextLevel;
    Grid->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
    while (Temp != NULL) {
      Grid->GridData->ZeroSolutionUnderSubgrid(Temp->GridData,
					       ZERO_UNDER_SUBGRID_FIELD);
      Temp = Temp->NextGridThisLevel;
    }
 
  } // end: (CommunicationDirection != SEND || SEND_RECEIVE)
 
  /* Deposit baryons to GravitatingMassField in this grid. */
 
  if (Grid->GridData->DepositBaryons(Grid->GridData, TimeMidStep) == FAIL) {
    fprintf(stderr, "Error in grid->DepositBaryons.\n");
    ENZO_FAIL("");
  }
 
  /* Recursively deposit baryons in children (at TimeMidStep). */
 
  if (Grid->NextGridNextLevel != NULL)
    if (DepositBaryonsChildren(Grid, Grid->NextGridNextLevel, TimeMidStep)
	== FAIL) {
      fprintf(stderr, "Error in DepositBaryonsChildren.\n");
      ENZO_FAIL("");
    }
 
  return SUCCESS;
}
 
 
/* this doesn't quite work yet (subgrid zeroing needs to be worked out). */
 
int DepositBaryonsChildren(HierarchyEntry *DepositGrid,
			   HierarchyEntry *Grid, FLOAT DepositTime)
{
 
  /* Set the field indicating if a cell is refined or not. */
 
  if (CommunicationDirection == COMMUNICATION_SEND ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
    HierarchyEntry *Temp = Grid->NextGridNextLevel;
    Grid->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
    while (Temp != NULL) {
      Grid->GridData->ZeroSolutionUnderSubgrid(Temp->GridData,
					       ZERO_UNDER_SUBGRID_FIELD);
      Temp = Temp->NextGridThisLevel;
    }
  }
 
  /* Deposit baryons in Grid into DepositGrid at the given time. */
 
  if (Grid->GridData->DepositBaryons(DepositGrid->GridData, DepositTime)
      == FAIL) {
    fprintf(stderr, "Error in grid->DepositBaryons.\n");
    ENZO_FAIL("");
  }
 
  /* Next grid on this level. */
 
  if (Grid->NextGridThisLevel != NULL)
    if (DepositBaryonsChildren(DepositGrid, Grid->NextGridThisLevel,
			       DepositTime) == FAIL) {
      fprintf(stderr, "Error in DepositBaryonsChildren(1).\n");
      ENZO_FAIL("");
    }
 
  /* Recursively deposit baryons in children. */
 
  if (Grid->NextGridNextLevel != NULL)
    if (DepositBaryonsChildren(DepositGrid, Grid->NextGridNextLevel,
			       DepositTime) == FAIL) {
      fprintf(stderr, "Error in DepositBaryonsChildren(2).\n");
      ENZO_FAIL("");
    }
 
  return SUCCESS;
}

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
/  DEPOSIT PARTICLES INTO PARTICLEMASSFIELD IN THIS GRID AND ALL SUBGRIDS
/
/  written by: Greg Bryan
/  date:       May, 1995
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

/* function prototypes */

int DepositParticleMassFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT Time);

int DepositParticleMassField(HierarchyEntry *Grid, FLOAT TimeMidStep)
{

  /* Get the time and dt for this grid.  Compute time+1/2 dt. */

  if (TimeMidStep < 0)
    TimeMidStep =     Grid->GridData->ReturnTime() + 
                  0.5*Grid->GridData->ReturnTimeStep();

  if (CommunicationDirection != COMMUNICATION_RECEIVE) {

    /* Initialize the gravitating mass field parameters (if necessary). */

    if (Grid->GridData->InitializeGravitatingMassFieldParticles(RefineBy) 
                                                                  == FAIL) {
      fprintf(stderr, 
	      "Error in grid->InitializeGravitatingMassFieldParticles.\n");
      return FAIL;
    }

    /* Clear the GravitatingMassFieldParticles. */

    if (Grid->GridData->ClearGravitatingMassFieldParticles() == FAIL) {
      fprintf(stderr, "Error in grid->ClearGravitatingMassFieldParticles.\n");
      return FAIL;
    }

  } // end: if (CommunicationDirection != COMMUNICATION_SEND)

  /* Deposit particles to GravitatingMassFieldParticles in this grid. */

  if (Grid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
				 GRAVITATING_MASS_FIELD_PARTICLES) == FAIL) {
    fprintf(stderr, "Error in grid->DepositParticlePositions.\n");
    return FAIL;
  }

  /* Recursively deposit particles in children (at TimeMidStep). */

  if (Grid->NextGridNextLevel != NULL)
    if (DepositParticleMassFieldChildren(Grid, Grid->NextGridNextLevel, 
					 TimeMidStep) 
	== FAIL) {
      fprintf(stderr, "Error in DepositParticleMassFieldChildren.\n");
      return FAIL;
    }

  return SUCCESS;
}




int DepositParticleMassFieldChildren(HierarchyEntry *DepositGrid, 
				     HierarchyEntry *Grid, FLOAT DepositTime)
{
 
  /* Deposit particles in Grid into DepositGrid at the given time. */

  if (Grid->GridData->DepositParticlePositions(DepositGrid->GridData, 
		     DepositTime, GRAVITATING_MASS_FIELD_PARTICLES) == FAIL) {
    fprintf(stderr, "Error in grid->DepositParticlePositions.\n");
    return FAIL;
  }

  /* Next grid on this level. */

  if (Grid->NextGridThisLevel != NULL)
    if (DepositParticleMassFieldChildren(DepositGrid, Grid->NextGridThisLevel, 
					 DepositTime) == FAIL) {
      fprintf(stderr, "Error in DepositParticleMassFieldChildren(1).\n");
      return FAIL;
    }

  /* Recursively deposit particles in children. */

  if (Grid->NextGridNextLevel != NULL)
    if (DepositParticleMassFieldChildren(DepositGrid, Grid->NextGridNextLevel, 
					 DepositTime) == FAIL) {
      fprintf(stderr, "Error in DepositParticleMassFieldChildren(2).\n");
      return FAIL;
    }


  return SUCCESS;
}

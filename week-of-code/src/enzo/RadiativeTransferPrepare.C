/***********************************************************************
/
/  PREPARE THE RADIATIVE TRANSFER MODULE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/ PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "Star.h"

int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove,
				     int level);
int StarParticleRadTransfer(LevelHierarchyEntry *LevelArray[], int level,
			    Star *AllStars);
int RestartPhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		   Star *AllStars);

int RadiativeTransferPrepare(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *&AllStars,
			     float dtLevelAbove)
{

  /* Return if this does not concern us */
  if (!(RadiativeTransfer)) return SUCCESS;

  FLOAT GridTime, dt;
  GridTime = LevelArray[level]->GridData->ReturnTime();
  dt = LevelArray[level]->GridData->ReturnTimeStep();

  //if (dtPhoton >= 0.0 && GridTime+dt >= PhotonTime) {
  if (GridTime+dt >= PhotonTime) {

    /* Determine the photon timestep */

    RadiativeTransferComputeTimestep(LevelArray, MetaData, dtLevelAbove, level);

    /* Convert star particles into radiation sources only if we're going
       into EvolvePhotons */

    StarParticleRadTransfer(LevelArray, level, AllStars);

  } // ENDIF
  
  /* If the first timestep after restart and we have radiation
     sources, go back a light crossing time of the box and run
     EvolvePhotons to populate the grids with the proper rates. */

  if (MetaData->FirstTimestepAfterRestart == TRUE)
    if (RestartPhotons(MetaData, LevelArray, AllStars) == FAIL) {
      fprintf(stderr, "Error in RestartPhotons.\n");
      ENZO_FAIL("");
    }

  return SUCCESS;

}

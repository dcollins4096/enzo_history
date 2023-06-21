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
/  CHECK FOR TIME ACTION
/
/  written by: Greg Bryan
/  date:       September, 2000
/  modified1:
/
/  PURPOSE:
/    This routine checks to see if the time has arrived for a given
/      "time-action" to occur (i.e. some action that is applied to
/       the data at a given time/redshift).
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
#include "CosmologyParameters.h"

		 
int CheckForTimeAction(LevelHierarchyEntry *LevelArray[], 
		       TopGridData &MetaData) 
{

  /* Declarations. */

  int i, level;

  /* Check for time against list of actions. */

  for (i = 0; i < MAX_TIME_ACTIONS; i++)
    if (MetaData.Time >= TimeActionTime[i] && TimeActionTime[i] > 0) {

      if (debug)
	printf("Applying TimeAction %d at t=%"GOUTSYM"\n", i, TimeActionTime[i]);

      /* Done, turn it off (-1 in redshift indicates off). */

      TimeActionTime[i] = 0;
      TimeActionRedshift[i] = -1;

      /* Now apply to all grids. */

      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	LevelHierarchyEntry *Temp = LevelArray[level];
	while (Temp != NULL) {
	  if (Temp->GridData->ApplyTimeAction(TimeActionType[i],
					   TimeActionParameter[i]) == FAIL) {
	    fprintf(stderr, "Errot in grid->ApplyTimeAction\n");
	    return FAIL;
	  }
	  Temp = Temp->NextGridThisLevel;
	}
      }

    }

  return SUCCESS;
}

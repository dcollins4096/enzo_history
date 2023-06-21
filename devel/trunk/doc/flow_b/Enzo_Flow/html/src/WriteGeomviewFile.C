/***********************************************************************
/
/  WRITE OUT A GEOMVIEW FILE OF MESH HIERARCHY
/
/  written by: James Bordner
/  date:       2003-09-13
/  modified1:
/
/  PURPOSE: This function walks through the grids and writes out 
/     a box for each one.
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
#include "message.h"

int WriteGeomviewHierarchy(char *basename,
			   int filenumber,
			   LevelHierarchyEntry *LevelArray[])
{

  int level;

  // Determine file name

  char filename[80];

  sprintf (filename, "%s.p%d.n%d",basename,MyProcessorNumber,filenumber);
  DEBUG_MESSAGE;
  printf ("filename=%s\n",filename);

  // Loop over all levels

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    // Loop over all grids in level

    LevelHierarchyEntry *LHE = LevelArray[level];

    while (LHE != NULL) {

      printf ("level %d\n",level);

      LHE = LHE->NextGridThisLevel;

    } // end loop over grids

  } // end loop over levels

  return SUCCESS;
}

int WriteGeomviewLevel(char *basename,
		       int filenumber,
		       int level,
		       LevelHierarchyEntry *LevelArray[])
{
  return SUCCESS;
}

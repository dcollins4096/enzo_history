/***********************************************************************
/
/  WRITE OUT THE DATA HIERARCHY (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// This function writes out the data hierarchy (TopGrid)
 
#include <string.h>
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
 
/* function prototypes */
 
 
 
int WriteDataHierarchy(FILE *fptr, HierarchyEntry *Grid,
		       char* base_name, int &GridID, FLOAT WriteTime)
{
 
  int OriginalID, NextGridThisLevelID, NextGridNextLevelID;
 
  /* Write out header info for this grid */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "\nGrid = %"ISYM"\n", GridID);
  OriginalID = GridID;
 
  /* Write out grid data for this grid (if WriteTime is < 0 then output
     at the grid's time, otherwise interpolate to WriteTime and output). */
 
  if (WriteTime < 0) {
    if (Grid->GridData->WriteGrid(fptr, base_name, GridID) == FAIL) {
      fprintf(stderr, "Error in grid->WriteGrid.\n");
      return FAIL;
    }
  }
  else
    if (Grid->GridData->WriteGridInterpolate(WriteTime, fptr, base_name,
					     GridID) == FAIL) {
      fprintf(stderr, "Error in grid->WriteGrid.\n");
      return FAIL;
    }
 
  /* Write out pointer information for the next grid this level */
 
  NextGridThisLevelID = GridID + 1;
  if (Grid->NextGridThisLevel == NULL) NextGridThisLevelID = 0;
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "Pointer: Grid[%"ISYM"]->NextGridThisLevel = %"ISYM"\n", OriginalID,
	    NextGridThisLevelID);
 
  if (NextGridThisLevelID != 0) {
    GridID++;
    if (WriteDataHierarchy(fptr, Grid->NextGridThisLevel, base_name, GridID,
			   WriteTime) == FAIL) {
      fprintf(stderr, "Error in WriteDataHierarchy(1).\n");
      return FAIL;
    }
  }
 
  /* Write out pointer information for the next grid next level */
 
  NextGridNextLevelID = GridID + 1;
  if (Grid->NextGridNextLevel == NULL) NextGridNextLevelID = 0;
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "Pointer: Grid[%"ISYM"]->NextGridNextLevel = %"ISYM"\n", OriginalID,
	    NextGridNextLevelID);
 
  if (NextGridNextLevelID != 0) {
    GridID++;
    if (WriteDataHierarchy(fptr, Grid->NextGridNextLevel, base_name, GridID,
			   WriteTime) == FAIL) {
      fprintf(stderr, "Error in WriteDataHierarchy(1).\n");
      return FAIL;
    }
  }
 
  return SUCCESS;
}

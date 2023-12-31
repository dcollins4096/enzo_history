/***********************************************************************
/
/  WRITE OUT THE DATA HIERARCHY (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2006
/              Pass through HDF5 file id to group routines
/
/  PURPOSE:
/
************************************************************************/
 
// This function writes out the data hierarchy (TopGrid)

#include <hdf5.h> 
#include <string.h>
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
 
/* function prototypes */
 
 
 
int Group_WriteDataHierarchy(FILE *fptr, TopGridData &MetaData, HierarchyEntry *Grid,
		            char* base_name, int &GridID, FLOAT WriteTime, hid_t file_id,
                    int CheckpointDump = FALSE)
{
 
  int OriginalID, NextGridThisLevelID, NextGridNextLevelID;
 
  /* Write out header info for this grid */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "\nGrid = %"ISYM"\n", GridID);
  OriginalID = GridID;
 
  /* Write out grid data for this grid (if WriteTime is < 0 then output
     at the grid's time, otherwise interpolate to WriteTime and output). */
 
  if ((WriteTime < 0) || (CheckpointDump == TRUE)) {
    if (Grid->GridData->Group_WriteGrid(fptr, base_name, GridID, file_id,
            CheckpointDump) == FAIL) {
      fprintf(stderr, "Error in grid->Group_WriteGrid.\n");
      ENZO_FAIL("");
    }
  }
  else
    if (Grid->GridData->Group_WriteGridInterpolate(WriteTime, fptr, base_name, GridID, file_id) == FAIL) {
      fprintf(stderr, "Error in grid->Group_WriteGridInterpolate.\n");
      ENZO_FAIL("");
    }
 
  /* Write out pointer information for the next grid this level */
 
  NextGridThisLevelID = GridID + 1;
  if (Grid->NextGridThisLevel == NULL) NextGridThisLevelID = 0;
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "Pointer: Grid[%"ISYM"]->NextGridThisLevel = %"ISYM"\n", OriginalID,
	    NextGridThisLevelID);
 
  if (NextGridThisLevelID != 0) {
    GridID++;
    if (Group_WriteDataHierarchy(fptr, MetaData, Grid->NextGridThisLevel,
            base_name, GridID, WriteTime, file_id, CheckpointDump) == FAIL) {
      fprintf(stderr, "Error in Group_WriteDataHierarchy(1).\n");
      ENZO_FAIL("");
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
    if (Group_WriteDataHierarchy(fptr, MetaData, Grid->NextGridNextLevel,
                base_name, GridID, WriteTime, file_id, CheckpointDump) == FAIL) {
      fprintf(stderr, "Error in Group_WriteDataHierarchy(1).\n");
      ENZO_FAIL("");
    }
  }
 
  return SUCCESS;
}

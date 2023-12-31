/***********************************************************************
/
/  WRITE OUT UNIGRID CUBES AT ALL TIMES
/
/  written by: Robert Harkness
/  date:       October, 2004
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
 
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
#include "CosmologyParameters.h"
 
// function prototypes
 
int FileMover(char *file_to_move);
int MakeDir(char *subdirectory);
int SysMkdir(char *startdir, char *directory);
 
int WriteDataCubes(HierarchyEntry *TopGrid, int TDdims[], char *gridbasename, int &GridID, FLOAT WriteTime);
 
int CommunicationCombineGrids(HierarchyEntry *OldHierarchy,
			      HierarchyEntry **NewHierarchyPointer,
			      FLOAT WriteTime);
void DeleteGridHierarchy(HierarchyEntry *GridEntry);
 
 
 
 
int WriteAllDataCubes(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, FLOAT WriteTime = -1)
{
 
  char cycle_id[9];
  char name[MAX_LINE_LENGTH];
 
  int GridID = 1;
  int GridJD = 1;
 
  /* If this is an interpolated time step, then temporary replace the time
     in MetaData. */
 
  FLOAT SavedTime = (WriteTime < 0) ? MetaData.Time : WriteTime;
 
  //  cycle number
  sprintf(cycle_id, "%8.8"ISYM, MetaData.CycleNumber);
 
  strcpy(name, "TS");
  strcat(name, cycle_id);
  //  fprintf(stderr, "Dump cycle: %"ISYM"  %s\n", MetaData.CycleNumber, name);
 
  //  Start I/O timing
 
#ifdef USE_MPI
//  io_start = MPI_Wtime();
#endif /* USE_MPI */
 
  /* Combine the top level grids into a single grid for output
     (TempTopGrid is the top of an entirely new hierarchy). */
 
  HierarchyEntry *TempTopGrid;
  CommunicationCombineGrids(TopGrid, &TempTopGrid, WriteTime);
 
  int TGdims[3];
 
  TGdims[0] = MetaData.TopGridDims[0];
  TGdims[1] = MetaData.TopGridDims[1];
  TGdims[2] = MetaData.TopGridDims[2];
 
  if (CubeDumpEnabled == 1) {
    if (WriteDataCubes(TempTopGrid, TGdims, name, GridJD, WriteTime) == FAIL) {
      fprintf(stderr, "Error in WriteDataCubes\n");
      return FAIL;
    }
  }
 
  // Clean up combined top level grid, and first two levels of hierarchy
 
  if (TempTopGrid != TopGrid) {
    if (TempTopGrid->NextGridNextLevel != NULL)
      DeleteGridHierarchy(TempTopGrid->NextGridNextLevel);
    delete TempTopGrid->GridData;
    delete TempTopGrid;
  }
 
  // Replace the time in metadata with the saved value (above)
 
  MetaData.Time = SavedTime;
 
  MPI_Barrier(MPI_COMM_WORLD);
 
//  Stop I/O timing
 
#ifdef USE_MPI
//  io_stop = MPI_Wtime();
#endif /* USE_MPI */
 
  return SUCCESS;
}
 
/*
void DeleteGridHierarchy(HierarchyEntry *GridEntry)
{
  if (GridEntry->NextGridThisLevel != NULL)
     DeleteGridHierarchy(GridEntry->NextGridThisLevel);
 
  delete GridEntry;
 
  return;
}
*/

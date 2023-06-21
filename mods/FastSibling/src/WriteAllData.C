/***********************************************************************
/
/  WRITE OUT ALL THE DATA (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// This function writes out the data hierarchy (TopGrid), the External
//   Boundary (Exterior), the TopGridData, and the global_data.

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
#include "CosmologyParameters.h"

/* function prototypes */

int WriteDataHierarchy(FILE *fptr, HierarchyEntry *TopGrid, 
		       char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteParameterFile(FILE *fptr, TopGridData &MetaData);
int WriteStarParticleData(FILE *fptr);
int WriteRadiationData(FILE *fptr);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
void DeleteGridHierarchy(HierarchyEntry *GridEntry);
int CommunicationCombineGrids(HierarchyEntry *OldHierarchy,
			      HierarchyEntry **NewHierarchyPointer,
			      FLOAT WriteTime);
void DeleteGridHierarchy(HierarchyEntry *GridEntry);
double ReturnWallTime();

char BCSuffix[]        = ".boundary";
char GridSuffix[]      = ".grid";
char HierarchySuffix[] = ".hierarchy";
char hdfsuffix[]       = ".hdf";
char RadiationSuffix[] = ".radiation";

int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData, 
		 ExternalBoundary *Exterior, FLOAT WriteTime = -1)
{

  /* declarations */

  char id[5], *cptr, name[MAX_LINE_LENGTH];
  char gridbasename[MAX_LINE_LENGTH], hierarchyname[MAX_LINE_LENGTH],
       radiationname[MAX_LINE_LENGTH];
  FILE *fptr;
  int GridID = 1;
  double time1 = ReturnWallTime();

  /* If this is an interpolated time step, then temporary replace the time
     in MetaData. */

  FLOAT SavedTime = (WriteTime < 0) ? MetaData.Time : WriteTime;

  /* Create main name */

  strcpy(name, basename);
  if (ComovingCoordinates && (cptr = strstr(name, "RRRR"))) {
    FLOAT a, dadt;
    CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
    sprintf(cptr, "%4.4d", nint(100*((1 + InitialRedshift)/a - 1)));
  } else {
    sprintf(id, "%4.4d", filenumber);
    if (filenumber >= 0)
      strcat(name, id);
  }

  /* Debugging info. */

  if (debug)
    printf("WriteAllData: writing file %s.\n", name);

  /* Set MetaData.BoundaryConditionName. */

  if (MetaData.BoundaryConditionName != NULL)
    delete MetaData.BoundaryConditionName;
  MetaData.BoundaryConditionName = new char[MAX_LINE_LENGTH];
  strcpy(MetaData.BoundaryConditionName, name);
  strcat(MetaData.BoundaryConditionName, BCSuffix);

  /* Output TopGrid data. */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(name, "w")) == NULL) {
      fprintf(stderr, "Error opening output file %s.\n", name);
      return FAIL;
    }
    if (WriteTime >= 0)
      fprintf(fptr, "# WARNING! Interpolated output: level = %d\n",
	      MetaData.OutputFirstTimeAtLevel-1);
    if (WriteParameterFile(fptr, MetaData) == FAIL) {
      fprintf(stderr, "Error in WriteParameterFile.\n");
      return FAIL;
    }
    fclose(fptr);
  }

  /* Output Boundary condition info. */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(MetaData.BoundaryConditionName, "w")) == NULL) {
      fprintf(stderr, "Error opening boundary condition file: %s\n", 
	      MetaData.BoundaryConditionName);
      return FAIL;
    }
    strcat(MetaData.BoundaryConditionName, hdfsuffix);
    if (Exterior->WriteExternalBoundary(fptr, MetaData.BoundaryConditionName) 
	== FAIL) {
      fprintf(stderr, "Error in WriteExternalBoundary.\n");
      return FAIL;
    }
    fclose(fptr);
  }

  /* Create hierarchy name and grid base name. */

  strcpy(hierarchyname, name);
  strcat(hierarchyname, HierarchySuffix);

  strcpy(gridbasename, name);
  strcat(gridbasename, GridSuffix);

  /* Combine the top level grids into a single grid for output
     (TempTopGrid is the top of an entirely new hierarchy). */

  HierarchyEntry *TempTopGrid;
  CommunicationCombineGrids(TopGrid, &TempTopGrid, WriteTime);

  /* Output Data Hierarchy. */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((fptr = fopen(hierarchyname, "w")) == NULL) {
      fprintf(stderr, "Error opening hierarchy file %s.\n", hierarchyname);
      return FAIL;
    }

  if (WriteDataHierarchy(fptr, TempTopGrid, gridbasename, GridID, WriteTime) 
      == FAIL) {
    fprintf(stderr, "Error in WriteDataHierarchy.\n");
    return FAIL;
  }

  /* Clean up combined top level grid, and first two levels of hierarchy. */

  if (TempTopGrid != TopGrid) {
    if (TempTopGrid->NextGridNextLevel != NULL)
      DeleteGridHierarchy(TempTopGrid->NextGridNextLevel);
    delete TempTopGrid->GridData;
    delete TempTopGrid;
  }

  /* Output StarParticle data (actually just number of stars). */

  if (WriteStarParticleData(fptr) == FAIL) {
    fprintf(stderr, "Error in WriteStarParticleData.\n");
    return FAIL;
  }

  /* Create radiation name and write radiation data. */

  if (RadiationFieldType >= 10 && RadiationFieldType <= 11 && 
      MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *Radfptr;
    strcpy(radiationname, name);
    strcat(radiationname, RadiationSuffix);
    if ((Radfptr = fopen(radiationname, "w")) == NULL) {
      fprintf(stderr, "Error opening radiation file %s.\n", radiationname);
      return FAIL;
    }
    if (WriteRadiationData(Radfptr) == FAIL) {
      fprintf(stderr, "Error in WriteRadiationData.\n");
      return FAIL;
    }
    fclose(Radfptr);
  }

  if (MyProcessorNumber == ROOT_PROCESSOR)
    fclose(fptr);

  /* Replace the time in metadata with the saved value (above). */

  MetaData.Time = SavedTime;

  /* Debugging info. */

  if (debug)
    printf("WriteAllData: finished writing data.\n");
  PerformanceTimers[7] += ReturnWallTime() - time1;

  return SUCCESS;
}

void DeleteGridHierarchy(HierarchyEntry *GridEntry)
{
  if (GridEntry->NextGridThisLevel != NULL)
     DeleteGridHierarchy(GridEntry->NextGridThisLevel);

  delete GridEntry;

  return;
}

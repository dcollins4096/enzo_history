/***********************************************************************
/
/  READ ALL THE DATA (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness
/              January, 2006
/
/  PURPOSE:
/
************************************************************************/
 
// This function reads in the data hierarchy (TopGrid), the External
//   Boundary (Exterior), the TopGridData, and the global_data.
 
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>

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
 
int ReadDataHierarchy(FILE *fptr, HierarchyEntry *TopGrid, int GridID,
		      HierarchyEntry *ParentGrid);
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int ReadStarParticleData(FILE *fptr);
int ReadRadiationData(FILE *fptr);
int AssignGridToTaskMap(Eint64 *GridIndex, Eint64 *Mem, int Ngrids);
 
extern char RadiationSuffix[];
extern char HierarchySuffix[];
extern char hdfsuffix[];
extern char TaskMapSuffix[];
 
int ReadAllData(char *name, HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior)
 
{
 
  /* declarations */
 
  char hierarchyname[MAX_LINE_LENGTH], radiationname[MAX_LINE_LENGTH];
  char taskmapname[MAX_LINE_LENGTH];

  FILE *fptr;
  FILE *tptr;

  int GridID = 1;
  int GridKD = 1;

  float dummy;
 
  /* Read TopGrid data. */
 
  if ((fptr = fopen(name, "r")) == NULL) {
    fprintf(stderr, "Error opening input file %s.\n", name);
    return FAIL;
  }
  if (ReadParameterFile(fptr, MetaData, &dummy) == FAIL) {
    fprintf(stderr, "Error in ReadParameterFile.\n");
    return FAIL;
  }
 
  /* Close main file. */
 
  fclose(fptr);
 
  /* Read Boundary condition info. */
 
  if ((fptr = fopen(MetaData.BoundaryConditionName, "r")) == NULL) {
    fprintf(stderr, "Error opening boundary condition file: %s\n",
	    MetaData.BoundaryConditionName);
    return FAIL;
  }
  if (Exterior->ReadExternalBoundary(fptr) == FAIL) {
    fprintf(stderr, "Error in ReadExternalBoundary (%s).\n",
	    MetaData.BoundaryConditionName);
    return FAIL;
  }
  strcat(MetaData.BoundaryConditionName, hdfsuffix);
  fclose(fptr);

  /* Create the task map name */
 
  strcpy(taskmapname, name);
  strcat(taskmapname, TaskMapSuffix);

  /* Read the task map */

  MPI_Barrier(MPI_COMM_WORLD);

  // fprintf(stderr, "All at sync point - read TaskMap and assign Grid->Task\n");

  if ((tptr = fopen(taskmapname, "r")) == NULL) {
    fprintf(stderr, "Error opening TaskMap file %s.\n", taskmapname);
    return FAIL;
  }

  Eint64 GridIndex[MAX_NUMBER_OF_TASKS], OldPN, Mem[MAX_NUMBER_OF_TASKS];
  int i, ntask;
  i = 0;
  while( fscanf(tptr, "Grid %4lld  PN %4lld  Memory %16lld\n", &GridIndex[i], &OldPN, &Mem[i]) != EOF) {
    // fprintf(stderr, "Grid %4lld  PN %4lld  Memory %16lld\n", GridIndex[i], OldPN, Mem[i]);
    i++;
  }
  ntask = i;
  // fprintf(stderr, "Read memory size for %"ISYM" tasks\n", ntask);

  if (AssignGridToTaskMap(GridIndex, Mem, ntask) == FAIL) {
    fprintf(stderr, "Error in AssignGridToTaskMap.\n");
    return FAIL;
  }

  fclose(tptr);

  /* Create hierarchy name and grid base name. */

  strcpy(hierarchyname, name);
  strcat(hierarchyname, HierarchySuffix);
 
  /* Read Data Hierarchy. */
 
  if ((fptr = fopen(hierarchyname, "r")) == NULL) {
    fprintf(stderr, "Error opening hierarchy file %s.\n", hierarchyname);
    return FAIL;
  }
  GridID = 1;
  if (ReadDataHierarchy(fptr, TopGrid, GridID, NULL) == FAIL) {
    fprintf(stderr, "Error in ReadDataHierarchy (%s).\n", hierarchyname);
    return FAIL;
  }

  /* Read StarParticle data. */
 
  if (ReadStarParticleData(fptr) == FAIL) {
    fprintf(stderr, "Error in ReadStarParticleData.\n");
    return FAIL;
  }
 
  /* Create radiation name and read radiation data. */
 
  if (RadiationFieldType >= 10 && RadiationFieldType <= 11) {
    FILE *Radfptr;
    strcpy(radiationname, name);
    strcat(radiationname, RadiationSuffix);
    if ((Radfptr = fopen(radiationname, "r")) == NULL) {
      fprintf(stderr, "Error opening radiation file %s.\n", name);
      return FAIL;
    }
    if (ReadRadiationData(Radfptr) == FAIL) {
      fprintf(stderr, "Error in ReadRadiationData.\n");
      return FAIL;
    }
    fclose(Radfptr);
  }
 
  fclose(fptr);
 
  return SUCCESS;
}

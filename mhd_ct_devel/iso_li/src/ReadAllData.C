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
/  READ ALL THE DATA (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// This function reads in the data hierarchy (TopGrid), the External
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

/* function prototypes */
#ifdef HAOXU 
  int ShockPoolRestart(FILE *fptr,  HierarchyEntry *TopGrid,
                       TopGridData &MetaData);
#endif

int ReadDataHierarchy(FILE *fptr, HierarchyEntry *TopGrid, int GridID,
		      HierarchyEntry *ParentGrid);
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int ReadStarParticleData(FILE *fptr);
int ReadRadiationData(FILE *fptr);

extern char RadiationSuffix[];
extern char HierarchySuffix[];
extern char hdfsuffix[];

int ReadAllData(char *name, HierarchyEntry *TopGrid, TopGridData &MetaData, 
		 ExternalBoundary *Exterior)

{

  /* declarations */

  char hierarchyname[MAX_LINE_LENGTH], radiationname[MAX_LINE_LENGTH];
  FILE *fptr;
  int GridID = 1;
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

#ifdef HAOXU
       /* For Shock Pool problem, read in the shock conditions */
  if ((fptr = fopen(name, "r")) == NULL) {
    fprintf(stderr, "Error opening input file %s.\n", name);
    return FAIL;
  }

  if (ProblemType == 3)
    if (ShockPoolRestart(fptr, TopGrid, MetaData) == FAIL) {
      fprintf(stderr, "Error in ShockPoolRestart.\n");
      return FAIL;
    }

#endif

  fclose(fptr);

  return SUCCESS;
}

/***********************************************************************
/
/  CHECKS THE ACCURACY OF THE GRAVITY SOLVER (FOR SPHERE)
/
/  written by: Greg Bryan
/  date:       September, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

/* function prototypes */


char TGSOutputFileName[] = "TestGravitySphereCheckResults.out";
  

int TestGravitySphereCheckResults(LevelHierarchyEntry *LevelArray[])
{

  /* declarations */

  int level;
  FILE *fptr;
  char name[MAX_LINE_LENGTH], proc[5];

  /* Open output file. */

  strcpy(name, TGSOutputFileName);
  if (NumberOfProcessors > 1) {
    sprintf(proc, "%4.4d", MyProcessorNumber);
    strcat(name, proc);
  }
  if ((fptr = fopen(name, "w")) == NULL) {
    fprintf(stderr, "TestGravitySphereCheckResults:\nError opening file %s.\n",
	    name);
    exit(FAIL);
  }

  /* Output header. */

  fprintf(fptr, "#     r         f_tang         f_radial        f_analytic\n");
  
  /* For each grid on each level, check the results. */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    LevelHierarchyEntry* Level = LevelArray[level];

    while (Level != NULL) {

      if (Level->GridData->TestGravitySphereCheckResults(fptr) == FAIL) {
	fprintf(stderr, "Error in grid->TestGravitySphereCheckResults\n");
	exit(FAIL);
      }

      /* Next grid on the this level. */

      Level = Level->NextGridThisLevel;

    }

  } // end: loop over levels

  /* Close output file. */

  fclose(fptr);

  return SUCCESS;
}

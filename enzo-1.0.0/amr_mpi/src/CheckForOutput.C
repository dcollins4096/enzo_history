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
/  CHECK FOR OUTPUT
/
/  written by: Greg Bryan
/  date:       January, 1996
/  modified1:
/
/  PURPOSE:
/    This routine checks a number of criteria for output and then calls
/      the appropriate routine.
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

/* function prototypes */

int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid, 
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);
		 
int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData, 
		   ExternalBoundary *Exterior, int &WroteData)
{

  /* Declarations. */

  char *Name;
  int i, Number;
  WroteData = FALSE;

  /* Check for output: time-based. */

  if (MetaData.Time >= MetaData.TimeLastDataDump + MetaData.dtDataDump
      && MetaData.dtDataDump > 0.0) {
    MetaData.TimeLastDataDump += MetaData.dtDataDump;
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
      }
    WroteData = TRUE;
  }

  /* Check for output: cycle-based. */

  if (MetaData.CycleNumber >= MetaData.CycleLastDataDump + 
                              MetaData.CycleSkipDataDump   &&
      MetaData.CycleSkipDataDump > 0) {
    MetaData.CycleLastDataDump += MetaData.CycleSkipDataDump;
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
      }
    WroteData = TRUE;
  }

  /* Check for output: redshift-based. */

  if (ComovingCoordinates)
    for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
      if (CosmologyOutputRedshift[i] != -1)
	if (MetaData.Time >= CosmologyOutputRedshiftTime[i]) {
	  CosmologyOutputRedshift[i] = -1; // done, turn it off
	  if (CosmologyOutputRedshiftName[i] == NULL) {
	    Name   = MetaData.RedshiftDumpName;
	    Number = i;   // append number to end of name
	  } 
	  else {
	    Name   = CosmologyOutputRedshiftName[i];
	    Number = -1;  // Don't append number (####) to end of name
	  }
	  if (WriteAllData(Name, Number, TopGrid, MetaData,Exterior) == FAIL) {
	    fprintf(stderr, "Error in WriteAllData.\n");
	    return FAIL;
	  }
	  WroteData = TRUE;
	}

  return SUCCESS;
}

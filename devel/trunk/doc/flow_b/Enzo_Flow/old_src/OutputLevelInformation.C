/***********************************************************************
/
/  COMPUTE AND OUTPUT SOME SUMMARY INFORMATION ABOUT HIERARCHY
/
/  written by: Greg Bryan
/  date:       September, 1996
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "hdf4.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"

#if defined(MALLOC_IRIS4)
#include <sys/types.h>
#include <malloc.h>
#endif /* MALLOC_IRIS4 */


int OutputLevelInformation(FILE *fptr, TopGridData &MetaData,
			   LevelHierarchyEntry *LevelArray[]) 
{

  /* Only for ROOT_PROCESSOR. */

  if (MyProcessorNumber != ROOT_PROCESSOR)
    return SUCCESS;

  /* Declarations */

  int level, maxdepth = 0, GridMemory, NumberOfCells;
  float AxialRatio, GridVolume;
  int NumberOfGrids[MAX_DEPTH_OF_HIERARCHY],
      Cells[MAX_DEPTH_OF_HIERARCHY], FlaggedCells[MAX_DEPTH_OF_HIERARCHY];
  float Coverage[MAX_DEPTH_OF_HIERARCHY], Memory[MAX_DEPTH_OF_HIERARCHY],
        MeanAxialRatio[MAX_DEPTH_OF_HIERARCHY],
        FractionFlagged[MAX_DEPTH_OF_HIERARCHY];
  LevelHierarchyEntry *Temp;

  /* Zero Hierarchy sums. */

  int   HierarchyNumberOfGrids = 0;
  float HierarchyAxialRatio    = 0, HierarchyMemory = 0;

  /* Loop over levels, counting grids & approx memory used. */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    /* Zero this level's summary info. */

    NumberOfGrids[level]  = 0;
    Memory[level]         = 0;
    Coverage[level]       = 0;
    MeanAxialRatio[level] = 0;
    Cells[level]          = 0;
    FlaggedCells[level]   = 0;

    /* Loop over grids on this level. */

    Temp = LevelArray[level];
    while (Temp != NULL) {
      
      /* Incorporate grid information. */

      maxdepth = level;

      Temp->GridData->CollectGridInformation(GridMemory, GridVolume, 
					     NumberOfCells, AxialRatio);
      NumberOfGrids[level]++;
      Memory[level]         += float(GridMemory);
      Coverage[level]       += GridVolume;
      MeanAxialRatio[level] += AxialRatio;

      /* Get the actual number of flagged cells. */

#ifdef UNUSED
      int dummy;
      if (MetaData->StaticHierarchy == FALSE) {
	Temp->GridData->ClearFlaggingField();
	if (Temp->GridData->SetFlaggingField(dummy, level-1) == FAIL) {
	  fprintf(stderr, "Error in grid->SetFlaggingField.\n");
	  return FAIL;
	}
	FlaggedCells[level] += Temp->GridData->FlagBufferZones();
	Cells[level]        += NumberOfCells;
	Temp->GridData->DeleteFlaggingField();
      }
#endif /* UNUSED */

      /* Next grid */

      Temp = Temp->NextGridThisLevel;
    }

    /* Add to overall sums. */

    HierarchyNumberOfGrids += NumberOfGrids[level];
    HierarchyMemory        += Memory[level];
    HierarchyAxialRatio    += MeanAxialRatio[level];

    if (NumberOfGrids[level] == 0) {
      MeanAxialRatio[level] = 0;
      FractionFlagged[level] = 0;
    } else {
      MeanAxialRatio[level] /= float(NumberOfGrids[level]);
      FractionFlagged[level] = float(FlaggedCells[level])/
	float(max(Cells[level], 1));
    }

  }

  /* Get the total memory, if possible. */

  float TotalMemoryDeclared = 0, TotalMemoryUsed = 0;
#if defined(MALLOC_IRIS4)
  struct mallinfo proc;
  proc = mallinfo();
  TotalMemoryDeclared = float(proc.arena)/1.049e6;
  TotalMemoryUsed     = float(proc.usmblks+proc.uordblks)/1.049e6;
#endif /* MALLOC_IRIS4 */

  /* Write output (memory in MB). */

  fprintf(fptr, "%"GOUTSYM"  %d %d %g %g %g %g   ", 
	  MetaData.Time, maxdepth, HierarchyNumberOfGrids,
	  float(HierarchyMemory)/1.049e6, 
	  HierarchyAxialRatio/float(HierarchyNumberOfGrids),
	  TotalMemoryUsed, TotalMemoryDeclared);
  for (level = 0; level <= MaximumRefinementLevel; level++)
    fprintf(fptr, "%d %g %g %g %g %d   ", 
	    NumberOfGrids[level], float(Memory[level])/1.049e6, 
	    Coverage[level], MeanAxialRatio[level], FractionFlagged[level],
	    Cells[level]);
  fprintf(fptr, "\n");

  return SUCCESS;
}

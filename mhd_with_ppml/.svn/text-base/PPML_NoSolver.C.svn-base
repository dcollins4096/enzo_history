

#include <math.h>
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

#ifdef SIB2
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#endif
int PPML_NoSolver(HierarchyEntry **Grids, int NumberOfGrids,int CycleCount,
		 int NumberOfSubgrids[],fluxes ***SubgridFluxesEstimate, int level,
		 TopGridData *MetaData,ExternalBoundary *Exterior,LevelHierarchyEntry *Level
#ifdef SIB2
		 ,SiblingGridList SiblingList[]
#endif //SIB2
		 ){

  fprintf(stderr," Not solving hydro equations\n");

  int grid1;

  //
  // Allocate subgrid fluxes for each Parent Grid on this level.  
  // This might get moved, but I don't think it makes a difference.
  //

  for( grid1=0; grid1<NumberOfGrids; grid1++){
    Grids[grid1]->GridData->AllocateFluxes(NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1]);
  }

}

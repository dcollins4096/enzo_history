

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "DaveTools.h"
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
int PPML_Wrapper(HierarchyEntry **Grids, int NumberOfGrids,int CycleCount,
		 int NumberOfSubgrids[],fluxes ***SubgridFluxesEstimate, int level,
		 TopGridData *MetaData,ExternalBoundary *Exterior,LevelHierarchyEntry *Level
#ifdef SIB2
		 ,SiblingGridList SiblingList[]
#endif //SIB2
		 ){


  fprintf(stderr," === PPML_Wrapper N %d L %d P %d\n",CycleCount,level, MyProcessorNumber);
  int grid1;

  for( grid1 = 0; grid1 < NumberOfGrids; grid1++){
    if( MidWayDumpCheck( 82 ) == TRUE ){
      char basename[30];
      sprintf(basename, "data82%d%d.grid",CycleCount, level);
      FILE *dummy = fopen(basename, "a");    
      Grids[grid1]->GridData->WriteGrid(dummy, basename, grid1+1);
      fclose(dummy);
    }
  }


  //
  // Allocate subgrid fluxes for each Parent Grid on this level.  
  // This might get moved, but I don't think it makes a difference.
  //

  JBMEM_MESSAGE(MyProcessorNumber, "Before AllocateFlux");
  for( grid1=0; grid1<NumberOfGrids; grid1++){
    Grids[grid1]->GridData->AllocateFluxes(NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1]);
  }
  JBMEM_MESSAGE(MyProcessorNumber, "After  AllocateFlux");    

  //<dbg>
  //fprintf(stderr,"EJECT FROM PPML_WRAPPER\n");
  //return SUCCESS;
  //</dbg>
  //
  // Monot 1 loop: Hyun monotonicity constraint.
  //
  if( PPML_SlopeLimiter[0] == 1 ){
    JBMEM_MESSAGE(MyProcessorNumber, "Before Mono1");
    for( grid1=0; grid1<NumberOfGrids; grid1++){
      if( Grids[grid1]->GridData->PPML_Mono1() == FAIL ){
	fprintf(stderr, "PPML_Wrapper: Error in Mono1\n");
	return FAIL;
      }
    }
    JBMEM_MESSAGE(MyProcessorNumber, "After Mono1");
  }

  //
  //SBC
  //

#ifdef SIB2
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, Level ) == FAIL)
      return FAIL;
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, Level ) == FAIL)
      return FAIL;
#endif

  //
  // Monot 2 loop: Bath Monotonicity constraint.
  //
    if( PPML_SlopeLimiter[1] == 1 ){
      JBMEM_MESSAGE(MyProcessorNumber, "Before Mono2");
      for( grid1=0; grid1<NumberOfGrids; grid1++){
	Grids[grid1]->GridData->PPML_Mono2();
      }
      JBMEM_MESSAGE(MyProcessorNumber, "After  Mono2");
    }
  //
  //SBC
  //

#ifdef SIB2
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, Level ) == FAIL)
      return FAIL;
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, Level ) == FAIL)
      return FAIL;
#endif


    //
    // Monot 3 loop: PPM Monotonicity constraint.
    //
    if( PPML_SlopeLimiter[2] == 1 ){
      JBMEM_MESSAGE(MyProcessorNumber, "Before Mono3");
      for( grid1=0; grid1<NumberOfGrids; grid1++){
	Grids[grid1]->GridData->PPML_Mono3();
      }
      JBMEM_MESSAGE(MyProcessorNumber, "After Mono3");
    }
    //
    //SBC
    //
    

#ifdef SIB2
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, Level ) == FAIL)
      return FAIL;
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, Level ) == FAIL)
      return FAIL;
#endif



  //
  // Flux Computation,
  // Flux Difference,
  // CT.
  //
  JBMEM_MESSAGE(MyProcessorNumber, "Before PPML Update");
  for( grid1=0; grid1<NumberOfGrids; grid1++){
    Grids[grid1]->GridData->PPML_Update(CycleCount, level, grid1);
  }
  JBMEM_MESSAGE(MyProcessorNumber, "After PPML Update");


  return SUCCESS;
}

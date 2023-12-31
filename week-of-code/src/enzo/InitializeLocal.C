/***********************************************************************
/
/  INITIALIZE LOCAL COMPONENTS OF A NEW SIMULATION
/
/  written by: Daniel Reynolds
/  date:       April 2006
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes the local components of a new simulation 
 
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
#include "StarParticleData.h"
 
// Function prototypes
 
int CosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);
 
int NestedCosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                          TopGridData &MetaData);
 
int TurbulenceSimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);
 
 
 

int InitializeLocal(int restart, HierarchyEntry &TopGrid, TopGridData &MetaData)
{

  // Call local problem initializer
  if (debug)
    printf("InitializeLocal: Starting problem initialization.\n");

  // For problem 30 if starting from scratch, using ParallelGridIO,
  // read in data only after partitioning the grid
 
  if (!restart) {
    if (debug)
      if (ParallelRootGridIO == TRUE && ProblemType == 30) {
	if (PartitionNestedGrids) {
	  printf("InitializeLocal: Re-initialize NestedCosmologySimulation\n");
	} else {
	  printf("InitializeLocal: Re-initialize CosmologySimulation\n");
	}
      }
    
    if (ParallelRootGridIO == TRUE && ProblemType == 30) {
      if (PartitionNestedGrids) {
	if (NestedCosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
	  fprintf(stderr, "Error in NestedCosmologySimulationReInitialize.\n");
	  ENZO_FAIL("");
	}
      } else {
	if (CosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
	  fprintf(stderr, "Error in CosmologySimulationReInitialize.\n");
	  ENZO_FAIL("");
	}
      }
    }
  }
 
  // For problem 60 if starting from scratch, using ParallelGridIO, 
  // read in data only after partitioning grid.

  if (!restart) {
    if (ParallelRootGridIO == TRUE && ProblemType == 60)
      if (TurbulenceSimulationReInitialize(&TopGrid, MetaData) == FAIL) {
	fprintf(stderr, "Error in TurbulenceSimulationReInitialize.\n");
	ENZO_FAIL("");
      }
  }

  // Insert new problem intializer here...
 
  if (debug)
    printf("InitializeLocal: Finished problem initialization.\n");
 

  return SUCCESS;
 
}

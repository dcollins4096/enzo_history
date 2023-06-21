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
 
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#endif
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#ifdef ISO_GRAV
#include "GravityPotentialBoundary.h"
#endif
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "StarParticleData.h"
#ifdef RAD_HYDRO
#include "gFLDProblem.h"
#endif
 
// Function prototypes
 
#ifdef ISO_GRAV
int ProtostellarCollapseInitializeGravity(TopGridData &MetaData,
					  GravityPotentialBoundary &PotentialBdry);
#endif
 
int CosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);
 
int NestedCosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                          TopGridData &MetaData);
 
int TurbulenceSimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);
 
#ifdef RAD_HYDRO
int RadHydroInitializeLocal(int restart, HierarchyEntry &TopGrid, 
			    TopGridData &MetaData, gFLDProblem &FLDsolver);
#endif
 
 
 

#ifdef ISO_GRAV
#ifdef RAD_HYDRO
int InitializeLocal(int restart, HierarchyEntry &TopGrid, TopGridData &MetaData, 
		    GravityPotentialBoundary &PotentialBdry, 
		    gFLDProblem &FLDsolver)
#else
int InitializeLocal(int restart, HierarchyEntry &TopGrid, TopGridData &MetaData, 
		    GravityPotentialBoundary &PotentialBdry)
#endif
#else
#ifdef RAD_HYDRO
int InitializeLocal(int restart, HierarchyEntry &TopGrid, 
		    TopGridData &MetaData, gFLDProblem &FLDsolver)
#else
int InitializeLocal(int restart, HierarchyEntry &TopGrid, TopGridData &MetaData)
#endif
#endif
{

#ifdef ISO_GRAV
  // If SelfGravity is on and using TopGridIsolated, have 
  // PotentialBdry prepare the gravity domain
  if ((SelfGravity==TRUE) && (MetaData.GravityBoundary==TopGridIsolated)) {
    if (PotentialBdry.SetupGravityDomain(TopGrid, MetaData) == FAIL) {
      fprintf(stderr, "ERROR: Could not set up Gravity Domain\n");
      return FAIL;
    }

    /* Read Gravitational Potential Boundary condition info from file,
       only if isolating, restart set, and name provided */
    if (restart && (MetaData.GravityBoundaryRestart == 1)) {
	
      // Open parameter file
      FILE *fptr;
      if ((fptr = fopen(MetaData.GravityBoundaryName, "r")) == NULL) {
	fprintf(stderr, "Error opening boundary condition file: %s.\n",
		MetaData.GravityBoundaryName);
	MetaData.GravityBoundaryRestart = 0;
      }

      // Read restart file
      if (PotentialBdry.ReadPotentialBoundary(fptr) == FAIL) {
	fprintf(stderr, "ERROR: Could not read Gravity Potential Boundary\n");
	return FAIL;
      }
      fclose(fptr);
    }
  }
#endif


#ifdef RAD_HYDRO
  // If performing a radiation-hydrodynamics simulation (200 series), 
  // set up the FLD problem based on grid structures, etc. 
  // (assumes that DetermineParallelism has already been run)
  if ((ProblemType > 199) && (ProblemType < 300)) {
    if (FLDsolver.Initialize(TopGrid, MetaData) == FAIL) {
      fprintf(stderr,"Error Initializing FLD solver\n");
      return FAIL;
    }
    if (debug)
      fprintf(stdout,"InitializeLocal: RadHydro problem with %"ISYM" species\n",
	      FLDsolver.GetNumChemicalSpecies());
  }
#endif
 
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
	  return FAIL;
	}
      } else {
	if (CosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
	  fprintf(stderr, "Error in CosmologySimulationReInitialize.\n");
	  return FAIL;
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
	return FAIL;
      }
  }

  // For problem 61, TopGridIsolated, set up the boundary conditions 
  // on the gravitational potential only after partitioning grid
#ifdef ISO_GRAV
  if ((MetaData.GravityBoundaryRestart == 0) || (!restart)) {
    if (ProblemType == 61 && MetaData.GravityBoundary == TopGridIsolated) {
      if (ProtostellarCollapseInitializeGravity(MetaData, 
						PotentialBdry) == FAIL) {
	fprintf(stderr, "Error in ProtostellarCollapseInitializeGravity.\n");
	return FAIL;
      }
    }
  }
#endif


  // For problems in the 200 series, we set up the radiation 
  // boundary conditions 
#ifdef RAD_HYDRO
  if ((ProblemType > 199) && (ProblemType < 300)) {
    if (RadHydroInitializeLocal(restart, TopGrid, MetaData, FLDsolver) == FAIL) {
      fprintf(stderr, "Error in RadHydroInitializeLocal.\n");
      return FAIL;
    }
  }
#endif

  // Insert new problem intializer here...
 
  if (debug)
    printf("InitializeLocal: Finished problem initialization.\n");
 

  return SUCCESS;
 
}

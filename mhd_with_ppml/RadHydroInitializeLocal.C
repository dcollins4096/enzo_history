/***********************************************************************
/
/  INITIALIZE LOCAL COMPONENTS OF RADIATION-HYDRODYNAMICS-CHEMICAL 
/  KINETICS SIMULATION (BOUNDARY CONDITIONS, ETC.)
/
/  written by: Daniel Reynolds
/  date:       November 2006
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
#ifdef RAD_HYDRO
 
#include <string.h>
#include <stdio.h>
 
#include "gFLDProblem_preincludes.h"
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
#include "gFLDProblem.h"
 
// Function prototypes



int RadHydroInitializeLocal(int restart, HierarchyEntry &TopGrid, 
			    TopGridData &MetaData, gFLDProblem &FLDsolver)
{

  if (debug)
    fprintf(stdout,"Entering RadHydroInitializeLocal routine\n");

  // set up boundary conditions on radiation field
  if (restart) {

    // get radiation boundary conditions from restart file
    FILE *RadBdryFile;
    if ((RadBdryFile = fopen(FLDsolver.BoundaryFName,"r")) == NULL) {
      fprintf(stderr,"Error opening Radiation boundary condition file.\n");
      return FAIL;
    }
    if (FLDsolver.ReadBoundary(RadBdryFile) == FAIL) {
      fprintf(stderr,"Error reading Radiation boundary condition file.\n");
      return FAIL;
    }
    fclose(RadBdryFile);
   
  }
  else {

    float ZERO = 0.0;
    float ONE  = 1.0;

    // set boundary conditions based on problem type
    // (default to homogeneous Dirichlet)
    switch (ProblemType) {
      
    // constant test problem, set BCs based on input.  
    // 0 implies periodic, otherwise set to homogeneous Dirichlet
    case 200:
      if (FLDsolver.BdryType[0][0] != 0) {
	if (FLDsolver.SetupBoundary(0,0,1,&ZERO) == FAIL) {
	  fprintf(stderr,"Error setting x0 left radiation BCs.\n");
	  return FAIL;
	}
	if (FLDsolver.SetupBoundary(0,1,1,&ZERO) == FAIL) {
	  fprintf(stderr,"Error setting x0 right radiation BCs.\n");
	  return FAIL;
	}
      }
      if (FLDsolver.BdryType[1][0] != 0) {
	if (FLDsolver.SetupBoundary(1,0,1,&ZERO) == FAIL) {
	  fprintf(stderr,"Error setting x1 left radiation BCs.\n");
	  return FAIL;
	}
	if (FLDsolver.SetupBoundary(1,1,1,&ZERO) == FAIL) {
	  fprintf(stderr,"Error setting x1 right radiation BCs.\n");
	  return FAIL;
	}
      }
      if (FLDsolver.BdryType[2][0] != 0) {
	if (FLDsolver.SetupBoundary(2,0,1,&ZERO) == FAIL) {
	  fprintf(stderr,"Error setting x2 left radiation BCs.\n");
	  return FAIL;
	}
	if (FLDsolver.SetupBoundary(2,1,1,&ZERO) == FAIL) {
	  fprintf(stderr,"Error setting x2 right radiation BCs.\n");
	  return FAIL;
	}
      }
      break;

    // Streaming test problem.
    case 201:
      // set x0 left to Dirichlet value of 1.0
      if (FLDsolver.SetupBoundary(0,0,1,&ONE) == FAIL) {
	fprintf(stderr,"Error setting x0 left radiation BCs.\n");
	return FAIL;
      }
      // set x0 right to Neumann value of 0.0
      if (FLDsolver.SetupBoundary(0,1,1,&ZERO) == FAIL) {
	fprintf(stderr,"Error setting x0 right radiation BCs.\n");
	return FAIL;
      }
      break;

      
    // Insert new problem intializers here...


    default:
      // set BC on all faces to homogeneous Dirichlet
      if (FLDsolver.SetupBoundary(0,0,1,&ZERO) == FAIL) {
	fprintf(stderr,"Error setting x0 left radiation BCs.\n");
	return FAIL;
      }
      if (FLDsolver.SetupBoundary(0,1,1,&ZERO) == FAIL) {
	fprintf(stderr,"Error setting x0 right radiation BCs.\n");
	return FAIL;
      }
      if (FLDsolver.SetupBoundary(1,0,1,&ZERO) == FAIL) {
	fprintf(stderr,"Error setting x1 left radiation BCs.\n");
	return FAIL;
      }
      if (FLDsolver.SetupBoundary(1,1,1,&ZERO) == FAIL) {
	fprintf(stderr,"Error setting x1 right radiation BCs.\n");
	return FAIL;
      }
      if (FLDsolver.SetupBoundary(2,0,1,&ZERO) == FAIL) {
	fprintf(stderr,"Error setting x2 left radiation BCs.\n");
	return FAIL;
      }
      if (FLDsolver.SetupBoundary(2,1,1,&ZERO) == FAIL) {
	fprintf(stderr,"Error setting x2 right radiation BCs.\n");
	return FAIL;
      }
      break;
    }
  }
 
  if (debug)
    printf("RadHydroInitializeLocal: Finished.\n");

  return SUCCESS;
 
}

#endif

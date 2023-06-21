/***********************************************************************
/
/  COMPUTE THE POTENTIAL FIELD
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:  Daniel R. Reynolds
/  date:       February, 2006
/  modified2:
/
/  PURPOSE:
/
************************************************************************/
/***********************************************************************
/
/  COMPUTE THE GRAVITATIONAL POTENTIAL FIELD
/
/  This file extends Greg Bryan's Original code that uses an FFT-based
/  algorithm to solve for the gravitational potential under periodic 
/  boundary conditions on the root (level 0) grid.
/
/  Additional functionality has been added to allow for isolating 
/  (Dirichlet) boundary conditions on the root grid.  This solve calls
/  the HYPRE library for solution of the poisson equation.
/
/  NOTE: both approaches compute and store relevant solver information 
/  during the first call to the routine.  Neither of these routines 
/  'clean up' after themselves upon exit of the program, i.e. memory 
/  is allocated but never freed, requiring that the compiler take care 
/  of the remaining data upon program completion.  A future version of 
/  this gravity solver module may include a C++ class for the solver, 
/  which is initialized at the same point as the Enzo grids, stores 
/  all necessary information internally and privately, and is cleared 
/  prior to exit of the overall Enzo program.
/
/  written by: Greg Bryan
/  date:       January, 1998
/
/  modified1:  Daniel R. Reynolds
/  date:       August, 2005
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#endif
 
/* Original includes for Enzo */
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#ifdef ISO_GRAV
#include "GravityPotentialBoundary.h"
// #include "InexactNewtonSundials.h"
#include "InexactNewton.h"
#include "EnzoVector.h"
#include "SelfGravityProblem.h"
#endif
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
 

/* Function prototypes */
int CommunicationParallelFFT(region *InRegion, int NumberOfInRegions,
			     region **OutRegion, int *NumberOfOutRegions,
			     int DomainDim[], int Rank,
			     int direction, int TransposeOnCompletion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

#ifdef ISO_GRAV
int ComputePotentialFieldLevelZeroIso(TopGridData *MetaData, 
				      HierarchyEntry *Grids[], 
				      int NumberOfGrids, 
				      GravityPotentialBoundary *PotBdry);
#endif

int ComputePotentialFieldLevelZeroPer(TopGridData *MetaData, 
				      HierarchyEntry *Grids[], 
				      int NumberOfGrids);


/******************************************************************/
/* ComputePotentialFieldLevelZero observes the root-grid boundary */
/* conditions and calls the corresponding potential field solver. */
/******************************************************************/
#ifdef ISO_GRAV
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], 
				   int NumberOfGrids,
				   GravityPotentialBoundary *PotBdry)
{

  /* call the appropriate solver depending on the
     desired root grid boundary conditions */
  if (MetaData->GravityBoundary == TopGridPeriodic) {

    if (debug)
      fprintf(stdout, "ComputePotentialFieldLevelZero: TopGridPeriodic \n");
    /* Periodic root grid BC's */
    if (ComputePotentialFieldLevelZeroPer(MetaData, Grids, 
					  NumberOfGrids) == FAIL) {
      fprintf(stderr, "Error in ComputePotentialFieldLevelZeroPer.\n");
      return FAIL;
    }    
  }

  else if (MetaData->GravityBoundary == TopGridIsolated) {
    
    if (debug)
      fprintf(stdout, "ComputePotentialFieldLevelZero: TopGridIsolated \n");
    if (MetaData->TopGridRank == 3) {
      /* Isolating root grid BC's */
      if (ComputePotentialFieldLevelZeroIso(MetaData, Grids,
					    NumberOfGrids,
					    PotBdry) == FAIL) {
	fprintf(stderr, "Error in ComputePotentialFieldLevelZeroIso.\n");
	return FAIL;
      }
    }
    else {
      fprintf(stderr, "Error in ComputePotentialFieldLevelZero: \n");
      fprintf(stderr, "  Isolating BC's are only allowd for 3D problems.\n");
    }
  }

  else {

    /* Unknown root grid BC's */
    fprintf(stderr, "Error in ComputePotentialFieldLevelZero: \n");
    fprintf(stderr, "  only Periodic and Isolating BC's are allowed.\n");
    return FAIL;
  }  // end: if (TopGridPeriodic)
  
  return SUCCESS;
}
#else
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], 
				   int NumberOfGrids)
{
  /* call the periodic solver */
  if (ComputePotentialFieldLevelZeroPer(MetaData, Grids, 
					NumberOfGrids) == FAIL) {
    fprintf(stderr, "Error in ComputePotentialFieldLevelZeroPer.\n");
    return FAIL;
  }    
  return SUCCESS;
}
#endif

  
#ifdef ISO_GRAV
/******************************************************************/
/*  ComputePotentialFieldLevelZeroIso performs a root-grid        */
/*  potential field solver using isolating (Dirichlet) boundary   */
/*  conditions, through calls to the external HYPRE library.      */
/******************************************************************/
int ComputePotentialFieldLevelZeroIso(TopGridData *MetaData,
				      HierarchyEntry *Grids[], 
				      int NumberOfGrids,
				      GravityPotentialBoundary *PotentialBdry)
{
  /* Static declarations (for solver setup) */
  static int FirstCall = TRUE;

  /* Declarations */
  int i;
  int mygrid = -1;

  /* Figure out which grid I own (as this MPI process) */
  for (i=0; i<NumberOfGrids; i++) {
    if (MyProcessorNumber == Grids[i]->GridData->ReturnProcessorNumber()) {
      mygrid = i;
      break;
    }
  }
  if (mygrid == -1) {
    fprintf(stderr, "ERROR: process %"ISYM" could not locate its grid.\n",
	    MyProcessorNumber);
    return FAIL;
  }


#ifdef USE_MPI
  //  check that MyProcessorNumber agrees with MPI process ID
  MPI_Arg MPI_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_id);
  if (MyProcessorNumber != MPI_id) {
    fprintf(stderr, "ERROR: Enzo proc ID %"ISYM" does not match MPI ID %i\n", 
	   MyProcessorNumber, MPI_id);
    return FAIL;
  }
#endif


  /* Ensure that PotentialBdry has been set up already */
  if (PotentialBdry->AmIReadyToGo() == FALSE) {
    fprintf(stderr,"ERROR: GravityPotentialBoundary not yet set up\n");
    return FAIL;
  }
  
  /* get process parallelism information from grid */
  /*   per* gives periodicity of each dimension (1=>periodic) */
  int ZERO = 0;
  int per0 = (PotentialBdry->GetBoundaryType(0) == ZERO) ? 1 : 0;
  int per1 = (PotentialBdry->GetBoundaryType(1) == ZERO) ? 1 : 0;
  int per2 = (PotentialBdry->GetBoundaryType(2) == ZERO) ? 1 : 0;
  
  /*   store whether this proc is on a given external boundary */
  bool x0LBdry = PotentialBdry->AmIOnBoundary(0,0);
  bool x0RBdry = PotentialBdry->AmIOnBoundary(0,1);
  bool x1LBdry = PotentialBdry->AmIOnBoundary(1,0);
  bool x1RBdry = PotentialBdry->AmIOnBoundary(1,1);
  bool x2LBdry = PotentialBdry->AmIOnBoundary(2,0);
  bool x2RBdry = PotentialBdry->AmIOnBoundary(2,1);

  /* get enzo local grid dimensions (gravitating_mass_field, etc.) */
  int edim0 = Grids[mygrid]->GridData->GetGravitatingMassFieldDimension(0);
  int edim1 = Grids[mygrid]->GridData->GetGravitatingMassFieldDimension(1);
  int edim2 = Grids[mygrid]->GridData->GetGravitatingMassFieldDimension(2);
  
  /* initialize re-usable items for solver, potential field */
  static float *potvecdata;
  static SelfGravityProblem *GravProb;
  static EnzoVector *vec, *potvec;
  static InexactNewtonSolver *INSolve;  
//   static InexactNewtonSolverSundials *INSolve;  

  /* ------------------------------------------------------------------- */
  /* If this is the first time this routine has been called, then set up */
  /* the required solver data structures describing the problem.         */
  if (FirstCall) {

    // initialize potential field data
    potvecdata = new float[edim0*edim1*edim2];

    // create a SelfGravityProblem object, vector, and Newton solver
//     if (debug)
//       fprintf(stdout,"creating SelfGravityProblem object\n");
    GravProb = new SelfGravityProblem(*(Grids[mygrid]), *MetaData);

//     if (debug)
//       fprintf(stdout,"setting SelfGravityProblem BCs\n");
    /* insert boundary conditionsinto GravProb */
    float *BCptr;
    /*    grid owns part of x0L boundary */
    if ((per0 == 0) && x0LBdry) {
      BCptr = PotentialBdry->GetBoundaryValues(0,0);
      if (BCptr != NULL) {
	if (GravProb->SetupBoundary(0,0,0,BCptr) == FAIL) {
	  fprintf(stderr, "Error in GravProb->SetupBoundary(0,0)\n");
	  return FAIL;
	}
      }
      else {
	fprintf(stderr,"Error: p%i, x0L boundary pointer NULL\n",MPI_id);
	return FAIL;
      }
    }
    /*    grid owns part of x0R boundary */
    if ((per0 == 0) && x0RBdry) {
      BCptr = PotentialBdry->GetBoundaryValues(0,1);
      if (BCptr != NULL) {
	if (GravProb->SetupBoundary(0,1,0,BCptr) == FAIL) {
	  fprintf(stderr, "Error in GravProb->SetupBoundary(0,1)\n");
	  return FAIL;
	}
      }
      else {
	fprintf(stderr,"Error: p%i, x0R boundary pointer NULL\n",MPI_id);
	return FAIL;
      }
    }
    /*    grid owns part of x1L boundary */
    if ((per1 == 0) && x1LBdry) {
      BCptr = PotentialBdry->GetBoundaryValues(1,0);
      if (BCptr != NULL) {
	if (GravProb->SetupBoundary(1,0,0,BCptr) == FAIL) {
	  fprintf(stderr, "Error in GravProb->SetupBoundary(1,0)\n");
	  return FAIL;
	}
      }
      else {
	fprintf(stderr,"Error: p%i, x1L boundary pointer NULL\n",MPI_id);
	return FAIL;
      }
    }
    /*    grid owns part of x1R boundary */
    if ((per1 == 0) && x1RBdry) {
      BCptr = PotentialBdry->GetBoundaryValues(1,1);
      if (BCptr != NULL) {
	if (GravProb->SetupBoundary(1,1,0,BCptr) == FAIL) {
	  fprintf(stderr, "Error in GravProb->SetupBoundary(1,1)\n");
	  return FAIL;
	}
      }
      else {
	fprintf(stderr,"Error: p%i, x1R boundary pointer NULL\n",MPI_id);
	return FAIL;
      }
    }
    /*    grid owns part of x2L boundary */
    if ((per2 == 0) && x2LBdry) {
      BCptr = PotentialBdry->GetBoundaryValues(2,0);
      if (BCptr != NULL) {
	if (GravProb->SetupBoundary(2,0,0,BCptr) == FAIL) {
	  fprintf(stderr, "Error in GravProb->SetupBoundary(2,0)\n");
	  return FAIL;
	}
      }
      else {
	fprintf(stderr,"Error: p%i, x2L boundary pointer NULL\n",MPI_id);
	return FAIL;
      }
    }
    /*    grid owns part of x2R boundary */
    if ((per2 == 0) && x2RBdry) {
      BCptr = PotentialBdry->GetBoundaryValues(2,1);
      if (BCptr != NULL) {
	if (GravProb->SetupBoundary(2,1,0,BCptr) == FAIL) {
	  fprintf(stderr, "Error in GravProb->SetupBoundary(2,1)\n");
	  return FAIL;
	}
      }
      else {
	fprintf(stderr,"Error: p%i, x2R boundary pointer NULL\n",MPI_id);
	return FAIL;
      }
    }

    // get example vector for SelfGravityProblem
//     if (debug)
//       fprintf(stdout,"creating SelfGravityProblem example vector\n");
    vec = GravProb->ExampleVector();

    // initialize solution vector
//     if (debug)
//       fprintf(stdout,"initializing solution vector\n");
    potvec = vec->clone(&potvecdata);
    Grids[mygrid]->GridData->SetPotentialField(potvecdata);
  

    // create solver object
//     if (debug)
//       fprintf(stdout,"creating InexactNewtonSolver object\n");
    INSolve = new InexactNewtonSolver(vec);
    INSolve->SetMaxIters(10);
    INSolve->SetInexactNewton(0, 1.0e-8, 1.0);
    INSolve->SetNewtonTolerance(1.0e-4);
    INSolve->SetNewtonNorm(0);
    INSolve->SetMinLinesearch(1.0e-10);
//     INSolve = new InexactNewtonSolverSundials(vec);
//     if (debug)
//       fprintf(stdout,"setting solver parameters\n");
//     INSolve->SetMaxIters(10);
//     INSolve->SetInexactNewton(3);
//     INSolve->SetConstantForce(1.0e-7);
//     INSolve->SetResidualTolerance(1.0e-4);
//     INSolve->SetScaledStepTolerance(1.0e-10);


    // unset FirstCall flag
    FirstCall = FALSE;

  } // end: if (FirstCall)


  // set up rhs for this time step
//   if (debug)
//     fprintf(stdout,"setting up SelfGravityProblem\n");
  if (GravProb->Setup(*(Grids[mygrid])) == FAIL) {
    fprintf(stderr, "Error in GravProb->Setup\n");
    return FAIL;
  }
  
  // set up initial guess/solution vector
//   if (debug)
//     fprintf(stdout,"setting up initial guess/solution vector\n");
  potvec->constant(0.0);
  
  // have INSolve solve the self-gravity problem
//   if (debug)
//     fprintf(stdout,"solving the self-gravity problem\n");
  if (INSolve->Solve(GravProb,potvec) == FAIL) {
    fprintf(stderr,"ERROR: INSolve failure\n");
    return FAIL;
  }
//   else if (debug)  fprintf(stdout,"INSolve success!\n");
  
//   // TEMPORARY: output solution to file
//   if (debug)  
//     fprintf(stdout,"Writing solution to file selfgrav_sol.vec\n");
//   potvec->writeall("selfgrav_sol.vec",0);
//   fflush(stdout);
  
  
  return SUCCESS;
}
#endif



/******************************************************************/
/*  ComputePotentialFieldLevelZeroPer performs a root-grid        */
/*  potential field solver using periodic boundary conditions,    */
/*  via an FFT-based solution strategy.  This solver just calls   */
/*  the pre-existing code that Greg Bryan wrote.                  */
/******************************************************************/
int ComputePotentialFieldLevelZeroPer(TopGridData *MetaData,
				      HierarchyEntry *Grids[], 
				      int NumberOfGrids)
{
  /* Static declarations (for Green's function). */
 
  static int FirstCall = TRUE, NumberOfGreensRegions;
  static region *GreensRegion;
 
  /* Declarations. */
 
  region *OutRegion = NULL;
  int NumberOfOutRegions, DomainDim[MAX_DIMENSION];
  int i, j, n, grid, grid2;
 
  /* Allocate space for grid info. */
 
  int NumberOfRegions = NumberOfGrids;
  region *InitialRegion = new region[NumberOfRegions];
 
  /* Compute adot/a at time = t+1/2dt (time-centered). */
 
  FLOAT a = 1, dadt, MidTime = Grids[0]->GridData->ReturnTime() +
                           0.5*Grids[0]->GridData->ReturnTimeStep();
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(MidTime, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
 
  /* ------------------------------------------------------------------- */
  /* If this is the first time this routine has been called, then generate
     the Green's function. */
 
  if (FirstCall) {
 
    if (MetaData->GravityBoundary == TopGridPeriodic) {
 
      /* Periodic -- Prepare in k-space. */
 
      NumberOfGreensRegions = NumberOfGrids;
      GreensRegion = new region[NumberOfGreensRegions];
      for (grid = 0; grid < NumberOfGrids; grid++)
	if (Grids[grid]->GridData->PreparePeriodicGreensFunction(
					     &(GreensRegion[grid])) == FAIL) {
	  fprintf(stderr, "Error in grid->PreparePeriodicGreensFunction.\n");
	  return FAIL;
	}
 
    } else {
 
      fprintf(stderr, "Isolated BC's not yet implemented.\n");
      return FAIL;
 
#ifdef UNUSED
 
      for (grid = 0; grid < NumberOfGrids; grid++) {
 
	/* Generate Green's function in real space (doesn't work!). */
	
	if (Grids[grid]->GridData->PrepareGreensFunction() == FAIL) {
	  fprintf(stderr, "Error in grid->PrepareGreensFunction.\n");
	  return FAIL;
	}
 
	/* Turn it into regions. */
	
	if (Grids[grid]->GridData->PrepareFFT(&InitialRegion[grid],
					      POTENTIAL_FIELD, DomainDim)
	    == FAIL) {
	  fprintf(stderr, "Error in grid->PrepareFFT.\n");
	  return FAIL;
	}
 
      } // end loop over grids
 
      /* Forward FFT Green's function. */
 
      if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
				   &GreensRegion, &NumberOfGreensRegions,
				   DomainDim, MetaData->TopGridRank,
				   FFT_FORWARD, FALSE) == FAIL) {
	fprintf(stderr, "Error in CommunicationParallelFFT.\n");
	return FAIL;
      }
 
#endif /* UNUSED */
 
    } // end: if (Periodic)
 
    FirstCall = FALSE;
 
  } // end: if (FirstCall)
 
  /* ------------------------------------------------------------------- */
  /* Generate FFT regions for density field. */
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->PrepareFFT(&InitialRegion[grid],
					  GRAVITATING_MASS_FIELD, DomainDim)
	== FAIL) {
      fprintf(stderr, "Error in grid->PrepareFFT.\n");
      return FAIL;
    }
 
  /* Forward FFT density field. */
 
  if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
			       &OutRegion, &NumberOfOutRegions,
			       DomainDim, MetaData->TopGridRank,
			       FFT_FORWARD, TRUE) == FAIL) {
    fprintf(stderr, "Error in CommunicationParallelFFT.\n");
    return FAIL;
  }
 
  /* Quick error check. */
 
  if (NumberOfOutRegions != NumberOfGreensRegions) {
    fprintf(stderr, "OutRegion(%"ISYM") != GreensRegion(%"ISYM")\n", NumberOfOutRegions,
	    NumberOfGreensRegions);
    return FAIL;
  }
 
  /* Compute coefficient for Greens function. */
 
  float coef = GravitationalConstant/a;
  //  for (int dim = 0; dim < MetaData->TopGridRank; dim++)
  //    coef *= (DomainRightEdge[dim] - DomainLeftEdge[dim])/float(DomainDim[dim]);
			
  /* Multiply density by Green's function to get potential. */
 
  for (i = 0; i < NumberOfGreensRegions; i++)
    if (OutRegion[i].Data != NULL) {
      int size = OutRegion[i].RegionDim[0]*OutRegion[i].RegionDim[1]*
	         OutRegion[i].RegionDim[2];
      for (n = 0, j = 0; j < size; j += 2, n++) {
	OutRegion[i].Data[j  ] *= coef*GreensRegion[i].Data[n];
	OutRegion[i].Data[j+1] *= coef*GreensRegion[i].Data[n];
      }
    }
 
  /* Inverse FFT potential field. */
 
  if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
			       &OutRegion, &NumberOfOutRegions,
			       DomainDim, MetaData->TopGridRank,
			       FFT_INVERSE, TRUE) == FAIL) {
    fprintf(stderr, "Error in CommunicationParallelFFT.\n");
    return FAIL;
  }
 
  /* Copy Potential in active region into while grid. */
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->FinishFFT(&InitialRegion[grid], POTENTIAL_FIELD,
			       DomainDim) == FAIL) {
      fprintf(stderr, "Error in grid->FinishFFT.\n");
      return FAIL;
    }
 
  /* Update boundary regions of potential
     (first set BCTempL/R which are fluid BC's because that's the format
      that CheckForOverlap takes). */
 
  boundary_type BCTempLeft[MAX_DIMENSION], BCTempRight[MAX_DIMENSION];
  if (Grids[0]->GridData->ReturnGravityBoundaryType() == TopGridPeriodic) {
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      BCTempLeft[dim] = BCTempRight[dim] = periodic;
  } else {
    fprintf(stderr, "Only periodic gravity BC's allowed.\n");
    return FAIL;
  }
 
  MPI_Barrier(MPI_COMM_WORLD);
  CommunicationDirection = COMMUNICATION_SEND;
  for (grid = 0; grid < NumberOfGrids; grid++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				      BCTempLeft, BCTempRight,
     	                              &grid::CopyPotentialField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	return FAIL;
      }
 
  MPI_Barrier(MPI_COMM_WORLD);
  CommunicationDirection = COMMUNICATION_RECEIVE;
  for (grid = 0; grid < NumberOfGrids; grid++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				      BCTempLeft, BCTempRight,
     	                              &grid::CopyPotentialField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	return FAIL;
      }
 
  MPI_Barrier(MPI_COMM_WORLD);
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  /* Clean up. */
 
  delete [] InitialRegion;
  if (OutRegion != InitialRegion)
    delete [] OutRegion;
 
  if (CopyGravPotential)
    for (grid = 0; grid < NumberOfGrids; grid++)
    {
      fprintf(stderr, "Call CP from ComputePotentialFieldLevelZero\n");
      Grids[grid]->GridData->CopyPotentialToBaryonField();
    }
 
  return SUCCESS;
}

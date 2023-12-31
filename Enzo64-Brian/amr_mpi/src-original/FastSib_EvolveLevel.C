/***********************************************************************
/
/  EVOLVE LEVEL FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  February, 1995 by GB
/              Overhauled to make sure that all the subgrid's of a grid
/              advance with in lock step (i.e. with the same timestep and
/              in order).  This was done to allow a subgrid to get it's
/              boundary values from another subgrid (with the same parent).
/              Previously, a subgrid' BVs were always interpolated from its
/              parent.
/  modified2:  August, 1995 by GB
/                1) All grids on a level are processed at the same time
/                 (rather than all the subgrids of one parent).
/                2) C routines are called to loop over subgrids
/                 (so parallelizing C compilers can be used).
/                3) Subgrid timesteps are not constant over top grid step.
/              June, 1999 by GB -- Clean up somewhat
/
/  modified3:  August, 2001 by Alexei Kritsuk
/                Added 2nd call of PrepareDensityField() to compute
/                grav. potential (to be written with other baryon fields).
/  modified4:  January, 2004 by Alexei Kritsuk
/                Added support for RandomForcing
/
/  PURPOSE:
/    This routine is the main grid evolution function.  It assumes that the
/    grids of level-1 have already been advanced by dt (passed
/    in the argument) and that their boundary values are properly set.
/    We then perform a complete update on all grids on level, including:
/       - for each grid: set the boundary values from parent/subgrids
/       - for each grid: get a list of its subgrids
/       - determine the timestep based on the minimum timestep for all grids
/       - subcycle over the grid timestep and for each grid:
/           - copy the fields to the old fields
/           - solve the hydro equations (and save fluxes around subgrid)
/           - set the boundary values from parent and/or other grids
/           - update time and check dt(min) for that grid against dt(cycle)
/           - call EvolveLevel(level+1)
/           - accumulate flux around this grid
/       - correct the solution on this grid based on subgrid solutions
/       - correct the solution on this grid based on improved subgrid fluxes
/
/    This routine essentially solves (completely) the grids of this level
/       and all finer levels and then corrects the solution of
/       grids on this level based on the improved subgrid solutions/fluxes.
/
/    Note: as a convenience, we store the current grid's fluxes (i.e. the
/          fluxes around the exterior of this grid) as the last entry in
/          the list of subgrids.
/
************************************************************************/
 
#include <stdlib.h>
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
#include "LevelHierarchy.h"
 
/* function prototypes */
 
void DeleteFluxes(fluxes *Fluxes);
int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid,
			     int NumberOfGrids, int level, float dt);
float CommunicationMinValue(float Value);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
 
#ifdef SIB3
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData, FLOAT When);
#else
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
                        int level, TopGridData *MetaData, FLOAT When);
#endif
 
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

#ifdef SIB2
int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    SiblingGridList SiblingList[],
			    int level, TopGridData *MetaData,
			    ExternalBoundary *Exterior,
			    LevelHierarchyEntry * Level,
			    int CycleNumber);
#else
int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    int level, TopGridData *MetaData, 
			    ExternalBoundary *Exterior,
			    LevelHierarchyEntry * Level,
			    int CycleNumber);
#endif

#ifdef FLUX_FIX
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[],
			 LevelHierarchyEntry *SUBlingList[],
			 TopGridData *MetaData);
#else
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[]);
#endif
 
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);
int WriteMovieData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime);
int WriteTracerParticleData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime);
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float * norm, float * pTopGridTimeStep);
 
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
 
#ifdef FLUX_FIX
int CreateSUBlingList(TopGridData *MetaData,
		      HierarchyEntry *Grids[],
		      int NumberOfGrids,
		      LevelHierarchyEntry ***SUBlingList);
int DeleteSUBlingList(int NumberOfGrids,
		      LevelHierarchyEntry **SUBlingList);
#endif

void my_exit(int status);
 
 
static int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
double LevelWallTime[MAX_DEPTH_OF_HIERARCHY];
double LevelZoneCycleCount[MAX_DEPTH_OF_HIERARCHY];
double LevelZoneCycleCountPerProc[MAX_DEPTH_OF_HIERARCHY];
 
static float norm = 0.0;            //AK
static float TopGridTimeStep = 0.0; //AK
 
 
/* EvolveGrid function */
 
int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior)
{
  /* Declarations */

  int dbx = 0;
 
  FLOAT When;
  float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid;
  int RefinementFactors[MAX_DIMENSION];
  int cycle = 0, counter = 0, grid1, subgrid, grid2;
  HierarchyEntry *NextGrid;
 
#if defined(USE_JBPERF) && defined(JB_PERF_LEVELS)
  Eint32 jb_level = level;
  jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

  /* Create an array (Grids) of all the grids. */

  JBPERF_START("evolve-level-01"); // GenerateGridArray ()

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];

  JBPERF_STOP("evolve-level-01"); // GenerateGridArray ()

  JBPERF_START("evolve-level-02"); // SetBoundaryConditions()

#ifdef FLUX_FIX
  /* Create a SUBling list of the subgrids */
 
  LevelHierarchyEntry **SUBlingList;
#endif

  /* Initialize the chaining mesh used in the FastSiblingLocator. */

  if (dbx) fprintf(stderr, "EL: Initialize FSL \n"); 
  ChainingMeshStructure ChainingMesh;
  FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
			       MetaData->TopGridDims);
  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];
 
  /* Add all the grids to the chaining mesh. */

  if (dbx) fprintf(stderr, "EL: FSL AddGrid entry \n");
#ifdef DC_OPT_SIBSUB
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) 
    Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh,
						      Grids[grid1]);
#else
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);
#endif
  if (dbx) fprintf(stderr, "EL: FSL AddGrid exit \n");
 
  /* For each grid, get a list of possible siblings from the chaining mesh. */
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (Grids[grid1]->GridData->FastSiblingLocatorFindSiblings(
                              &ChainingMesh, &SiblingList[grid1],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition) == FAIL) {
      fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
      return FAIL;
    }
 
  /* Clean up the chaining mesh. */
 
  FastSiblingLocatorFinalize(&ChainingMesh);
 
  /* ================================================================== */
  /* For each grid: a) interpolate boundaries from its parent.
                    b) copy any overlapping zones.  */
 
#ifdef SIB2
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			    level, MetaData, Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                            Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#endif
 
  JBPERF_STOP("evolve-level-02"); // SetBoundaryConditions()

  /* Clear the boundary fluxes for all Grids (this will be accumulated over
     the subcycles below (i.e. during one current grid step) and used to by the
     current grid to correct the zones surrounding this subgrid (step #18). */
 
  JBPERF_START("evolve-level-03"); // ClearBoundaryFluxes()

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->ClearBoundaryFluxes();
 
  JBPERF_STOP("evolve-level-03"); // ClearBoundaryFluxes()

  /* ================================================================== */
  /* Loop over grid timesteps until the elapsed time equals the timestep
     from the level above (or loop once for the top level). */
 
  while (dtThisLevelSoFar < dtLevelAbove) {
 
    /* Determine the timestep for this iteration of the loop. */
 
    JBPERF_START("evolve-level-04"); // SetTimeStep()

    if (level == 0) {
 
      /* For root level, use dtLevelAbove. */
 
      dtThisLevel      = dtLevelAbove;
      dtThisLevelSoFar = dtLevelAbove;
 
    } else {
 
      /* Compute the mininum timestep for all grids. */
 
      dtThisLevel = huge_number;
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
	dtGrid      = Grids[grid1]->GridData->ComputeTimeStep();
	dtThisLevel = min(dtThisLevel, dtGrid);
      }
      dtThisLevel = CommunicationMinValue(dtThisLevel);
 
      /* Advance dtThisLevelSoFar (don't go over dtLevelAbove). */
 
      if (dtThisLevelSoFar+dtThisLevel*1.05 >= dtLevelAbove) {
	dtThisLevel      = dtLevelAbove - dtThisLevelSoFar;
	dtThisLevelSoFar = dtLevelAbove;
      }
      else
	dtThisLevelSoFar += dtThisLevel;
 
    }
    if (debug) printf("Level[%"ISYM"]: dt = %"GSYM"(%"GSYM"/%"GSYM")\n", level, dtThisLevel,
		      dtThisLevelSoFar, dtLevelAbove);
 
    /* Set all grid's timestep to this minimum dt. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->SetTimeStep(dtThisLevel);
 
    JBPERF_STOP("evolve-level-04"); // SetTimeStep()

    /* For each grid, compute the number of it's subgrids. */
 
    JBPERF_START("evolve-level-05"); // compute number of subgrids

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      NextGrid = Grids[grid1]->NextGridNextLevel;
      counter = 0;
      while (NextGrid != NULL) {
	NextGrid = NextGrid->NextGridThisLevel;
	if (++counter > MAX_NUMBER_OF_SUBGRIDS) {
	  fprintf(stderr, "More subgrids than MAX_NUMBER_OF_SUBGRIDS.\n");
	  return FAIL;
	}
      }
      NumberOfSubgrids[grid1] = counter + 1;
    }
 
    JBPERF_STOP("evolve-level-05"); // compute number of subgrids

    /* For each grid, create the subgrid list. */
 
    JBPERF_START("evolve-level-06"); // create subgrid list

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Allocate the subgrid fluxes for this grid. */
 
      SubgridFluxesEstimate[grid1] = new fluxes *[NumberOfSubgrids[grid1]];
 
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid1]; subgrid++)
	SubgridFluxesEstimate[grid1][subgrid] = NULL;
 
      /* Collect the flux data and store it in the newly minted fluxes.
	 Or rather that's what we should do.  Instead, we create fluxes one
	 by one in this awkward array of pointers to pointers.  This should be
	 changed so that all the routines take arrays of flux rather than
	 arrays of pointers to flux.  Dumb. */
 
      counter = 0;

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
      if (MyProcessorNumber ==
          Grids[grid1]->GridData->ReturnProcessorNumber()) {
#endif
 
	NextGrid = Grids[grid1]->NextGridNextLevel;
	while (NextGrid != NULL) {
	  SubgridFluxesEstimate[grid1][counter] = new fluxes;
	  Grids[grid1]->GridData->ComputeRefinementFactors
	                              (NextGrid->GridData, RefinementFactors);
	  NextGrid->GridData->ReturnFluxDims
             (*(SubgridFluxesEstimate[grid1][counter++]), RefinementFactors);
	  NextGrid = NextGrid->NextGridThisLevel;
	}
 
	/* Add the external boundary of this subgrid to the subgrid list. This
	   makes it easy to keep adding up the fluxes of this grid, but we must
	   keep in mind that the last subgrid should be ignored elsewhere. */
 
	SubgridFluxesEstimate[grid1][counter] = new fluxes;
	Grids[grid1]->GridData->ComputeRefinementFactors
                                   (Grids[grid1]->GridData, RefinementFactors);
	Grids[grid1]->GridData->ReturnFluxDims
               (*(SubgridFluxesEstimate[grid1][counter]), RefinementFactors);

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
      }
#endif
 
    } // end loop over grids (create Subgrid list)
 
    JBPERF_STOP("evolve-level-06"); // create subgrid list

    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Enter PrepareDensityField\n", MyProcessorNumber);
 
    JBPERF_START("evolve-level-07"); // PrepareDensityField()

    When = 0.5;
 
#ifdef SIB3
    if (SelfGravity)
      if (PrepareDensityField(LevelArray, SiblingList,
			      level, MetaData, When) == FAIL) {
	fprintf(stderr, "Error in PrepareDensityField.\n");
	return FAIL;
      }
#else
    if (SelfGravity)
      if (PrepareDensityField(LevelArray, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        return FAIL;
      }
#endif
 
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Exit PrepareDensityField\n", MyProcessorNumber);
 
    /* Prepare normalization for random forcing. Involves top grid only. */
 
    if (RandomForcing && MetaData->CycleNumber > 0 && level == 0)
      if ( ComputeRandomForcingNormalization(LevelArray, 0, MetaData,
                                             &norm, &TopGridTimeStep)
           == FAIL ) {
        fprintf(stderr, "Error in ComputeRandomForcingNormalization.\n");
        return FAIL;
      }
 
    JBPERF_STOP("evolve-level-07"); // PrepareDensityField()

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Call analysis routines. */
 
      JBPERF_START_LOW("evolve-level-08"); // Call analysis routines

      if (ProblemType == 24)
	Grids[grid1]->GridData->SphericalInfallGetProfile(level, 1);
      if (ProblemType == 30)
	Grids[grid1]->GridData->AnalyzeTrackPeaks(level, 0);
      if (ProblemType == 27)
	if (Grids[grid1]->GridData->ReturnProcessorNumber()==MyProcessorNumber){
	  float AM[3], MeanVelocity[3], DMVelocity[3];
	  FLOAT Center[] = {0,0,0}, CenterOfMass[3], DMCofM[3];
	  Grids[grid1]->GridData->CalculateAngularMomentum(Center, AM,
			   MeanVelocity, DMVelocity, CenterOfMass, DMCofM);
	  fprintf(stdout, "level = %"ISYM" %"ISYM" %"ISYM"  Vel %"FSYM" %"FSYM" %"FSYM"  DMVel %"FSYM" %"FSYM" %"FSYM"  CofM %"PSYM" %"PSYM" %"PSYM"  DMCofM %"FSYM" %"FSYM" %"FSYM"\n",
		level, LevelCycleCount[level], grid1, MeanVelocity[0],
		MeanVelocity[1], MeanVelocity[2],
		DMVelocity[0], DMVelocity[1], DMVelocity[2],
		-CenterOfMass[0], -CenterOfMass[1], -CenterOfMass[2],
		DMCofM[0], DMCofM[1], DMCofM[2]);
	}
 
      JBPERF_STOP_LOW("evolve-level-08"); // Call analysis routines

      /* Gravity: compute acceleration field for grid and particles. */
 
      JBPERF_START("evolve-level-09"); // Compute self-gravity acceleration

      if (SelfGravity) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {
 
	  /* Compute the potential. */
 
	  if (level > 0)
	    if (Grids[grid1]->GridData->SolveForPotential(Dummy, level)
		== FAIL) {
	      fprintf(stderr, "Error in grid->SolveForPotential.\n");
	      return FAIL;
	    }
	  if (Grids[grid1]->GridData->ComputeAccelerations(level) == FAIL) {
	    fprintf(stderr, "Error in grid->ComputeAccelerations.\n");
	    return FAIL;
	  }
	}
	  /* otherwise, interpolate potential from coarser grid, which is
	     now done in PrepareDensity. */
 
      } // end: if (SelfGravity)
 
      JBPERF_STOP("evolve-level-09"); // Compute self-gravity acceleration

      /* Gravity: compute field due to preset sources. */
 
      JBPERF_START_LOW("evolve-level-10"); // ComputeAccelerationFieldExternal()

      if (UniformGravity || PointSourceGravity)
	if (Grids[grid1]->GridData->ComputeAccelerationFieldExternal() ==FAIL) {
	  fprintf(stderr,"Error in grid->ComputeAccelerationFieldExternal.\n");
	  return FAIL;
	}
 
      JBPERF_STOP_LOW("evolve-level-10"); // ComputeAccelerationFieldExternal()

      /* Check for energy conservation. */
/*
      if (ComputePotential)
	if (CheckEnergyConservation(Grids, grid, NumberOfGrids, level,
				    dtThisLevel) == FAIL) {
	  fprintf(stderr, "Error in CheckEnergyConservation.\n");
	  return FAIL;
	}
*/

    } // End of loop over grids

    //This ensures that all subgrids agree in the boundary.
    //Not a big deal for hydro, but essential for DivB = 0 in MHD runs.
    //Only called on level > 0 because the root grid is dealt with differently than SG's.

#ifdef SAB
    if ( (SelfGravity || UniformGravity || PointSourceGravity) && level > 0) {
#ifdef SIB2
      if( SetAccelerationBoundary(Grids, NumberOfGrids,
				  SiblingList,
				  level, MetaData,
				  Exterior, LevelArray[level], LevelCycleCount[level]) == FAIL ) {
	fprintf(stderr,"Error with AccelerationBoundary.\n");
	return FAIL;
      }
#else
      if( SetAccelerationBoundary(Grids, NumberOfGrids,
				  level, MetaData,
				  Exterior, LevelArray[level], LevelCycleCount[level]) == FAIL ) {
	fprintf(stderr,"Error with AccelerationBoundary.\n");
	return FAIL;
      }
#endif
    }
#endif

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
 
      JBPERF_START("evolve-level-11"); // CopyBaryonFieldToOldBaryonField()

	if (Grids[grid1]->GridData->CopyBaryonFieldToOldBaryonField() == FAIL) {
	  fprintf(stderr, "Error in grid->CopyBaryonFieldToOldBaryonField.\n");
	  return FAIL;
	}
 
      JBPERF_STOP("evolve-level-11"); // CopyBaryonFieldToOldBaryonField()

      /* Add RandomForcing fields to velocities after the copying of current
         fields to old. I also update the total energy accordingly here.
         It makes no sense to force on the very first time step. */
 
      JBPERF_START_LOW("evolve-level-12"); // AddRandomForcing()

      if (RandomForcing && MetaData->CycleNumber > 0) //AK
        if(Grids[grid1]->GridData->AddRandomForcing(&norm,
                                                   TopGridTimeStep) == FAIL)
          fprintf(stderr, "Error in AddRandomForcing.\n");
 
      JBPERF_STOP_LOW("evolve-level-12"); // AddRandomForcing()

      /* Call hydro solver and save fluxes around subgrids. */
 
      JBPERF_START("evolve-level-13"); // SolveHydroEquations()

//      fprintf(stderr, "%"ISYM": Calling Hydro\n", MyProcessorNumber);
 
      if (Grids[grid1]->GridData->SolveHydroEquations(LevelCycleCount[level],
	 NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1], level) == FAIL) {
	fprintf(stderr, "Error in grid->SolveHydroEquations.\n");
	return FAIL;
      }
 
      JBPERF_STOP("evolve-level-13"); // SolveHydroEquations()

//      fprintf(stderr, "%"ISYM": Called Hydro\n", MyProcessorNumber);
 
      /* Solve the species rate equations. */
 
//      fprintf(stderr, "%"ISYM": Calling MultiSpecies\n", MyProcessorNumber);
 
      JBPERF_START("evolve-level-14"); // SolveRateEquations()

      if (MultiSpecies)
	if (Grids[grid1]->GridData->SolveRateEquations() == FAIL) {
	  fprintf(stderr, "Error in grid->SolveRateEquations.\n");
	  return FAIL;
	}
 
      JBPERF_STOP("evolve-level-14"); // SolveRateEquations()

//      fprintf(stderr, "%"ISYM": Called MultiSpecies\n", MyProcessorNumber);
 
      /* Include radiative cooling/heating. */
 
//      fprintf(stderr, "%"ISYM": Calling RadiativeCooling\n", MyProcessorNumber);
 
      JBPERF_START("evolve-level-15"); // SolveRadiativeCooling()

      if (RadiativeCooling)
	if (Grids[grid1]->GridData->SolveRadiativeCooling() == FAIL) {
	  fprintf(stderr, "Error in grid->SolveRadiativeCooling.\n");
	  return FAIL;
	}
 
      JBPERF_STOP("evolve-level-15"); // SolveRadiativeCooling()

//      fprintf(stderr, "%"ISYM": Called RadiativeCooling\n", MyProcessorNumber);
 
      /* Update particle positions (if present). */
 
//      fprintf(stderr, "%"ISYM": Calling UpdatePP\n", MyProcessorNumber);
 
      JBPERF_START("evolve-level-16"); // UpdateParticlePositions()

      if (UpdateParticlePositions(Grids[grid1]->GridData) == FAIL) {
	fprintf(stderr, "Error in UpdateParticlePositions.\n");
	return FAIL;
      }
 
      JBPERF_STOP("evolve-level-16"); // UpdateParticlePositions()

//      fprintf(stderr, "%"ISYM": Called UpdatePP\n", MyProcessorNumber);
 
      /* Include 'star' particle creation and feedback.
         (first, set the under_subgrid field). */
 
      JBPERF_START_LOW("evolve-level-17"); // star particle creation/feedback

      if (StarParticleCreation || StarParticleFeedback) {
	Grids[grid1]->GridData->ZeroSolutionUnderSubgrid(NULL,
						 ZERO_UNDER_SUBGRID_FIELD);
	LevelHierarchyEntry *Temp2 = LevelArray[level+1];
	while (Temp2 != NULL) {
	  Grids[grid1]->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData,
					 ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }
      if (StarParticleCreation || StarParticleFeedback) {
	if (Grids[grid1]->GridData->StarParticleHandler(level) == FAIL) {
	  fprintf(stderr, "Error in grid->StarParticleWrapper");
	  return FAIL;
	}
      }
 
      JBPERF_STOP_LOW("evolve-level-17"); // star particle creation/feedback

      /* Gravity: clean up AccelerationField. */

      JBPERF_START_LOW("evolve-level-18"); // clean up AccelerationField

      // David Collins removes this for MHD Amr
      if (SelfGravity || UniformGravity || PointSourceGravity) {
	if (level != MaximumGravityRefinementLevel ||
	    MaximumGravityRefinementLevel == MaximumRefinementLevel)
	  Grids[grid1]->GridData->DeleteAccelerationField();
	Grids[grid1]->GridData->DeleteParticleAcceleration();
      }
 
      JBPERF_STOP_LOW("evolve-level-18"); // clean up AccelerationField

      /* Update current problem time of this subgrid. */
 
      JBPERF_START_LOW("evolve-level-19"); // SetTimeNextTimestep()

      Grids[grid1]->GridData->SetTimeNextTimestep();
 
      JBPERF_STOP_LOW("evolve-level-19"); // SetTimeNextTimestep()

      /* If using comoving co-ordinates, do the expansion terms now. */
 
      JBPERF_START("evolve-level-20"); // ComovingExpansionTerms()

      if (ComovingCoordinates)
	Grids[grid1]->GridData->ComovingExpansionTerms();
 
      JBPERF_STOP("evolve-level-20"); // ComovingExpansionTerms()

    }  // end loop over grids
 
    /* For each grid: a) interpolate boundaries from the parent grid.
                      b) copy any overlapping zones from siblings. */
 
    JBPERF_START("evolve-level-21"); // SetBoundaryConditions()

#ifdef SIB2
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, LevelArray[level]) == FAIL)
      return FAIL;
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, LevelArray[level]) == FAIL)
      return FAIL;
#endif

    JBPERF_STOP("evolve-level-21"); // SetBoundaryConditions()

    /* Update the star particle counters. */
 
    JBPERF_START("evolve-level-22"); // CommunicationUpdateStarParticleCount()

    if (StarParticleCreation)
      if (CommunicationUpdateStarParticleCount(Grids, MetaData,
					       NumberOfGrids) == FAIL)
	return FAIL;
 
    JBPERF_STOP("evolve-level-22"); // CommunicationUpdateStarParticleCount()

    /* Check for movie output (only check if this is bottom of hierarchy). */
 
    JBPERF_START("evolve-level-23"); // WriteMovieData()

    if (LevelArray[level+1] == NULL)
      if (LevelArray[level]->GridData->ReturnTime() >=
	  MetaData->TimeLastMovieDump + MetaData->dtMovieDump &&
	  MetaData->dtMovieDump > 0.0) {
	MetaData->TimeLastMovieDump += MetaData->dtMovieDump;
	if (WriteMovieData(MetaData->MovieDumpName,
			  MetaData->MovieDumpNumber++, LevelArray, MetaData,
			  LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	  fprintf(stderr, "Error in WriteMovieData.\n");
	  return FAIL;
	}
      }
 
    JBPERF_STOP("evolve-level-23"); // WriteMovieData()

    /* Check for tracer particle output (only if this bottom of hierarchy). */
 
    JBPERF_START("evolve-level-24"); // WriteTracerParticleData()

    if (LevelArray[level+1] == NULL)
      if (LevelArray[level]->GridData->ReturnTime() >=
	  MetaData->TimeLastTracerParticleDump +
	  MetaData->dtTracerParticleDump &&
	  MetaData->dtTracerParticleDump > 0.0) {
	MetaData->TimeLastTracerParticleDump += MetaData->dtTracerParticleDump;
	if (WriteTracerParticleData(MetaData->TracerParticleDumpName,
				    MetaData->TracerParticleDumpNumber++,
				    LevelArray, MetaData,
			  LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	  fprintf(stderr, "Error in WriteTracerParticleData.\n");
	  return FAIL;
	}
      }
 
    JBPERF_STOP("evolve-level-24"); // WriteTracerParticleData()

    /* If cosmology, then compute grav. potential for output if needed. */
 
    JBPERF_START("evolve-level-25"); // PrepareDensityField()

    if (ComovingCoordinates && SelfGravity && WritePotential) {
      CopyGravPotential = TRUE;
      When = 0.0;
 
#ifdef SIB3
      if (PrepareDensityField(LevelArray, SiblingList, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        return FAIL;
      }
#else
      if (PrepareDensityField(LevelArray, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        return FAIL;
      }
#endif
 
      CopyGravPotential = FALSE;
 
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
        int Dummy;
        if (level <= MaximumGravityRefinementLevel) {
 
          /* Compute the potential. */
 
          if (level > 0)
            if (Grids[grid1]->GridData->SolveForPotential(Dummy, level)
                == FAIL) {
              fprintf(stderr, "Error in grid->SolveForPotential.\n");
              return FAIL;
            }
          // fprintf(stderr, "Call CP from EvolveLevel\n");
          Grids[grid1]->GridData->CopyPotentialToBaryonField();
        }
        /* otherwise output empty potential field. */
 
      } //  end loop over grids
    } // if WritePotential
 
    JBPERF_STOP("evolve-level-25"); // PrepareDensityField()

    /* Check for new level output (only if this is bottom of hierarchy). */
 
    JBPERF_START("evolve-level-26"); // WriteAllData()

    if (MetaData->OutputFirstTimeAtLevel > 0 &&
	level >= MetaData->OutputFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {
      MetaData->OutputFirstTimeAtLevel = level+1;
      LevelHierarchyEntry *Temp2 = LevelArray[0];
      while (Temp2->NextGridThisLevel != NULL)
	Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
      if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior,
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in WriteMovieData.\n");
	return FAIL;
      }
    }
 
    JBPERF_STOP("evolve-level-26"); // WriteAllData()

    /* Check for stop (unpleasant to exit from here, but...). */
 
    if (MetaData->StopFirstTimeAtLevel > 0 &&
	level >= MetaData->StopFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {
      fprintf(stderr, "Stopping due to request on level %"ISYM"\n", level);
      my_exit(EXIT_SUCCESS);
    }
 
    /* For each grid, delete the GravitatingMassFieldParticles. */
 
    JBPERF_START("evolve-level-27"); // DeleteGravitatingMassFieldParticles()

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->DeleteGravitatingMassFieldParticles();
 
    JBPERF_STOP("evolve-level-27"); // DeleteGravitatingMassFieldParticles()

    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */
 
//    LevelWallTime[level] += ReturnWallTime() - time1;

    if (dbx) fprintf(stderr, "EL Level %d going to Level %d\n", level, level+1);
    if (LevelArray[level+1] != NULL)
      if (EvolveLevel(MetaData, LevelArray, level+1, dtThisLevel, Exterior)
	  == FAIL) {
	fprintf(stderr, "Error in EvolveLevel (%"ISYM").\n", level);
	return FAIL;
      }

#if defined(USE_JBPERF) && defined(JB_PERF_LEVELS)
    jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

    if (dbx) fprintf(stderr, "EL Level %d returns from Level %d\n", level, level+1);

//    time1 = ReturnWallTime();
 
    /* ------------------------------------------------------- */
    /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */
 
    JBPERF_START("evolve-level-28"); // UpdateFromFinerGrids()

#ifdef FLUX_FIX

    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
    for(int list=0; list < NumberOfGrids; list++)
      SUBlingList[list] = NULL;

 
    if (FluxCorrection) {

      /* Fill in the SUBling list */

      if (dbx) fprintf(stderr, "EL: CSL entry \n");
      if (CreateSUBlingList(MetaData, Grids,
                              NumberOfGrids, &SUBlingList) == FAIL) {
        fprintf(stderr, "Error in CreateSUBlingList.\n");
        return FAIL;
      }
      if (dbx) fprintf(stderr, "EL: CSL exit \n");
    }

/* 
    LevelHierarchyEntry *NextMonkey;
 
    for(grid1 = 0; grid1 < NumberOfGrids; grid1++){
      NextMonkey = SUBlingList[grid1];
      while (NextMonkey != NULL) {
        // fprintf(stderr, "SGcheckEL%"ISYM": SUBling[%"ISYM"]->Grid pointer %p\n",
        //         MyProcessorNumber, grid1, NextMonkey->GridData);
        NextMonkey=NextMonkey->NextGridThisLevel;
      }
    }
*/

#endif
 
#ifdef FLUX_FIX
    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate,
			     SUBlingList,
			     MetaData) == FAIL)
      return FAIL;
#else
    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate) == FAIL)
      return FAIL;
#endif

    JBPERF_STOP("evolve-level-28"); // UpdateFromFinerGrids()

    if (dbx) fprintf(stderr, "OK after UpdateFromFinerGrids \n");

#ifdef FLUX_FIX
    if ( FluxCorrection ) {
      /* Clean up SUBlings */
      if (DeleteSUBlingList( NumberOfGrids, SUBlingList ) == FAIL) {
        fprintf(stderr, "Error in DeleteSUBlingList.\n");
        return FAIL;
      }
    }
#endif

    if (dbx) fprintf(stderr, "OK after DeleteSUBlingList \n");
 
    /* ------------------------------------------------------- */
    /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
       fluxes for this subgrid .
       (Note: this must be done after CorrectForRefinedFluxes). */
 
    JBPERF_START("evolve-level-29"); // AddToBoundaryFluxes()

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
      if (MyProcessorNumber ==
          Grids[grid1]->GridData->ReturnProcessorNumber()) {
#endif
 
      if (FluxCorrection)
	if (Grids[grid1]->GridData->AddToBoundaryFluxes
	    (SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1])
	    == FAIL) {
	  fprintf(stderr, "Error in grid->AddToBoundaryFluxes.\n");
	  return FAIL;
	}
 
      /* Delete fluxes pointed to by SubgridFluxesEstimate[subgrid]. */
 
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid1]; subgrid++) {
	DeleteFluxes(SubgridFluxesEstimate[grid1][subgrid]);
	delete       SubgridFluxesEstimate[grid1][subgrid];
      }
      delete [] SubgridFluxesEstimate[grid1];

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
      }
#endif
 
    } // end of loop over grids
 
    JBPERF_STOP("evolve-level-29"); // AddToBoundaryFluxes()

    /* Recompute radiation field, if requested. */
 
    JBPERF_START("evolve-level-30"); // RadiationFieldUpdate()

    if (RadiationFieldType >= 10 && RadiationFieldType <= 11 &&
	level <= RadiationFieldLevelRecompute)
      if (RadiationFieldUpdate(LevelArray, level, MetaData) == FAIL) {
	fprintf(stderr, "Error in RecomputeRadiationField.\n");
	return FAIL;
      }
 
    JBPERF_STOP("evolve-level-30"); // RadiationFieldUpdate()

    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */
 
    JBPERF_START("evolve-level-31"); // RebuildHierarchy()

//    LevelWallTime[level] += ReturnWallTime() - time1;
    if (dtThisLevelSoFar < dtLevelAbove)
      if (RebuildHierarchy(MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
//    time1 = ReturnWallTime();
 
    /* Count up number of grids on this level. */
 
    int GridMemory, NumberOfCells, CellsTotal;
    float AxialRatio, GridVolume;
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->CollectGridInformation(GridMemory, GridVolume,
					  NumberOfCells, AxialRatio, CellsTotal);
      LevelZoneCycleCount[level] += NumberOfCells;
      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber())
	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }
 
    JBPERF_STOP("evolve-level-31"); // RebuildHierarchy()

    cycle++;
    LevelCycleCount[level]++;
 
  } // end of loop over subcycles
 
  if (debug)
    fprintf(stderr, "EvolveLevel[%"ISYM"]: NumberOfSubCycles = %"ISYM" (%"ISYM" total)\n", level,
           cycle, LevelCycleCount[level]);
 
  /* If possible & desired, report on memory usage. */
 
  //  if (debug)
  ReportMemoryUsage("Memory usage report: Evolve Level");
 
#if defined(USE_JBPERF) && defined(JB_PERF_LEVELS)
  jbPerf.attribute ("level",0,JB_NULL);
#endif

  /* Clean up. */
 
#ifdef UNUSED
  if (level > MaximumGravityRefinementLevel &&
      level == MaximumRefinementLevel)
    ZEUSQuadraticArtificialViscosity /= 1;
#endif
 
  delete [] NumberOfSubgrids;
  delete [] Grids;
  delete [] SubgridFluxesEstimate;
 
  /* Clean up the sibling list. */
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    delete [] SiblingList[grid1].GridList;
  delete [] SiblingList;

  if (dbx) fprintf(stderr, "Return from EL Level %d\n", level);
 
  return SUCCESS;
 
}

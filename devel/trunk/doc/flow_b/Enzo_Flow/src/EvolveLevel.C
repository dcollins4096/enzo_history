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
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			int level, TopGridData *MetaData);
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior);
int UpdateFromFinerGrids(HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[], 
			 fluxes **SubgridFluxesEstimate[]);
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);
int WriteMovieData(char *basename, int filenumber, 
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData, 
		   FLOAT WriteTime);
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid, 
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);
void my_exit(int status);


static int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];


/* EvolveGrid function */

int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior)
{
  /* Declarations */

  float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid;
  int RefinementFactors[MAX_DIMENSION];
  int cycle = 0, counter = 0, grid, subgrid;
  HierarchyEntry *NextGrid;

  //=================================================================

  JBPERF_INIT;

  JBPERF_START_LEVEL("EL00"); // EvolveLevel

  //=================================================================

  /* Create an array (Grids) of all the grids. */

  JBPERF_START_LEVEL("EL01"); // GenerateGridArray ()

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];

  JBPERF_STOP_LEVEL("EL01"); // Create grid array

  /* ================================================================== */
  /* For each grid: a) interpolate boundaries from its parent.
                    b) copy any overlapping zones.  */

  JBPERF_START_LEVEL("EL02"); // SetBoundaryConditions ()

  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, 
			    Exterior) == FAIL)
    return FAIL;

  JBPERF_STOP_LEVEL("EL02"); // SetBoundaryConditions ()

  //=================================================================
  /* Clear the boundary fluxes for all Grids (this will be accumulated over
     the subcycles below (i.e. during one current grid step) and used to by the
     current grid to correct the zones surrounding this subgrid (step #18). */

  JBPERF_START_LEVEL("EL03"); // ClearBoundaryFluxes ()

  for (grid = 0; grid < NumberOfGrids; grid++)
    Grids[grid]->GridData->ClearBoundaryFluxes();

  JBPERF_STOP_LEVEL("EL03"); // ClearBoundaryFluxes ()

  /* ================================================================== */
  /* Loop over grid timesteps until the elapsed time equals the timestep
     from the level above (or loop once for the top level). */

  while (dtThisLevelSoFar < dtLevelAbove) {

    //=================================================================
    
    /* Determine the timestep for this iteration of the loop. */

    JBPERF_START_LEVEL("EL04"); // SetTimeStep ()

    if (level == 0) {

      /* For root level, use dtLevelAbove. */

      dtThisLevel      = dtLevelAbove;
      dtThisLevelSoFar = dtLevelAbove;

    } else {

      /* Compute the mininum timestep for all grids. */

      dtThisLevel = huge_number;
      for (grid = 0; grid < NumberOfGrids; grid++) {
	dtGrid      = Grids[grid]->GridData->ComputeTimeStep();
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
    if (debug) printf("Level[%d]: dt = %g(%g/%g)\n", level, dtThisLevel,
		      dtThisLevelSoFar, dtLevelAbove);

    /* Set all grid's timestep to this minimum dt. */

    for (grid = 0; grid < NumberOfGrids; grid++)
      Grids[grid]->GridData->SetTimeStep(dtThisLevel);

    JBPERF_STOP_LEVEL("EL04"); // SetTimeStep ()

    //=================================================================
    /* For each grid, compute the number of it's subgrids. */

    JBPERF_START_LEVEL("EL05"); // compute number of subgrids

    for (grid = 0; grid < NumberOfGrids; grid++) {
      NextGrid = Grids[grid]->NextGridNextLevel;
      counter = 0;
      while (NextGrid != NULL) {
	NextGrid = NextGrid->NextGridThisLevel;
	if (++counter > MAX_NUMBER_OF_SUBGRIDS) {
	  fprintf(stderr, "More subgrids than MAX_NUMBER_OF_SUBGRIDS.\n");
	  return FAIL;
	}
      }
      NumberOfSubgrids[grid] = counter + 1;
    }

    JBPERF_STOP_LEVEL("EL05"); // compute number of subgrids

    //=================================================================
    /* For each grid, create the subgrid list. */

    JBPERF_START_LEVEL("EL06"); // Create subgrid list

    for (grid = 0; grid < NumberOfGrids; grid++) {

      /* Allocate the subgrid fluxes for this grid. */

      SubgridFluxesEstimate[grid] = new fluxes *[NumberOfSubgrids[grid]];

      /* Collect the flux data and store it in the newly minted fluxes.
	 Or rather that's what we should do.  Instead, we create fluxes one
	 by one in this awkward array of pointers to pointers.  This should be
	 changed so that all the routines take arrays of flux rather than
	 arrays of pointers to flux.  Dumb. */

      counter = 0;
      NextGrid = Grids[grid]->NextGridNextLevel;
      while (NextGrid != NULL) {
	SubgridFluxesEstimate[grid][counter] = new fluxes;
	Grids[grid]->GridData->ComputeRefinementFactors
	                              (NextGrid->GridData, RefinementFactors);
	NextGrid->GridData->ReturnFluxDims
             (*(SubgridFluxesEstimate[grid][counter++]), RefinementFactors);
	NextGrid = NextGrid->NextGridThisLevel;
      }

      /* Add the external boundary of this subgrid to the subgrid list. This
	 makes it easy to keep adding up the fluxes of this grid, but we
	 must keep in mind that the last subgrid should be ignored elsewhere.*/

      SubgridFluxesEstimate[grid][counter] = new fluxes;
      Grids[grid]->GridData->ComputeRefinementFactors
                                   (Grids[grid]->GridData, RefinementFactors);
      Grids[grid]->GridData->ReturnFluxDims
               (*(SubgridFluxesEstimate[grid][counter]), RefinementFactors);

    } // end loop over grids (create Subgrid list)

    JBPERF_STOP_LEVEL("EL06"); // Create subgrid list

    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */

    JBPERF_START_LEVEL("EL07"); // PrepareDensityField ()

    if (SelfGravity)
      if (PrepareDensityField(LevelArray, level, MetaData) == FAIL) {
	fprintf(stderr, "Error in PrepareDensityField.\n");
	return FAIL;
      }

    JBPERF_STOP_LEVEL("EL07"); // PrepareDensityField ()

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */

    for (grid = 0; grid < NumberOfGrids; grid++) {

      //=================================================================
      /* Call analysis routines. */

      JBPERF_START_LEVEL("EL08"); // Call analysis routines

      if (ProblemType == 24) 
	Grids[grid]->GridData->SphericalInfallGetProfile(level, 1);
      if (ProblemType == 30)
	Grids[grid]->GridData->AnalyzeTrackPeaks(level, 0);
      if (ProblemType == 27) 
	if (Grids[grid]->GridData->ReturnProcessorNumber()==MyProcessorNumber){
	  float AM[3], MeanVelocity[3], DMVelocity[3];
	  FLOAT Center[] = {0,0,0}, CenterOfMass[3], DMCofM[3];
	  Grids[grid]->GridData->CalculateAngularMomentum(Center, AM, 
			   MeanVelocity, DMVelocity, CenterOfMass, DMCofM);
	  fprintf(stdout, "level = %d %d %d  Vel %f %f %f  DMVel %f %f %f  CofM %"PSYM" %"PSYM" %"PSYM"  DMCofM %f %f %f\n",
		level, LevelCycleCount[level], grid, MeanVelocity[0],
		MeanVelocity[1], MeanVelocity[2],
		DMVelocity[0], DMVelocity[1], DMVelocity[2],
		-CenterOfMass[0], -CenterOfMass[1], -CenterOfMass[2],
		DMCofM[0], DMCofM[1], DMCofM[2]);
	}

      JBPERF_STOP_LEVEL("EL08"); // Call analysis routines

      //=================================================================
      /* Gravity: compute acceleration field for grid and particles. */

      JBPERF_START_LEVEL("EL09"); // Compute self-gravity

      if (SelfGravity) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {

	  /* Compute the potential. */

	  if (level > 0)
	    if (Grids[grid]->GridData->SolveForPotential(Dummy, level) 
		== FAIL) {
	      fprintf(stderr, "Error in grid->SolveForPotential.\n");
	      return FAIL;
	    }
	  if (Grids[grid]->GridData->ComputeAccelerations(level) == FAIL) {
	    fprintf(stderr, "Error in grid->ComputeAccelerations.\n");
	    return FAIL;
	  }
	}
	  /* otherwise, interpolate potential from coarser grid, which is
	     now done in PrepareDensity. */

      } // end: if (SelfGravity)

      JBPERF_STOP_LEVEL("EL09"); // Compute self-gravity

      //=================================================================
      /* Gravity: compute field due to preset sources. */

      JBPERF_START_LEVEL("EL10"); // Compute external gravity

      if (UniformGravity || PointSourceGravity)
	if (Grids[grid]->GridData->ComputeAccelerationFieldExternal() ==FAIL) {
	  fprintf(stderr,"Error in grid->ComputeAccelerationFieldExternal.\n");
	  return FAIL;
	}

      JBPERF_STOP_LEVEL("EL10"); // Compute external gravity

      /* Check for energy conservation. */
/*
      if (ComputePotential)
	if (CheckEnergyConservation(Grids, grid, NumberOfGrids, level, 
				    dtThisLevel) == FAIL) {
	  fprintf(stderr, "Error in CheckEnergyConservation.\n");
	  return FAIL;
	}
*/

      //=================================================================
      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
    
      JBPERF_START_LEVEL("EL11"); // CopyBaryonFieldToOldBaryonField ()

      if (Grids[grid]->GridData->CopyBaryonFieldToOldBaryonField() == FAIL) {
	fprintf(stderr, "Error in grid->CopyBaryonFieldToOldBaryonField.\n");
	return FAIL;
      }

      JBPERF_STOP_LEVEL("EL11"); // CopyBaryonFieldToOldBaryonField ()

      //=================================================================
      /* Call hydro solver and save fluxes around subgrids. */

      JBPERF_START_LEVEL("EL12"); // SolveHydroEquations ()

      if (Grids[grid]->GridData->SolveHydroEquations(LevelCycleCount[level],
	 NumberOfSubgrids[grid], SubgridFluxesEstimate[grid], level) == FAIL) {
	fprintf(stderr, "Error in grid->SolveHydroEquations.\n");
	return FAIL;
      }

      JBPERF_STOP_LEVEL("EL12"); // SolveHydroEquations ()

      //=================================================================
      /* Solve the species rate equations. */

      JBPERF_START_LEVEL("EL13"); // SolveRateEquations ()

      if (MultiSpecies)
	if (Grids[grid]->GridData->SolveRateEquations() == FAIL) {
	  fprintf(stderr, "Error in grid->SolveRateEquations.\n");
	  return FAIL;
	}

      JBPERF_STOP_LEVEL("EL13"); // SolveRateEquations ()

      //=================================================================
      /* Include radiative cooling/heating. */

      JBPERF_START_LEVEL("EL14"); // SolveRadiativeCooling ()

      if (RadiativeCooling)
	if (Grids[grid]->GridData->SolveRadiativeCooling() == FAIL) {
	  fprintf(stderr, "Error in grid->SolveRadiativeCooling.\n");
	  return FAIL;
	}

      JBPERF_STOP_LEVEL("EL14"); // SolveRadiativeCooling ()

      //=================================================================
      /* Update particle positions (if present). */

      JBPERF_START_LEVEL("EL15"); // UpdateParticlePositions ()

      if (UpdateParticlePositions(Grids[grid]->GridData) == FAIL) {
	fprintf(stderr, "Error in UpdateParticlePositions.\n");
	return FAIL;
      }

      JBPERF_STOP_LEVEL("EL15"); // UpdateParticlePositions ()

      //=================================================================
      /* Include 'star' particle creation and feedback. 
         (first, set the under_subgrid field). */

      JBPERF_START_LEVEL("EL16"); // Include star particle creation and feedback

      if (StarParticleCreation || StarParticleFeedback) {
	Grids[grid]->GridData->ZeroSolutionUnderSubgrid(NULL, 
						 ZERO_UNDER_SUBGRID_FIELD);
	LevelHierarchyEntry *Temp2 = LevelArray[level+1];
	while (Temp2 != NULL) {
	  Grids[grid]->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
					 ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }
      if (StarParticleCreation || StarParticleFeedback) {
	if (Grids[grid]->GridData->StarParticleHandler(level) == FAIL) {
	  fprintf(stderr, "Error in grid->StarParticleWrapper");
	  return FAIL;
	}
      }

      JBPERF_STOP_LEVEL("EL16"); // Include star particle creation and feedback

      //=================================================================
      /* Gravity: clean up AccelerationField. */

      JBPERF_START_LEVEL("EL17"); // clean up gravity acceleration field

      if (SelfGravity || UniformGravity || PointSourceGravity) {
	if (level != MaximumGravityRefinementLevel ||
	    MaximumGravityRefinementLevel == MaximumRefinementLevel)
	  Grids[grid]->GridData->DeleteAccelerationField();
	Grids[grid]->GridData->DeleteParticleAcceleration();
      }

      JBPERF_STOP_LEVEL("EL17"); // clean up gravity acceleration field

      //=================================================================
      /* Update current problem time of this subgrid. */

      JBPERF_START_LEVEL("EL18"); // SetNextTimestep ()

      Grids[grid]->GridData->SetTimeNextTimestep();

      JBPERF_STOP_LEVEL("EL18"); // SetNextTimestep ()

      //=================================================================
      /* If using comoving co-ordinates, do the expansion terms now. */

      JBPERF_START_LEVEL("EL19"); // ComovingExpansionTerms ()

      if (ComovingCoordinates)
	Grids[grid]->GridData->ComovingExpansionTerms();

      JBPERF_STOP_LEVEL("EL19"); // ComovingExpansionTerms ()

    }  // end loop over grids

    //=================================================================
    /* For each grid: a) interpolate boundaries from the parent grid.
                      b) copy any overlapping zones from siblings. */

    JBPERF_START_LEVEL("EL20"); // SetBoundaryConditions ()

    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, 
			      Exterior) == FAIL)
      return FAIL;

    JBPERF_STOP_LEVEL("EL20"); // SetBoundaryConditions ()

    //=================================================================
    /* Update the star particle counters. */

    JBPERF_START_LEVEL("EL21"); // CommunicationUpdateStarParticleCount ()

    if (StarParticleCreation)
      if (CommunicationUpdateStarParticleCount(Grids, MetaData,
					       NumberOfGrids) == FAIL)
	return FAIL;

    JBPERF_STOP_LEVEL("EL21"); // CommunicationUpdateStarParticleCount ()

    //=================================================================
    /* Check for movie output (only check if this is bottom of hierarchy). */

    JBPERF_START_LEVEL("EL22"); // WriteMovieData ()

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

    JBPERF_STOP_LEVEL("EL22"); // WriteMovieData ()

    //=================================================================
    /* Check for new level output (only if this is bottom of hierarchy). */

    JBPERF_START_LEVEL("EL23"); // WriteAllData ()

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

    JBPERF_STOP_LEVEL("EL23"); // WriteAllData ()

    /* Check for stop (unpleasant to exit from here, but...). */

    if (MetaData->StopFirstTimeAtLevel > 0 &&
	level >= MetaData->StopFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {
      fprintf(stderr, "Stopping due to request on level %d\n", level);
      my_exit(EXIT_SUCCESS);
    }

    //=================================================================
    /* For each grid, delete the GravitatingMassFieldParticles. */

    JBPERF_START_LEVEL("EL24"); // DeleteGravitatingMassFieldParticles ()

    for (grid = 0; grid < NumberOfGrids; grid++)
      Grids[grid]->GridData->DeleteGravitatingMassFieldParticles();

    JBPERF_STOP_LEVEL("EL24"); // DeleteGravitatingMassFieldParticles ()

    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */

    JBPERF_START_LEVEL("EL25"); // EvolveLevel recursion

    if (LevelArray[level+1] != NULL)
      if (EvolveLevel(MetaData, LevelArray, level+1, dtThisLevel, Exterior) 
	  == FAIL) {
	fprintf(stderr, "Error in EvolveLevel (%d).\n", level);
	return FAIL;
      }

    JBPERF_STOP_LEVEL("EL25"); // EvolveLevel recursion

    /* ------------------------------------------------------- */
    /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */

    JBPERF_START_LEVEL("EL26"); // UpdateFromFinerGrids ()

    if (UpdateFromFinerGrids(Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate) == FAIL)
      return FAIL;

    JBPERF_STOP_LEVEL("EL26"); // UpdateFromFinerGrids ()

    /* ------------------------------------------------------- */
    /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
       fluxes for this subgrid .
       (Note: this must be done after CorrectForRefinedFluxes). */

    JBPERF_START_LEVEL("EL27"); // Add saved fluxes to this subgrids exterior fluxes 

    for (grid = 0; grid < NumberOfGrids; grid++) {

      if (FluxCorrection)
	if (Grids[grid]->GridData->AddToBoundaryFluxes
	    (SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1])
	    == FAIL) {
	  fprintf(stderr, "Error in grid->AddToBoundaryFluxes.\n");
	  return FAIL;
	}

      /* Delete fluxes pointed to by SubgridFluxesEstimate[subgrid]. */

      for (subgrid = 0; subgrid < NumberOfSubgrids[grid]; subgrid++) {
	DeleteFluxes(SubgridFluxesEstimate[grid][subgrid]);
	delete       SubgridFluxesEstimate[grid][subgrid];
      }
      delete [] SubgridFluxesEstimate[grid];
      
    } // end of loop over grids

    JBPERF_STOP_LEVEL("EL27"); // Add saved fluxes to this subgrids exterior fluxes 

    //=================================================================
    /* Recompute radiation field, if requested. */

    JBPERF_START_LEVEL("EL28"); // RadiationFieldUpdate ()

    if (RadiationFieldType >= 10 && RadiationFieldType <= 11 && 
	level <= RadiationFieldLevelRecompute)
      if (RadiationFieldUpdate(LevelArray, level, MetaData) == FAIL) {
	fprintf(stderr, "Error in RecomputeRadiationField.\n");
	return FAIL;
      }

    JBPERF_STOP_LEVEL("EL28"); // RadiationFieldUpdate ()

    //=================================================================
    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */

    JBPERF_START_LEVEL("EL29"); // RebuildHierarchy ()

    if (dtThisLevelSoFar < dtLevelAbove)
      if (RebuildHierarchy(MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }

    JBPERF_STOP_LEVEL("EL29"); // RebuildHierarchy ()

    cycle++;
    LevelCycleCount[level]++;

  } // end of loop over subcycles

  if (debug)
    printf("EvolveLevel[%d]: NumberOfSubCycles = %d (%d total)\n", level, 
           cycle, LevelCycleCount[level]);

  /* If possible & desired, report on memory usage. */

  ReportMemoryUsage("Memory usage report: Evolve Level");

  /* Clean up. */

  delete [] NumberOfSubgrids;
  delete [] Grids;
  delete [] SubgridFluxesEstimate;

  //=================================================================
  JBPERF_STOP_LEVEL("EL00") // EvolveLevel
  //=================================================================

  return SUCCESS;

}

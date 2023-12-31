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
/  modified5:  February, 2006 by Daniel Reynolds
/                Added PotentialBdry to EvolveLevel and 
/                PrepareDensityField calls, so that it can be used
/                within computing isolating BCs for self-gravity.
/  modified6:  January, 2007 by Robert Harkness
/                Group and in-core i/o
/  modified7:  December, 2007 by Robert Harkness
/                Optional StaticSiblingList for root grid
/  modified8:  April, 2009 by John Wise
/                Added star particle class and radiative transfer
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
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "performance.h"
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
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
 
/* function prototypes */
 

int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid,
			     int NumberOfGrids, int level, float dt);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int WriteStreamData(LevelHierarchyEntry *LevelArray[], int level,
		    TopGridData *MetaData, int *CycleCount, int open=FALSE);
int CallProblemSpecificRoutines(TopGridData * MetaData, HierarchyEntry *ThisGrid,
				int GridNum, float *norm, float TopGridTimeStep, 
				int level, int LevelCycleCount[]);  //moo

#ifdef FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData, FLOAT When);
#else  // !FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
                        int level, TopGridData *MetaData, FLOAT When);
#endif  // end FAST_SIB
 
#ifdef FAST_SIB
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#endif



#ifdef SAB
#ifdef FAST_SIB
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
int CreateFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]);		 
int FinalizeFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]);		 
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);


int OutputFromEvolveLevel(LevelHierarchyEntry *LevelArray[],TopGridData *MetaData,
		      int level, ExternalBoundary *Exterior);
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float * norm, float * pTopGridTimeStep);
int CreateSiblingList(HierarchyEntry ** Grids, int NumberOfGrids, SiblingGridList *SiblingList, 
		      int StaticLevelZero,TopGridData * MetaData,int level);

#ifdef FLUX_FIX
int CreateSUBlingList(TopGridData *MetaData,
		      HierarchyEntry *Grids[],
		      int NumberOfGrids,
		      LevelHierarchyEntry ***SUBlingList);
int DeleteSUBlingList(int NumberOfGrids,
		      LevelHierarchyEntry **SUBlingList);
#endif

int StarParticleInitialize(LevelHierarchyEntry *LevelArray[], int ThisLevel,
			   TopGridData *MetaData, Star *&AllStars);
int StarParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			 int level, Star *&AllStars);
int AdjustRefineRegion(LevelHierarchyEntry *LevelArray[], 
		       TopGridData *MetaData, int EL_level);

#ifdef TRANSFER
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *AllStars, FLOAT GridTime, int level, int LoopTime = TRUE);
int RadiativeTransferPrepare(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *&AllStars,
			     float dtLevelAbove);
#endif

int SetLevelTimeStep(HierarchyEntry *Grids[],
        int NumberOfGrids, int level,
        float *dtThisLevelSoFar, float *dtThisLevel,
        float dtLevelAbove);

void my_exit(int status);
 
int CallPython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level);
 
int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
int MovieCycleCount[MAX_DEPTH_OF_HIERARCHY];
double LevelWallTime[MAX_DEPTH_OF_HIERARCHY];
double LevelZoneCycleCount[MAX_DEPTH_OF_HIERARCHY];
double LevelZoneCycleCountPerProc[MAX_DEPTH_OF_HIERARCHY];
 
static float norm = 0.0;            //AK
static float TopGridTimeStep = 0.0; //AK
#ifdef STATIC_SIBLING_LIST
static int StaticLevelZero = 1;
#else
static int StaticLevelZero = 0;
#endif

int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior)
{
  /* Declarations */

  int dbx = 0;
 
  FLOAT When, GridTime;
  //float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid, dtActual, dtLimit;
  float dtThisLevelSoFar = 0.0, dtThisLevel;
  int cycle = 0, counter = 0, grid1, subgrid, grid2;
  HierarchyEntry *NextGrid;
  int dummy_int;
 
  // Update lcaperf "level" attribute

  Eint32 jb_level = level;
#ifdef USE_JBPERF
  jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

  /* Create an array (Grids) of all the grids. */

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];

#ifdef FLUX_FIX
  /* Create a SUBling list of the subgrids */
  LevelHierarchyEntry **SUBlingList;
#endif


  /* Initialize the chaining mesh used in the FastSiblingLocator. */

  if (dbx) fprintf(stderr, "EL: Initialize FSL \n"); 
  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];
  CreateSiblingList(Grids, NumberOfGrids, SiblingList, StaticLevelZero,MetaData,level);
  
  /* On the top grid, adjust the refine region so that only the finest
     particles are included.  We don't want the more massive particles
     to contaminate the high-resolution region. */

  AdjustRefineRegion(LevelArray, MetaData, level);

  /* ================================================================== */
  /* For each grid: a) interpolate boundaries from its parent.
                    b) copy any overlapping zones.  */
 
#ifdef FAST_SIB
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			    level, MetaData, Exterior, LevelArray[level]) == FAIL)
    ENZO_FAIL("");
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                            Exterior, LevelArray[level]) == FAIL)
    ENZO_FAIL("");
#endif
 


  /* Clear the boundary fluxes for all Grids (this will be accumulated over
     the subcycles below (i.e. during one current grid step) and used to by the
     current grid to correct the zones surrounding this subgrid (step #18). */
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->ClearBoundaryFluxes();
 
  /* After we calculate the ghost zones, we can initialize streaming
     data files (only on level 0) */

  if (MetaData->FirstTimestepAfterRestart == TRUE && level == 0)
    WriteStreamData(LevelArray, level, MetaData, MovieCycleCount);

  /* ================================================================== */
  /* Loop over grid timesteps until the elapsed time equals the timestep
     from the level above (or loop once for the top level). */
 
  while (dtThisLevelSoFar < dtLevelAbove) {
 
    SetLevelTimeStep(Grids, NumberOfGrids, level, 
        &dtThisLevelSoFar, &dtThisLevel, dtLevelAbove);

    /* Streaming movie output (write after all parent grids are
       updated) */

    WriteStreamData(LevelArray, level, MetaData, MovieCycleCount);

    /* Initialize the star particles */

    Star *AllStars = NULL;
    StarParticleInitialize(LevelArray, level, MetaData, AllStars);

    /* Initialize the radiative transfer */

#ifdef TRANSFER
    RadiativeTransferPrepare(LevelArray, level, MetaData, AllStars, 
			     dtLevelAbove);
#endif /* TRANSFER */
 
    CreateFluxes(Grids,SubgridFluxesEstimate,NumberOfGrids,NumberOfSubgrids);

    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */

    When = 0.5;
 
#ifdef FAST_SIB
     PrepareDensityField(LevelArray, SiblingList, level, MetaData, When);
#else   // !FAST_SIB
     PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB
 
 
    /* Prepare normalization for random forcing. Involves top grid only. */
 
    ComputeRandomForcingNormalization(LevelArray, 0, MetaData,
				      &norm, &TopGridTimeStep);

    /* Solve the radiative transfer */
	
#ifdef TRANSFER
    GridTime = Grids[0]->GridData->ReturnTime();// + dtThisLevel;
    EvolvePhotons(MetaData, LevelArray, AllStars, GridTime, level);
#endif /* TRANSFER */
 
    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      CallProblemSpecificRoutines(MetaData, Grids[grid1], grid1, &norm, 
				  TopGridTimeStep, level, LevelCycleCount);

      /* Gravity: compute acceleration field for grid and particles. */
 
      if (SelfGravity) {
	if (level <= MaximumGravityRefinementLevel) {
 
	  /* Compute the potential. */
 
	  if (level > 0)
	    Grids[grid1]->GridData->SolveForPotential(level);
	  Grids[grid1]->GridData->ComputeAccelerations(level);
	}
	  /* otherwise, interpolate potential from coarser grid, which is
	     now done in PrepareDensity. */
 
      } // end: if (SelfGravity)
 
      /* Gravity: compute field due to preset sources. */
 
      Grids[grid1]->GridData->ComputeAccelerationFieldExternal();
 
      /* Radiation Pressure: add to acceleration field */

#ifdef TRANSFER
      Grids[grid1]->GridData->AddRadiationPressureAcceleration();
#endif /* TRANSFER */

      /* Check for energy conservation. */
/*
      if (ComputePotential)
	if (CheckEnergyConservation(Grids, grid, NumberOfGrids, level,
				    dtThisLevel) == FAIL) {
	  fprintf(stderr, "Error in CheckEnergyConservation.\n");
	  ENZO_FAIL("");
	}
*/
#ifdef SAB
    } // End of loop over grids
    
    //Ensure the consistency of the AccelerationField
    SetAccelerationBoundary(Grids, NumberOfGrids,SiblingList,level, MetaData,
			    Exterior, LevelArray[level], LevelCycleCount[level]);
    
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
#endif //SAB.
      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
 
      Grids[grid1]->GridData->CopyBaryonFieldToOldBaryonField();

      /* Call hydro solver and save fluxes around subgrids. */

      Grids[grid1]->GridData->SolveHydroEquations(LevelCycleCount[level],
	    NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1], level);

      /* Solve the cooling and species rate equations. */
 
      Grids[grid1]->GridData->MultiSpeciesHandler();

      /* Update particle positions (if present). */
 
      UpdateParticlePositions(Grids[grid1]->GridData);

      /* Include 'star' particle creation and feedback. */

      Grids[grid1]->GridData->StarParticleHandler
	(Grids[grid1]->NextGridNextLevel, level);
 
      /* Gravity: clean up AccelerationField. */

	 if (level != MaximumGravityRefinementLevel ||
	     MaximumGravityRefinementLevel == MaximumRefinementLevel)
	     Grids[grid1]->GridData->DeleteAccelerationField();

      Grids[grid1]->GridData->DeleteParticleAcceleration();
 
      /* Update current problem time of this subgrid. */
 
      Grids[grid1]->GridData->SetTimeNextTimestep();
 
      /* If using comoving co-ordinates, do the expansion terms now. */
 
      if (ComovingCoordinates)
	Grids[grid1]->GridData->ComovingExpansionTerms();
 
    }  // end loop over grids
 
    /* For each grid: a) interpolate boundaries from the parent grid.
                      b) copy any overlapping zones from siblings. */
 
#ifdef FAST_SIB
    SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			  level, MetaData, Exterior, LevelArray[level]);
#else
    SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
			  Exterior, LevelArray[level]);
#endif

    /* Finalize (accretion, feedback, etc.) star particles */

    StarParticleFinalize(Grids, MetaData, NumberOfGrids, LevelArray,
			 level, AllStars);

    /* If cosmology, then compute grav. potential for output if needed. */

    //dcc cut second potential cut: Duplicate?
 
    if (ComovingCoordinates && SelfGravity && WritePotential) {
      CopyGravPotential = TRUE;
      When = 0.0;
 
#ifdef FAST_SIB
      PrepareDensityField(LevelArray, SiblingList, level, MetaData, When);
#else   // !FAST_SIB
      PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB
 
      CopyGravPotential = FALSE;
 
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
        if (level <= MaximumGravityRefinementLevel) {
 
          /* Compute the potential. */
 
          if (level > 0)
            Grids[grid1]->GridData->SolveForPotential(level);
          Grids[grid1]->GridData->CopyPotentialToBaryonField();
        }
      } //  end loop over grids
    } // if WritePotential
 

    /* For each grid, delete the GravitatingMassFieldParticles. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->DeleteGravitatingMassFieldParticles();


    /* Run the Divergence Cleaing                */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->PoissonSolver(level);
    

 
    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */
 
    MetaData->FirstTimestepAfterRestart = FALSE;

    if (LevelArray[level+1] != NULL) {
      if (EvolveLevel(MetaData, LevelArray, level+1, dtThisLevel, Exterior) == FAIL) {
	fprintf(stderr, "Error in EvolveLevel (%"ISYM").\n", level);
	ENZO_FAIL("");
      }
    }

  

#ifdef USE_JBPERF
    // Update lcaperf "level" attribute

    jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

    OutputFromEvolveLevel(LevelArray,MetaData,level,Exterior);
    CallPython(LevelArray, MetaData, level);

    /* Update SubcycleNumber and the timestep counter for the
       streaming data if this is the bottom of the hierarchy -- Note
       that this not unique based on which level is the highest, it
       just keeps going */

    if (LevelArray[level+1] == NULL) {
      MetaData->SubcycleNumber++;
      MetaData->TimestepCounter++;
    }

    /* ------------------------------------------------------- */
    /* For each grid,
     * (a) project the subgrid's solution into this grid (step #18)
     * (b) correct for the difference between this grid's fluxes and the
     *     subgrid's fluxes. (step #19)
     */
 
#ifdef FLUX_FIX
    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
    CreateSUBlingList(MetaData, Grids,NumberOfGrids, &SUBlingList);
#endif

#ifdef FLUX_FIX
    UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate,SUBlingList,MetaData);
#else
    UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			 SubgridFluxesEstimate);
#endif

#ifdef FLUX_FIX
    DeleteSUBlingList( NumberOfGrids, SUBlingList );
#endif

  /* ------------------------------------------------------- */
  /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
     fluxes for this subgrid .
     (Note: this must be done after CorrectForRefinedFluxes). */

    FinalizeFluxes(Grids,SubgridFluxesEstimate,NumberOfGrids,NumberOfSubgrids);

    /* Recompute radiation field, if requested. */
 
    if (RadiationFieldType >= 10 && RadiationFieldType <= 11 &&
	level <= RadiationFieldLevelRecompute)
      if (RadiationFieldUpdate(LevelArray, level, MetaData) == FAIL) {
	fprintf(stderr, "Error in RecomputeRadiationField.\n");
	ENZO_FAIL("");
      }
 
    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */
 
    if (dtThisLevelSoFar < dtLevelAbove)
      RebuildHierarchy(MetaData, LevelArray, level);

    /* Count up number of grids on this level. */

    int GridMemory, NumberOfCells, CellsTotal, Particles;
    float AxialRatio, GridVolume;
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->CollectGridInformation
        (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
      LevelZoneCycleCount[level] += NumberOfCells;
      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber())
	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }
 
    

    cycle++;
    LevelCycleCount[level]++;
 
  } // end of loop over subcycles
 
  if (debug)
    fprintf(stdout, "EvolveLevel[%"ISYM"]: NumberOfSubCycles = %"ISYM" (%"ISYM" total)\n", level,
           cycle, LevelCycleCount[level]);
 
  /* If possible & desired, report on memory usage. */
 
  ReportMemoryUsage("Memory usage report: Evolve Level");
 
#ifdef USE_JBPERF
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

  if (( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      delete [] SiblingList[grid1].GridList;
    delete [] SiblingList;
  }

  return SUCCESS;
 
}

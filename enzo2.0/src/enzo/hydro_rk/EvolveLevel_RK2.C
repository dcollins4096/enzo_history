/***********************************************************************
/
/  EVOLVE LEVEL USING RUNGE-KUTTA2 METHOD
/
/  written by: Peng Wang & Tom Abel 
/  date:       May, 2007
/  modified1:
/
/  PURPOSE:
/         Evolve Hydro & MHD using 2nd order Runge-Kutta method.
/
************************************************************************/
#include "preincludes.h"

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "Grid.h"
#include "LevelHierarchy.h"
#include "../hydro_rk/tools.h"
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif

/* function prototypes */

int StarParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
			   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			   int ThisLevel, Star *&AllStars,
			   int TotalStarParticleCountPrevious[]);
int StarParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			 int level, Star *&AllStars,
			 int TotalStarParticleCountPrevious[]);
int AdjustRefineRegion(LevelHierarchyEntry *LevelArray[], 
		       TopGridData *MetaData, int EL_level);
int ComputeDednerWaveSpeeds(TopGridData *MetaData,LevelHierarchyEntry *LevelArray[], 
			    int level, FLOAT dt0);
#ifdef TRANSFER
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *AllStars, FLOAT GridTime, int level, int LoopTime = TRUE);
int RadiativeTransferPrepare(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *&AllStars,
			     float dtLevelAbove);
#endif

int CreateSiblingList(HierarchyEntry ** Grids, int NumberOfGrids, 
		      SiblingGridList *SiblingList, int StaticLevelZero, 
		      TopGridData* MetaData, int level);
void DeleteFluxes(fluxes *Fluxes);
int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid, 
			     int NumberOfGrids, int level, float dt);
#ifdef USE_MPI
int CommunicationReduceValues(float *Values, int Number, MPI_Op ReduceOperation);
#endif
int CommunicationAllSumValues(float *Values, int Number);
float CommunicationMinValue(float Value);
float CommunicationMaxValue(float Value);
int CommunicationBarrier();
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CallProblemSpecificRoutines(TopGridData * MetaData, HierarchyEntry *ThisGrid,
				int GridNum, float *norm, float TopGridTimeStep, 
				int level, int LevelCycleCount[]);  //moo

#ifdef FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData, FLOAT When);
int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    SiblingGridList SiblingList[],
			    int level, TopGridData *MetaData,
			    ExternalBoundary *Exterior,
			    LevelHierarchyEntry * Level,
			    int CycleNumber);
#else  // !FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
                        int level, TopGridData *MetaData, FLOAT When);
int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    int level, TopGridData *MetaData, 
			    ExternalBoundary *Exterior,
			    LevelHierarchyEntry * Level,
			    int CycleNumber);
#endif  // end FAST_SIB

int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior);
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




int OutputFromEvolveLevel(LevelHierarchyEntry *LevelArray[],TopGridData *MetaData,
			  int level, ExternalBoundary *Exterior
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
			  );

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

#ifdef FLUX_FIX
#ifdef FAST_SIB 
int CreateSUBlingList(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level,
		      SiblingGridList SiblingList[],
		      LevelHierarchyEntry ***SUBlingList);
#else
int CreateSUBlingList(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level,
		      LevelHierarchyEntry ***SUBlingList);
#endif /* FAST_SIB */
int DeleteSUBlingList(int NumberOfGrids,
		      LevelHierarchyEntry **SUBlingList);
#endif

// int UpdateFromFinerGrids(HierarchyEntry *Grids[], int NumberOfGrids,
// 			 int NumberOfSubgrids[], 
// 			 fluxes **SubgridFluxesEstimate[]);
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);
int FinalizeFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]);		 
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);
int WriteStreamData(LevelHierarchyEntry *LevelArray[], int level,
                    TopGridData *MetaData, int *CycleCount, int open=FALSE);
//int WriteMovieData(char *basename, int filenumber, 
//		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData, 
//		   FLOAT WriteTime);
int WriteTracerParticleData(char *basename, int filenumber, 
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData, 
		   FLOAT WriteTime);
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
                                 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
int FastSiblingLocatorInitializeStaticChainingMesh(ChainingMeshStructure *Mesh, int Rank,
						   int TopGridDims[]); 

double ReturnWallTime();
int CallPython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level);
int SetLevelTimeStep(HierarchyEntry *Grids[], int NumberOfGrids, int level, 
		     float *dtThisLevelSoFar, float *dtThisLevel, 
		     float dtLevelAbove);

void my_exit(int status);
int CallPython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level);

/* Counters for performance and cycle counting. */

static int MovieCycleCount[MAX_DEPTH_OF_HIERARCHY];
static double LevelWallTime[MAX_DEPTH_OF_HIERARCHY];
static double LevelZoneCycleCount[MAX_DEPTH_OF_HIERARCHY];
static double LevelZoneCycleCountPerProc[MAX_DEPTH_OF_HIERARCHY];
 
static float norm = 0.0;            //AK
static float TopGridTimeStep = 0.0; //AK

static int StaticSiblingListInitialized = 0;

#ifdef STATIC_SIBLING_LIST
static SiblingGridList StaticSiblingList[MAX_NUMBER_OF_SUBGRIDS];
static int StaticLevelZero = 1;
#else
static int StaticLevelZero = 0;
#endif

extern int RK2SecondStepBaryonDeposit;

/* EvolveGrid function */

int EvolveLevel_RK2(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		    int level, float dtLevelAbove, ExternalBoundary *Exterior, 
#ifdef TRANSFER
		    ImplicitProblemABC *ImplicitSolver,
#endif
		    FLOAT dt0)
{

  float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid;
  int RefinementFactors[MAX_DIMENSION];
  int cycle = 0, counter = 0, grid1, subgrid, iLevel;
  HierarchyEntry *NextGrid;
  double time1 = ReturnWallTime();

  // Update lcaperf "level" attribute

  Eint32 jb_level = level;
#ifdef USE_JBPERF
  jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

#ifdef FLUX_FIX
  /* Create a SUBling list of the subgrids */
  LevelHierarchyEntry **SUBlingList;
#endif

  FLOAT When;

  
  /* Create an array (Grids) of all the grids. */

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];
  int *TotalStarParticleCountPrevious = new int[NumberOfGrids];

  /* Initialize the chaining mesh used in the FastSiblingLocator. */


  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];
  CreateSiblingList(Grids, NumberOfGrids, SiblingList, StaticLevelZero, 
		    MetaData, level);

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
 
  /* Count the number of colours in the first grid (to define NColor) */

  Grids[0]->GridData->SetNumberOfColours();
  //  fprintf(stdout, "EvolveLevel_RK2: NColor = %"ISYM", NSpecies = %"ISYM"\n", NColor, NSpecies); 

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
    StarParticleInitialize(Grids, MetaData, NumberOfGrids, LevelArray,
			   level, AllStars, TotalStarParticleCountPrevious);

 
#ifdef TRANSFER
    RadiativeTransferPrepare(LevelArray, level, MetaData, AllStars, dtLevelAbove);
#endif /* TRANSFER */

    /* For each grid, compute the number of it's subgrids. */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      NextGrid = Grids[grid1]->NextGridNextLevel;
      counter = 0;
      while (NextGrid != NULL) {
	NextGrid = NextGrid->NextGridThisLevel;
	if (++counter > MAX_NUMBER_OF_SUBGRIDS)
	  ENZO_FAIL("More subgrids than MAX_NUMBER_OF_SUBGRIDS.");
      }
      NumberOfSubgrids[grid1] = counter + 1;
    }


    /* For each grid, create the subgrid list. */

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
      if (MyProcessorNumber == 
	  Grids[grid1]->GridData->ReturnProcessorNumber()) {
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

      } // end: (is this grid on this processor)

    } // end loop over grids (create Subgrid list)

    /* compute wave speed Reference: Matsumoto, PASJ, 2007, 59, 905 */

    if (HydroMethod == MHD_RK) 
      ComputeDednerWaveSpeeds(MetaData, LevelArray, level, dt0);
	
    if (debug && HydroMethod == MHD_RK) 
      fprintf(stderr, "wave speeds: timestep: %"GSYM"  C_h: %"GSYM"  C_p: %"GSYM"\n ", 
	       dt0, C_h, C_p);


    When = 0.5;
    RK2SecondStepBaryonDeposit = 0;
    if (SelfGravity) {
#ifdef FAST_SIB
      PrepareDensityField(LevelArray, SiblingList, level, MetaData, When);
#else   // !FAST_SIB
      PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB
    }

    /* Solve the radiative transfer */

#ifdef TRANSFER
    FLOAT GridTime = Grids[0]->GridData->ReturnTime();
    EvolvePhotons(MetaData, LevelArray, AllStars, GridTime, level);
#endif /* TRANSFER */

    /* Compute particle-particle acceleration */

    /*    if (NBodyDirectSummation == TRUE) 
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->ComputeParticleParticleAcceleration(level); */

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. Predictor Step*/

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      
      CallProblemSpecificRoutines(MetaData, Grids[grid1], grid1, &norm, 
				  TopGridTimeStep, level, LevelCycleCount);

      /* Gravity: compute acceleration field for grid and particles. */
      if (SelfGravity) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {
	  if (level > 0) 
	    Grids[grid1]->GridData->SolveForPotential(level) ;
	  Grids[grid1]->GridData->ComputeAccelerations(level) ;
	}
	// otherwise, use interpolated potential from coarser grid, which is
	// done in PrepareDensity.
      } // end: if (SelfGravity)

      Grids[grid1]->GridData->ComputeAccelerationFieldExternal() ;

#ifdef TRANSFER
      /* Radiation Pressure: add to acceleration field */
      Grids[grid1]->GridData->AddRadiationPressureAcceleration();
#endif /* TRANSFER */

      Grids[grid1]->GridData->CopyBaryonFieldToOldBaryonField();  

      if (UseHydro) 
	if (HydroMethod == HD_RK)
	  Grids[grid1]->GridData->RungeKutta2_1stStep
	    (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);
	else if (HydroMethod == MHD_RK) 
	  Grids[grid1]->GridData->MHDRK2_1stStep
	    (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);
	
      /* Do this here so that we can get the correct time interpolated boundary condition */
      Grids[grid1]->GridData->SetTimeNextTimestep();
      
    }  // end loop over grids

#ifdef FAST_SIB
      SetBoundaryConditions(Grids, NumberOfGrids, SiblingList, level, MetaData, Exterior, LevelArray[level]);
#else
      SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, Exterior, LevelArray[level]);
#endif

    // Recompute potential and accelerations with time centered baryon Field
    // this also does the particles again at the moment so could be made more efficient.

    RK2SecondStepBaryonDeposit = 0; // set this to (0/1) to (not use/use) this extra step
    //    printf("SECOND STEP\n");
    if (RK2SecondStepBaryonDeposit && SelfGravity && UseHydro) {  
      When = 0.5;
#ifdef FAST_SIB
      PrepareDensityField(LevelArray, SiblingList, level, MetaData, When);
#else  
      PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB
    }

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Gravity: compute acceleration field for grid and particles. */
      if (RK2SecondStepBaryonDeposit) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {
	  if (level > 0) 
	    Grids[grid1]->GridData->SolveForPotential(level) ;
	  Grids[grid1]->GridData->ComputeAccelerations(level) ;
	}
	Grids[grid1]->GridData->ComputeAccelerationFieldExternal() ;
      } // end: if (SelfGravity)

      if (UseHydro) {
	if (HydroMethod == HD_RK)
	  Grids[grid1]->GridData->RungeKutta2_2ndStep
	    (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);

	else if (HydroMethod == MHD_RK) {
	  
	  Grids[grid1]->GridData->MHDRK2_2ndStep
	    (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);
	  
	  if (UseAmbipolarDiffusion) 
	    Grids[grid1]->GridData->AddAmbipolarDiffusion();
	  
	  if (UseResistivity) 
	    Grids[grid1]->GridData->AddResistivity();
	
	  time1 = ReturnWallTime();
	  
	 
	
	} // ENDIF MHD_RK
      } // ENDIF UseHydro

      /* Add viscosity */

      if (UseViscosity) {
	Grids[grid1]->GridData->AddViscosity();

      }


      /* Solve the cooling and species rate equations. */
 
      Grids[grid1]->GridData->MultiSpeciesHandler();

      /* Update particle positions (if present). */

      UpdateParticlePositions(Grids[grid1]->GridData);

      /* Include 'star' particle creation and feedback. */

      Grids[grid1]->GridData->StarParticleHandler
	(Grids[grid1]->NextGridNextLevel, level
#ifdef EMISSIVITY 
			  , dtLevelAbove
#endif
        );
 
      /* Gravity: clean up AccelerationField. */

	 if (level != MaximumGravityRefinementLevel ||
	     MaximumGravityRefinementLevel == MaximumRefinementLevel)
	     Grids[grid1]->GridData->DeleteAccelerationField();



      Grids[grid1]->GridData->DeleteParticleAcceleration();
 
      if (UseFloor) 
	Grids[grid1]->GridData->SetFloor();

      /* If using comoving co-ordinates, do the expansion terms now. */

      if (ComovingCoordinates)
	Grids[grid1]->GridData->ComovingExpansionTerms();

    }  // end loop over grids

    if (UseDivergenceCleaning != 0){

#ifdef FAST_SIB
      SetBoundaryConditions(Grids, NumberOfGrids, SiblingList, level, MetaData, Exterior, LevelArray[level]);
#else
      SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, Exterior, LevelArray[level]);
#endif
      
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

	Grids[grid1]->GridData->PoissonSolver(level);
      }
      
    }
    
#ifdef FAST_SIB
    SetBoundaryConditions(Grids, NumberOfGrids, SiblingList, level, MetaData, Exterior, LevelArray[level]);
#else
    SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, Exterior, LevelArray[level]);
#endif

    /* Finalize (accretion, feedback, etc.) star particles */
 
    StarParticleFinalize(Grids, MetaData, NumberOfGrids, LevelArray,
			 level, AllStars, TotalStarParticleCountPrevious);


    OutputFromEvolveLevel(LevelArray,MetaData,level,Exterior
#ifdef TRANSFER
			  , ImplicitSolver
#endif
			  );
    CallPython(LevelArray, MetaData, level);

    /* For each grid, delete the GravitatingMassFieldParticles. */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->DeleteGravitatingMassFieldParticles();


    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */

    MetaData->FirstTimestepAfterRestart = FALSE;

    //    LevelWallTime[level] += ReturnWallTime() - time1;
    if (LevelArray[level+1] != NULL) {
      if (EvolveLevel_RK2(MetaData, LevelArray, level+1, dtThisLevel, Exterior, 
#ifdef TRANSFER
			  ImplicitSolver, 
#endif
			  dt0) 
	  == FAIL) {
	fprintf(stderr, "Error in EvolveLevel_RK2 (%"ISYM").\n", level);
	ENZO_FAIL("");
      }
    }
    time1 = ReturnWallTime();

#ifdef USE_JBPERF
    // Update lcaperf "level" attribute
    jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

    OutputFromEvolveLevel(LevelArray, MetaData, level, Exterior
#ifdef TRANSFER
			  , ImplicitSolver
#endif
			  );
    CallPython(LevelArray, MetaData, level);

    /* Update SubcycleNumber and the timestep counter for the
       streaming data if this is the bottom of the hierarchy -- Note
       that this not unique based on which level is the highest, it
       just keeps going */

    if (LevelArray[level+1] == NULL) {
      MetaData->SubcycleNumber++;
      MetaData->MovieTimestepCounter++;
    }

    /* ------------------------------------------------------- */
    /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */


#ifdef FLUX_FIX
    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
#ifdef FAST_SIB
    CreateSUBlingList(MetaData, LevelArray, level, SiblingList,
		      &SUBlingList);
#else
    CreateSUBlingList(MetaData, LevelArray, level, &SUBlingList);
#endif /* FAST_SIB */
#endif

#ifdef FLUX_FIX
    UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			 SubgridFluxesEstimate,
			 SUBlingList,
			 MetaData);
#else
    UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids, SubgridFluxesEstimate);
#endif
    
#ifdef FLUX_FIX        /* Clean up SUBlings */
    DeleteSUBlingList( NumberOfGrids, SUBlingList );
#endif

     
    /* ------------------------------------------------------- */
    /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
       fluxes for this subgrid .
       (Note: this must be done after CorrectForRefinedFluxes). */

    FinalizeFluxes(Grids,SubgridFluxesEstimate,NumberOfGrids,NumberOfSubgrids);
    
    /* Recompute radiation field, if requested. */

    RadiationFieldUpdate(LevelArray, level, MetaData);

    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */

    //    LevelWallTime[level] += ReturnWallTime() - time1;
    if (dtThisLevelSoFar < dtLevelAbove) 
      RebuildHierarchy(MetaData, LevelArray, level);

    time1 = ReturnWallTime();

    /* Count up number of grids on this level. */

    int GridMemory, NumberOfCells, CellsTotal, Particles;
    float AxialRatio, GridVolume;
#ifdef UNUSED
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->CollectGridInformation
        (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
      LevelZoneCycleCount[level] += NumberOfCells;
      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber())
      	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }
#endif

    cycle++;
    LevelCycleCount[level]++;

  } // while (dtThisLevelSoFar < dtLevelAbove)

  if (debug)
    printf("EvolveLevelRK2[%"ISYM"]: NumberOfSubCycles = %"ISYM" (%"ISYM" total)\n", level, 
           cycle, LevelCycleCount[level]);

  /* If possible & desired, report on memory usage. */

  ReportMemoryUsage("Memory usage report: Evolve Level");

  /* Clean up. */

  delete [] NumberOfSubgrids;
  delete [] Grids;
  delete [] SubgridFluxesEstimate;

  /* Clean up the sibling list. */

  if (( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {
    for (int grid1 = 0; grid1 < NumberOfGrids; grid1++)
      delete [] SiblingList[grid1].GridList;
    delete [] SiblingList;
  }

  return SUCCESS;

}

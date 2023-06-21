/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
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
//take this out soon.
#ifdef USE_MPI
#include <mpi.h>
#ifdef USE_MPE
#include <mpe.h>
#endif /* USE_MPE */
#endif /* USE_MPI */

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
#include "pout.h"
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
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
                        SiblingGridList SiblingList[],
#endif
                        int level, TopGridData *MetaData);

int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
                          SiblingGridList SiblingList[],
#endif
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);

int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
			    SiblingGridList SiblingList[],
#endif
			  int level, TopGridData *MetaData, 
			    ExternalBoundary *Exterior, LevelHierarchyEntry * Level,
			    int CycleNumber);


int UpdateFromFinerGrids(HierarchyEntry *Grids[], LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData,
			 int NumberOfGrids,
			 int NumberOfSubgrids[], 
			 fluxes **SubgridFluxesEstimate[]
#ifdef JB_OPT_FLUXES_FIX
			 , LevelHierarchyEntry *SUBlingList[]
#endif
			 , int Cycle //dcc this is only for debugging.
);

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

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
                                 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
#endif

#ifdef JB_OPT_FLUXES_FIX
int CreateSUBlingList(TopGridData *MetaData,
                      HierarchyEntry *Grids[],
                      int NumberOfGrids,
                      LevelHierarchyEntry ***SUBlingList);
int DeleteSUBlingList( int NumberOfGrids,
                       LevelHierarchyEntry **SUBlingList);
#ifdef DC_OPT_SIBSUB_II
int CreateSUBlingListFast(TopGridData *MetaData,
                      HierarchyEntry *Grids[],
                      int NumberOfGrids,
                      LevelHierarchyEntry ***SUBlingList);
#endif


#endif

int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float * norm, float *bulkMomentum,
				      float * pTopGridTimeStep);

int ComputeGlobalStats(LevelHierarchyEntry *LevelArray[],
		       int level, TopGridData *MetaData, float *);

int CommunicationCheckOutstanding(); 
#ifndef DCC_EXTERNAL_CYCLE_COUNT
static int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
#endif //!DCC_EXTERNAL_CYCLE_COUNT
static float norm = 0.0;            //AK
static float bulkMomentum[3] = {0.0, 0.0, 0.0}; //AK
static float TopGridTimeStep = 0.0; //AK


//dcc.  This was an attempt to understand why my job was crashing in DeleteFluxes.
//      It shed no light on the subject.
#define DC_ONLYALLOCATEMYFLUXES

/* EvolveGrid function */

int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior)
{

    //<dbg>
    if (WriteInThisF(500)){
      fprintf(stderr,"LBbug: EL dump\n");
    
      //MetaData->OutputFirstTimeAtLevel = level+1;
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

    //</dbg>


  int oot = (MyProcessorNumber == ROOT_PROCESSOR ) ? TRUE : FALSE;
  if( oot ) {
    fprintf(stderr, "\n========== EL: level %d ==========\n", level);
  }
  wall_time("Start EvolveLevel");

  /* Declarations */

  float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid;
  int RefinementFactors[MAX_DIMENSION];
  int cycle = 0, counter = 0, grid, subgrid;
  HierarchyEntry *NextGrid;

  JBPERF_INIT;

  //  JBPERF_START_LEVEL("EL00"); // EvolveLevel

  /* Create an array (Grids) of all the grids. */

  JBPERF_START_LEVEL("EL01"); // GenerateGridArray ()

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];

#ifdef JB_OPT_FLUXES_FIX
  /* Create a SUBling list of the subgrids */

  LevelHierarchyEntry **SUBlingList;
#endif

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
  /* Initialize the chaining mesh used in the FastSiblingLocator. */

  ChainingMeshStructure ChainingMesh;
  FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
                               MetaData->TopGridDims);
  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];

  /* Add all the grids to the chaining mesh. */

  
  for (grid = 0; grid < NumberOfGrids; grid++)
    Grids[grid]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh
						     );

  /* For each grid, get a list of possible siblings from the chaining mesh. */

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->FastSiblingLocatorFindSiblings(
							      &ChainingMesh, &SiblingList[grid],
							      MetaData->LeftFaceBoundaryCondition,
							      MetaData->RightFaceBoundaryCondition) == FAIL) {
      fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
      return FAIL;
    }

  /* Clean up the chaining mesh. */

  FastSiblingLocatorFinalize(&ChainingMesh);

#endif


  JBPERF_STOP_LEVEL("EL01"); // Create grid array   

  /* ================================================================== */
  /* For each grid: a) interpolate boundaries from its parent.
                    b) copy any overlapping zones.  */

  JBPERF_START_LEVEL("EL02"); // SetBoundaryConditions ()

  JBMEM_MESSAGE(MyProcessorNumber,"jb: EnterEvolveLevel");

  //wall_time("Start SBC1");
  if (SetBoundaryConditions(Grids, NumberOfGrids,
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
                            SiblingList,
#endif
			    level, MetaData, 
			    Exterior, LevelArray[level]) == FAIL)
    return FAIL;

  JBMEM_MESSAGE(MyProcessorNumber,"jb: SetBoundaryConditions1");

  CommunicationCheckOutstanding();
  JBPERF_STOP_LEVEL("EL02"); // SetBoundaryConditions ()
  //wall_time("End SBC1");
  //=================================================================
  /* Clear the boundary fluxes for all Grids (this will be accumulated over
     the subcycles below (i.e. during one current grid step) and used to by the
     current grid to correct the zones surrounding this subgrid (step #18). */

  JBPERF_START_LEVEL("EL03"); // ClearBoundaryFluxes ()

  for (grid = 0; grid < NumberOfGrids; grid++){
    
    Grids[grid]->GridData->ClearBoundaryFluxes();

    //only call on subgrids: Wastefull if called on root grid.
    if(MHD_Used == TRUE && level > 0 ){
      Grids[grid]->GridData->ClearAvgElectricField();
    }
  }

  JBPERF_STOP_LEVEL("EL03"); // ClearBoundaryFluxes ()

  

  /* ================================================================== */
  /* Loop over grid timesteps until the elapsed time equals the timestep
     from the level above (or loop once for the top level). */

  while (dtThisLevelSoFar < dtLevelAbove) {

    //<dbg> For fixed time steps.
    dccCounter9= LevelCycleCount[level];
    //</dbg>
    
    JBMEM_MESSAGE(MyProcessorNumber,"jb: StartSubcycleLoop");
    //=================================================================
    
    /* Determine the timestep for this iteration of the loop. */

    JBPERF_START_LEVEL("EL04"); // SetTimeStep ()



    //MPI_Barrier(MPI_COMM_WORLD);

    if (level == 0) {

      /* For root level, use dtLevelAbove. */

      dtThisLevel      = dtLevelAbove;
      dtThisLevelSoFar = dtLevelAbove;

    } else {

      /* Compute the mininum timestep for all grids. */
      ////wall_time("Start dt");
      ////wall_time("Start dt");
      dtThisLevel = huge_number;
      for (grid = 0; grid < NumberOfGrids; grid++) {
	dtGrid      = Grids[grid]->GridData->ComputeTimeStep(level);
	dtThisLevel = min(dtThisLevel, dtGrid);
      }
      ////wall_time("End dt");
      ////wall_time("End dt");
      //////wall_time("Start dtComm");
      //////wall_time("Start dtComm");
      //////wall_time("Start dtComm");
      //////wall_time("Start dtComm");
      dtThisLevel = CommunicationMinValue(dtThisLevel);
      //////wall_time("End dtComm");
      //////wall_time("End dtComm");
      //////wall_time("End dtComm");
      //////wall_time("End dtComm");
      /* Advance dtThisLevelSoFar (don't go over dtLevelAbove). */
      
      if (dtThisLevelSoFar+dtThisLevel*1.05 >= dtLevelAbove) {
	dtThisLevel      = dtLevelAbove - dtThisLevelSoFar;
	dtThisLevelSoFar = dtLevelAbove;
      }
      else
	dtThisLevelSoFar += dtThisLevel;

    }
    if (debug) printf("Level[%d]: dt = %g(This So Far: %g Above: %g Frac %g)\n", level, dtThisLevel,
		      dtThisLevelSoFar, dtLevelAbove,dtThisLevelSoFar/dtLevelAbove );

    JBMEM_MESSAGE(MyProcessorNumber,"jb: AfterComputeTimestep");


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

    JBMEM_MESSAGE(MyProcessorNumber,"jb: BeforeCreateSubgridList");
    for (grid = 0; grid < NumberOfGrids; grid++) {


      /* Allocate the subgrid fluxes for this grid. */

      SubgridFluxesEstimate[grid] = new fluxes *[NumberOfSubgrids[grid]];

#ifdef DC_ONLYALLOCATEMYFLUXES
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid]; subgrid++)
        SubgridFluxesEstimate[grid][subgrid] = NULL;
#endif


      /* Collect the flux data and store it in the newly minted fluxes.
	 Or rather that's what we should do.  Instead, we create fluxes one
	 by one in this awkward array of pointers to pointers.  This should be
	 changed so that all the routines take arrays of flux rather than
	 arrays of pointers to flux.  Dumb. */

      counter = 0;
#ifdef DC_ONLYALLOCATEMYFLUXES
      if (MyProcessorNumber ==
          Grids[grid]->GridData->ReturnProcessorNumber()) {
#endif

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

#ifdef DC_ONLYALLOCATEMYFLUXES
      }
#endif

    } // end loop over grids (create Subgrid list)

    JBMEM_MESSAGE(MyProcessorNumber,"jb: AfterCreateSubgridList");
    JBPERF_STOP_LEVEL("EL06"); // Create subgrid list

    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */

    //MPI_Barrier(MPI_COMM_WORLD);
    JBPERF_START_LEVEL("EL07"); // PrepareDensityField ()
    wall_time("Start PrepareDensity");
    ////wall_time("Start PrepareDensity");
    if (SelfGravity)
      if (PrepareDensityField(LevelArray,
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
                              SiblingList,
#endif
                              level, MetaData) == FAIL) {
	fprintf(stderr, "Error in PrepareDensityField.\n");
	return FAIL;
      }
  
    wall_time("End PrepareDensity");
    ////wall_time("End PrepareDensity");
    if (RandomForcing && MetaData->CycleNumber > 0 && level == 0)
      if ( ComputeRandomForcingNormalization(LevelArray, 0, MetaData,
                                             &norm, bulkMomentum,&TopGridTimeStep)
           == FAIL ) {
        fprintf(stderr, "Error in ComputeRandomForcingNormalization.\n");
        return FAIL;
      }

    //some global statistics

    if( ComputeGlobalStats(LevelArray,level,MetaData, bulkMomentum) == FAIL )
      {fprintf(stderr,"Error in ComputeGlobalStatistics\n"); return FAIL;}


    JBPERF_STOP_LEVEL("EL07"); // PrepareDensityField ()

    //dcf for all grids on this level
    for (grid = 0; grid < NumberOfGrids; grid++) {

      //=================================================================
      /* Call analysis routines. */

      //MPI_Barrier(MPI_COMM_WORLD);

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
      
      
      if(WriteInThisF(280)==TRUE ){
	for(grid=0;grid<NumberOfGrids; grid++){
	  
	  if( Grids[grid]->GridData->CheckForNans("before write data2") == FAIL ) {
	    fprintf(stderr,"CFN 1 \n");
	    return FAIL;
	  }
	  char basename[30];
	  sprintf(basename,"data28%d%d%d.grid", LevelCycleCount[level], level, grid);
	  FILE * dccptr = fopen(basename,"a");
	  //Grids[grid]->GridData->WriteGrid(dccptr,basename,0);
	  //Grids[grid]->GridData->WriteGrid(dccptr,basename,MyProcessorNumber);
	  Grids[grid]->GridData->WriteGrid(dccptr,basename,1);
	  fclose(dccptr);
	  
	  if( Grids[grid]->GridData->CheckForNans("after write data2") == FAIL ) {
	    fprintf(stderr,"CFN 2\n");
	    return FAIL;
	  }
	}
      }

#ifdef ATHENA      
      //<dcc>
      //Acceleration Field needs to get zeroed.  This is mostly for point sources.
      if (SelfGravity || UniformGravity || PointSourceGravity) {
	if( Grids[grid]->GridData->ZeroAcceleration() == FAIL ){
	  fprintf(stderr,"EvolveLevel: failure in ZeroAccelration\n");
	  return FAIL;
	}
      }
      //</dcc>
#endif //ATHENA

      
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
      
      if(WriteInThisF(290)==TRUE ){
	for(grid=0;grid<NumberOfGrids; grid++){
	  
	  if( Grids[grid]->GridData->CheckForNans("before write data2") == FAIL ) {
	    fprintf(stderr,"CFN 1 \n");
	    return FAIL;
	  }
	  char basename[30];
	  sprintf(basename,"data29%d%d%d.grid", LevelCycleCount[level], level, grid);
	  FILE * dccptr = fopen(basename,"a");
	  //Grids[grid]->GridData->WriteGrid(dccptr,basename,0);
	  //Grids[grid]->GridData->WriteGrid(dccptr,basename,MyProcessorNumber);
	  Grids[grid]->GridData->WriteGrid(dccptr,basename,1);
	  fclose(dccptr);
	  
	  if( Grids[grid]->GridData->CheckForNans("after write data2") == FAIL ) {
	    fprintf(stderr,"CFN 2\n");
	    return FAIL;
	  }
	}
      }
      
      
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

      if(level==0) 
	Grids[grid]->GridData->TotalMass("ComputeAcceleration");
#ifndef ACCELOFF

    }//grid loop
    //dcf endfor
    //This ensures that all subgrids agree in the boundary.
    //Not a big deal for hydro, but essential for DivB = 0 in MHD runs.
    //Only called on level > 0 because the root grid is dealt with differently than SG's.
    //10/3/06, dcc: With the new solver, the root grid acceleration was misbehaving.
    //              Necessitated a call to the root grid, too.
    if(0==0){
#ifdef ATHENA
      if (SelfGravity || UniformGravity || PointSourceGravity)
#else
      if ( (SelfGravity || UniformGravity || PointSourceGravity)  && level > 0) 
#endif
      if( SetAccelerationBoundary(Grids, NumberOfGrids, 
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
				  SiblingList,
#endif
				  level, MetaData, 
				  Exterior, LevelArray[level], LevelCycleCount[level]) == FAIL ){
	fprintf(stderr,"Error with AccelerationBoundary.\n");
	return FAIL;
      }
    }else{fprintf(stderr,"kludge: no SAB\n");}


#ifndef OLD_CENTER

    //dcc 01/24/06.  
    // The new unsplit MagneticField Centering requires a boundary condition set.
    // In truth, only CopyZonesFromGrid needs to be called, and only on the CenteredB field.
    // However, creating such a routine will invlove slightly altering CopyZonesFromGrid and 
    // CommunicationSendRegion to deal only with CenteredB, which is a good 2 days of work.
    // Also note that other centering methods probably dont' need an entire boundary check anyhow.

    if( MHD_CenteringMethod != MHD_Volumetric ) {
      if( MHD_Used == TRUE ){
	for (grid = 0; grid < NumberOfGrids; grid++) {
	  //If you do move CMF into SBC, make a comment that CMF is the important routine.
	  if( Grids[grid]->GridData->CenterMagneticField() == FAIL ){
	    fprintf( stderr, "error in Center MagneticField\n");
	    return FAIL;
	  }
	}//Grids
      }else{
	if (SetBoundaryConditions(Grids, NumberOfGrids,
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
				  SiblingList,
#endif
				  level, MetaData, 
				  Exterior, LevelArray[level]) == FAIL)
	  return FAIL;
      }//centering method
      }//mhd_used
#endif //OLD_CENTER
    
    
    //
    // re-start grid loop.
    //
    
    //dcf for all grids on this level
    for (grid = 0; grid < NumberOfGrids; grid++) {
#endif // ACCELOFF      
    
      //=================================================================
      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
    
      JBPERF_START_LEVEL("EL11"); // CopyBaryonFieldToOldBaryonField ()
      JBMEM_MESSAGE(MyProcessorNumber,"jb: BeforeCopyBaryonField");

      if (Grids[grid]->GridData->CopyBaryonFieldToOldBaryonField() == FAIL) {
	fprintf(stderr, "Error in grid->CopyBaryonFieldToOldBaryonField.\n");
	return FAIL;
      }
      JBMEM_MESSAGE(MyProcessorNumber,"jb: AfterCopyBaryonField");

      //Random forcing terms.
      //Athena is 2nd order in time, so forcing must be done in 


      if (RandomForcing && MetaData->CycleNumber > 0
#ifdef ATHENA
	  && HydroMethod != Athena 
#endif //ATHENA	  
	  ){ //AK
	if(Grids[grid]->GridData->AddRandomForcing(&norm, bulkMomentum,
						   TopGridTimeStep) == FAIL)
	  fprintf(stderr, "Error in AddRandomForcing.\n");
      } //random
      
      JBPERF_STOP_LEVEL("EL11"); // CopyBaryonFieldToOldBaryonField ()
      
      //=================================================================
      /* Call hydro solver and save fluxes around subgrids. */

      //MPI_Barrier(MPI_COMM_WORLD);
      
      JBPERF_START_LEVEL("EL12"); // SolveHydroEquations Or MHD ()
      JBMEM_MESSAGE(MyProcessorNumber,"jb: BeforeSMHD");

      if( MHD_Used == TRUE )
	{

#ifdef ATHENA
	  if( HydroMethod == Athena ){
	    if( Grids[grid]->GridData->NewSMHD(LevelCycleCount[level], 
	       NumberOfSubgrids[grid], SubgridFluxesEstimate[grid], Exterior, level, grid,
					       &norm, TopGridTimeStep) == FAIL ){
	      fprintf(stderr, "Error in grid->NewSMHD.\n");
	      return FAIL;
	    }
	  }else{
	  if( Grids[grid]->GridData->SolveMHDEquations(LevelCycleCount[level], 
       NumberOfSubgrids[grid], SubgridFluxesEstimate[grid], Exterior, level, grid) == FAIL ){
	    fprintf(stderr, "Error in grid->SolveMHDEquations.\n");
	    return FAIL;
	  }
	  }//hydro method switch.
#else  //ATHENA
	  if( Grids[grid]->GridData->SolveMHDEquations(LevelCycleCount[level], 
       NumberOfSubgrids[grid], SubgridFluxesEstimate[grid], Exterior, level, grid) == FAIL ){
	    fprintf(stderr, "Error in grid->SolveMHDEquations.\n");
	    return FAIL;
	  }

#endif //ATHENA
	}else{
	  if (Grids[grid]->GridData->SolveHydroEquations(LevelCycleCount[level],
	      NumberOfSubgrids[grid], SubgridFluxesEstimate[grid], level, grid) == FAIL) {
	    fprintf(stderr, "Error in grid->SolveHydroEquations.\n");
	    return FAIL;
	  }
	}

      if(level==0)
	Grids[grid]->GridData->TotalMass("after solve hydro");

      JBMEM_MESSAGE(MyProcessorNumber,"jb: AfterSMHD");

      JBPERF_STOP_LEVEL("EL12"); // SolveHydroEquations Or MHDEquations ()

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
	
	if (Grids[grid]->GridData->StarParticleHandler(level) == FAIL) {
	  fprintf(stderr, "Error in grid->StarParticleWrapper");
	  return FAIL;
	}
      }
         
      JBPERF_STOP_LEVEL("EL16"); // Include star particle creation and feedback

      //=================================================================
      /* Gravity: clean up AccelerationField. */

      JBPERF_START_LEVEL("EL17"); // clean up gravity acceleration field

      //<dcc> This needs to be removed for the gravity hack.
      /*
      if (SelfGravity || UniformGravity || PointSourceGravity) {
	if (level != MaximumGravityRefinementLevel ||
	    MaximumGravityRefinementLevel == MaximumRefinementLevel)
	  Grids[grid]->GridData->DeleteAccelerationField();
	Grids[grid]->GridData->DeleteParticleAcceleration();
      }
      */
      //</dcc> end removal.
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

#ifdef BIERMANN_
      // Computer the Biermann Bettery term here
      if(BiermannBattery){
	Grids[grid]->GridData->ComputeBiermannTerms();
	Grids[grid]->GridData->DeleteBiermannTerms(); 
      }
#endif //BIERMANN

    } 
    //<dbg>
    if( SuggestFailure == TRUE ){
      fprintf(stderr,"Failure has been suggested by some routine in the grid loop.\n");
      return FAIL;
    }
    //</dbg>
    ///////////////
    /////////////// end grid loop
    ///////////////
    ///////////////
    //dcf endfor

    if(0==0){
      //=================================================================
      /* For each grid: a) interpolate boundaries from the parent grid.
	 b) copy any overlapping zones from siblings. */
      
      JBPERF_START_LEVEL("EL20"); // SetBoundaryConditions ()
      //wall_time("Start SBC2");
      if (SetBoundaryConditions(Grids, NumberOfGrids, 
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
				SiblingList,
#endif
				level, MetaData, 

				Exterior, LevelArray[level]) == FAIL)
	return FAIL;

      JBPERF_STOP_LEVEL("EL20"); // SetBoundaryConditions ()
      //wall_time("End SBC2");
    }
    



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

    if (LevelArray[level+1] == NULL){
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
    }
    JBPERF_STOP_LEVEL("EL22"); // WriteMovieData ()
     
    //=================================================================
    /* Check for new level output (only if this is bottom of hierarchy). */

    JBPERF_START_LEVEL("EL23"); // WriteAllData ()


    //<dbg>
    if (WriteInThisF(501)){
      fprintf(stderr,"LBbug: EL dump\n");
    
      //MetaData->OutputFirstTimeAtLevel = level+1;
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

    //</dbg>
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

    //    JBPERF_START_LEVEL("EL25"); // EvolveLevel recursion


    if (LevelArray[level+1] != NULL){
      JBMEM_MESSAGE(MyProcessorNumber,"jb: RecursiveToEL");

      if (EvolveLevel(MetaData, LevelArray, level+1, dtThisLevel, Exterior) 
	  == FAIL) {
	fprintf(stderr, "Error in EvolveLevel (%d).\n", level);
	return FAIL;
      }
      JBMEM_MESSAGE(MyProcessorNumber,"jb: PopBackFromRecursion");

    }
 
  
#ifdef JB_OPT_FLUXES_FIX

    JBMEM_MESSAGE(MyProcessorNumber,"jb: Before CreateSUBling");

    //Pout("CSshd: Call");
    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];

    for(int list=0;list<NumberOfGrids;list++)
      SUBlingList[list] = NULL;

    if (FluxCorrection) {
      /* Fill in the SUBling list */

#ifdef DC_OPT_SIBSUB_II

      if (CreateSUBlingListFast(MetaData, Grids,
				NumberOfGrids, &SUBlingList) == FAIL){
        fprintf(stderr, "Error in CreateSUBlingList.\n");
        return FAIL;
      }

#else

      if (CreateSUBlingList(MetaData, Grids,
			    NumberOfGrids, &SUBlingList) == FAIL){
        fprintf(stderr, "Error in CreateSUBlingList.\n");
        return FAIL;}
#endif

      /*
	if(level == 0 )
	for(grid=0;grid<NumberOfGrids;grid++)
	fprintf(stderr," sublinglist %d %p\n", grid, SUBlingList[grid]);
      */
	
    }//FluxCorrection
    JBMEM_MESSAGE(MyProcessorNumber,"jb: After CreateSUBling");    
#endif /* JB_OPT_FLUXES_FIX */


    //<dcc>
    if(WriteInThisF(19)==TRUE ){
      for(grid=0;grid<NumberOfGrids; grid++){

	char basename[30];
	sprintf(basename,"data19%d%d%d.grid", LevelCycleCount[level], level, grid);
	FILE * dccptr = fopen(basename,"a");
	//Grids[grid]->GridData->WriteGrid(dccptr,basename,0);
	//Grids[grid]->GridData->WriteGrid(dccptr,basename,MyProcessorNumber);
	Grids[grid]->GridData->WriteGrid(dccptr,basename,1);
	fclose(dccptr);

      }      
    }
    //</dcc>

    /* ------------------------------------------------------- */
    /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */
 
    JBPERF_START_LEVEL("EL26"); // UpdateFromFinerGrids ()

	
      if(level==0)
	Grids[0]->GridData->TotalMass("Before UFG");
    ////wall_time("Start UFG");    
    ////wall_time("Start UFG");    
    if(UpdateFromFinerGrids(Grids, LevelArray, level,MetaData, NumberOfGrids, NumberOfSubgrids, 
			    SubgridFluxesEstimate
#ifdef JB_OPT_FLUXES_FIX
			     , SUBlingList
#endif
			    , LevelCycleCount[level] //dcc this is only for debugging.
) == FAIL)
      return FAIL;

      if(level==0)
	Grids[0]->GridData->TotalMass("after UFG");

    JBMEM_MESSAGE(MyProcessorNumber,"jb: Before DeleteSUBling");
#ifdef JB_OPT_FLUXES_FIX
    if( FluxCorrection==1){
      /* Clean up SUBlings */
      if (DeleteSUBlingList( NumberOfGrids, SUBlingList ) == FAIL){
        fprintf(stderr, "Error in DeleteSUBlingList.\n");
        return FAIL;
      }
    }
#endif /* JB_OPT_FLUXES_FIX */
    ////wall_time("End UFG"); 
    ////wall_time("End UFG"); 
    JBMEM_MESSAGE(MyProcessorNumber,"jb: After DeleteSUBling");
    
    if(WriteInThisF(20)==TRUE ){
      for(grid=0;grid<NumberOfGrids; grid++){

	if( Grids[grid]->GridData->CheckForNans("before write data2") == FAIL ) {
	  fprintf(stderr,"CFN 1 \n");
	  return FAIL;
	}
	char basename[30];
	sprintf(basename,"data20%d%d.grid", LevelCycleCount[level], level);
	FILE * dccptr = fopen(basename,"a");
	//Grids[grid]->GridData->WriteGrid(dccptr,basename,0);
	//Grids[grid]->GridData->WriteGrid(dccptr,basename,MyProcessorNumber);
	Grids[grid]->GridData->WriteGrid(dccptr,basename,grid+1);
	fclose(dccptr);

	if( Grids[grid]->GridData->CheckForNans("after write data2") == FAIL ) {
	  fprintf(stderr,"CFN 2\n");
	  return FAIL;
	}
      }
    }
    
    JBMEM_MESSAGE(MyProcessorNumber,"jb: AterUFG");
    
    JBPERF_STOP_LEVEL("EL26"); // UpdateFromFinerGrids ()

    //Now that the electric field is updated to its fullest, create MagneticField
    ////wall_time("Start UMF");
    ////wall_time("Start UMF");
    if(MHD_Used == TRUE && MHD_ProjectE == TRUE){
      
      JBPERF_START_LEVEL("EL30"); // UpdateMagneticField
      
      for(grid=0;grid<NumberOfGrids; grid++){
	
	//dcc moved a write grid from here to before the ProjectE conditional
	
	if( Grids[grid]->GridData->MHD_UpdateMagneticField(level, LevelArray[level+1]) == FAIL )
	  return FAIL;
	
	if(WriteInThisF(21)==TRUE){
	  
	  char basename[30];
	  sprintf(basename,"data21%d%d.grid", LevelCycleCount[level],level);
	  FILE * dccptr = fopen(basename,"a");
	  Grids[grid]->GridData->WriteGrid(dccptr,basename,grid+1);
	  fclose(dccptr);
	}
	
      }
      JBPERF_STOP_LEVEL("EL30"); // UpdateMagneticField
    }//MHD True
    ////wall_time("End UMF");
    ////wall_time("End UMF");
    //=================================================================
    /* For each grid: a) interpolate boundaries from the parent grid.
       b) copy any overlapping zones from siblings. */
    
    
    JBPERF_START_LEVEL("EL20"); // SetBoundaryConditions ()
    //I ultimately want this here, I think.
    JBMEM_MESSAGE(MyProcessorNumber,"jb: AfterUMF");
    
    
    //fprintf(stderr," === EL: Second SBC === \n");
    //wall_time("Start SBC3");
       //<mhd change>
    if (SetBoundaryConditions(Grids, NumberOfGrids, 
#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
			      SiblingList,
#endif
			      level, MetaData, 
			      Exterior, LevelArray[level]) == FAIL)
      return FAIL;
    JBMEM_MESSAGE(MyProcessorNumber,"jb: AfterSBC2");
    //wall_time("End SBC3");
    JBPERF_STOP_LEVEL("EL20"); // SetBoundaryConditions ()
   
    /* ------------------------------------------------------- */
    /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
       fluxes for this subgrid .
       (Note: this must be done after CorrectForRefinedFluxes). */

    JBPERF_START_LEVEL("EL27"); // Add saved fluxes to this subgrids exterior fluxes 

    for (grid = 0; grid < NumberOfGrids; grid++) {
#ifdef DC_ONLYALLOCATEMYFLUXES
      if (MyProcessorNumber ==
	  Grids[grid]->GridData->ReturnProcessorNumber()) {
#endif
      
      if (FluxCorrection==1){
	
	
	if (Grids[grid]->GridData->AddToBoundaryFluxes
	    (SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1])
	    == FAIL) {
	  fprintf(stderr, "Error in grid->AddToBoundaryFluxes.\n");
	  return FAIL;
	}
      }
      /* Delete fluxes pointed to by SubgridFluxesEstimate[subgrid]. */
             
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid]; subgrid++) {
       	
	DeleteFluxes(SubgridFluxesEstimate[grid][subgrid]);

	delete       SubgridFluxesEstimate[grid][subgrid];
      }
      delete [] SubgridFluxesEstimate[grid];
#ifdef DC_ONLYALLOCATEMYFLUXES
      }//Processor==...      
#endif
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


    //    JBPERF_STOP_LEVEL("EL25"); // EvolveLevel recursion
    //wall_time("Start RH_EL");
    if (dtThisLevelSoFar < dtLevelAbove){
      JBMEM_MESSAGE(MyProcessorNumber,"jb: ELRH");
#ifndef HAOXU
      fprintf(stderr, " EL: Rebuild Hierarchy, level %d \n", level);
#endif
      if (RebuildHierarchy(MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      JBMEM_MESSAGE(MyProcessorNumber,"jb: AfterELRH");
    }
    JBPERF_STOP_LEVEL("EL29"); // RebuildHierarchy ()
    //wall_time("End RH_EL");  

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

#ifdef JB_OPT_FAST_NEIGHBOR_SEARCH
  /* Clean up the sibling list. */

  for (grid = 0; grid < NumberOfGrids; grid++)
    delete [] SiblingList[grid].GridList;
  delete [] SiblingList;
#endif


  //=================================================================
  //  JBPERF_STOP_LEVEL("EL00") // EvolveLevel
  //=================================================================

  if( oot ) {
    //fprintf(stderr, "\n========== EL: End  ==========\n");
  }  

  
  wall_time("End EvolveLevel");
  return SUCCESS;

}

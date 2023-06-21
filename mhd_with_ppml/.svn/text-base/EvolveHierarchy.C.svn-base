
//
//
//
//
//
//
//
//
//
//
// Not the working copy.
//
//
//
//
//
//
//
//
//
//
//
//
//
/***********************************************************************
/
/  EVOLVE HIERARCHY FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  February, 1995 by GB
/              Changed to reflect changes in EvolveGrid & EvolveTopGrid.
/  modified2:  July, 1995 by GB
/              Changed to reflect new routine EvolveLevel.
/  modified3:  February, 2006 by Daniel Reynolds
/              Updated call interface to ComputePotentialFieldLevelZero
/
/  PURPOSE:
/    This routine is responsible for the evolution of the grid hierarchy.
/    It assumes the hierarchy is already constructed and the grids
/    initialized.  The routine then loops over time until one of the
/    stopping criteria is reached.  This routine also handles data dumps,
/    history dumps and restarts dumps (although the later two have not
/    yet been implemented).
/
************************************************************************/
 
#include <stdio.h>
 
#ifdef RAD_HYDRO
#include "gFLDProblem_preincludes.h"
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif
 
#include "performance.h"
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
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#ifdef RAD_HYDRO
#include "gFLDProblem.h"
#endif
 
// function prototypes
 
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

#ifdef ISO_GRAV
#ifdef RAD_HYDRO
int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior,
		GravityPotentialBoundary *PotentialBdry,
		gFLDProblem *FLDsolver);
#else
int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior,
		GravityPotentialBoundary *PotentialBdry);
#endif
#else
#ifdef RAD_HYDRO
int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior,
		gFLDProblem *FLDsolver);
#else
int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior);
#endif
#endif

#ifdef ISO_GRAV
int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, 
		 GravityPotentialBoundary *PotentialBdry, FLOAT WriteTime = -1);
#else
int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, FLOAT WriteTime = -1);
#endif
int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
int TestGravityCheckResults(LevelHierarchyEntry *LevelArray[]);
int TestGravitySphereCheckResults(LevelHierarchyEntry *LevelArray[]);
#ifdef ISO_GRAV
int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior, 
		   GravityPotentialBoundary *PotentialBdry, 
		   int &WroteData);
#else
int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior, int &WroteData);
#endif
int CheckForTimeAction(LevelHierarchyEntry *LevelArray[],
		       TopGridData &MetaData);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int OutputLevelInformation(FILE *fptr, TopGridData &MetaData,
			   LevelHierarchyEntry *LevelArray[]);
#ifdef ISO_GRAV
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], 
				   int NumberOfGrids,
				   GravityPotentialBoundary *PotentialBdry);
#else
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids);
#endif
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
				LevelHierarchyEntry *LevelArray[], int level);
float CommunicationMinValue(float Value);
int ReduceFragmentation(HierarchyEntry &TopGrid, TopGridData &MetaData,
			ExternalBoundary *Exterior,
			LevelHierarchyEntry *LevelArray[]);

#ifdef MEM_TRACE
Eint64 mused(void);
#endif
 
#define NO_REDUCE_FRAGMENTATION
 
 
 
 
#ifdef ISO_GRAV
#ifdef RAD_HYDRO
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
		    ExternalBoundary *Exterior,
		    GravityPotentialBoundary *PotentialBdry,
		    gFLDProblem &FLDsolver,
		    LevelHierarchyEntry *LevelArray[], float Initialdt)
#else
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
		    ExternalBoundary *Exterior,
		    GravityPotentialBoundary *PotentialBdry,
		    LevelHierarchyEntry *LevelArray[], float Initialdt)
#endif
#ifdef RAD_HYDRO
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
                    ExternalBoundary *Exterior,
		    gFLDProblem &FLDsolver,
		    LevelHierarchyEntry *LevelArray[], float Initialdt)
#else
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
                    ExternalBoundary *Exterior,
		    LevelHierarchyEntry *LevelArray[], float Initialdt)
#endif
#endif
{
 
  float dt;
 
  int i, Stop = FALSE, WroteData;
  double tev0, tev1;
 
  if (MetaData.Time        >= MetaData.StopTime ) Stop = TRUE;
  if (MetaData.CycleNumber >= MetaData.StopCycle) Stop = TRUE;
 
#ifdef USE_MPI
  MetaData.CPUTime = MPI_Wtime();
#endif
 
  /* Attach RandomForcingFields to BaryonFields temporarily to apply BCs */
 
  if (RandomForcing) { //AK
    LevelHierarchyEntry *Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->AppendForcingToBaryonFields();
      Temp = Temp->NextGridThisLevel;
    }
    Exterior->AppendForcingToBaryonFields();
  }
 
  /* Set top grid boundary conditions. */
 
  LevelHierarchyEntry *Temp = LevelArray[0];
 
  while (Temp != NULL) {
    if (Temp->GridData->SetExternalBoundaryValues(Exterior) == FAIL) {
      fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
      return FAIL;
    }
    if (CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, 0)
	== FAIL) {
      fprintf(stderr, "Error in CopyOverlappingZones.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }
 
 
  /* Remove RandomForcingFields from BaryonFields when BCs are set. */
 
  if (RandomForcing) { //AK
    LevelHierarchyEntry *Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->DetachForcingFromBaryonFields();
      Temp = Temp->NextGridThisLevel;
    }
    Exterior->DetachForcingFromBaryonFields();
  }
 
  /* Check for output. */
 
#ifdef ISO_GRAV
  if (CheckForOutput(&TopGrid, MetaData, Exterior, PotentialBdry, WroteData) == FAIL) {
#else
  if (CheckForOutput(&TopGrid, MetaData, Exterior, WroteData) == FAIL) {
#endif
    fprintf(stderr, "Error in CheckForOutput.\n");
    return FAIL;
  }
 
  /* Compute the acceleration field so ComputeTimeStep can find dtAccel.
     (Actually, this is a huge pain, so only do it if the
      problem really requires it). */
 
/*
  if (ProblemType == 21) {
    PrepareGravitatingMassField(&TopGrid, &MetaData, LevelArray, 0);
#ifdef ISO_GRAV
    ComputePotentialFieldLevelZero(&MetaData, Grids, NumberOfGrids, PotentialBdry);
#else
    ComputePotentialFieldLevelZero(&MetaData, Grids, NumberOfGrids);
#endif
    TopGrid.GridData->ComputeAccelerationField(GRIDS);
  }
*/
 
  /* Do the first grid regeneration. */
 
  if (RebuildHierarchy(&MetaData, LevelArray, 0) == FAIL) {
    fprintf(stderr, "Error in RebuildHierarchy.\n");
    return FAIL;
  }
 
  /* Open the OutputLevelInformation file. */
 
  FILE *LevelInfofptr;
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    LevelInfofptr = fopen("OutputLevelInformation.out", "w");
    fclose(LevelInfofptr);
  }
 
#ifdef USE_JBPERF
  Eint32 jb_iter;
#endif

  /* ====== MAIN LOOP ===== */
 
  while (!Stop) {

#ifdef USE_JBPERF
    jb_iter = MetaData.CycleNumber;
    static bool isFirstCall = true;
    if ((jb_iter % JB_ITER_PER_SEGMENT)==0 || isFirstCall) jbPerf.begin("EL");
    isFirstCall = false;
    jbPerf.attribute ("timestep",&jb_iter, JB_INT);
    jbPerf.start("EL");
    int time_sim = 1000000*MetaData.Time;
    jbPerf.assign("time-sim",time_sim);
#endif

    /* Output level information to log file. */
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      LevelInfofptr = fopen("OutputLevelInformation.out", "a");
    }

    // OutputLevelInformation() only needs to be called by all processors
    // when jbPerf is enabled.

    OutputLevelInformation(LevelInfofptr, MetaData, LevelArray);

    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fclose(LevelInfofptr);
    }
 
    /* Compute minimum timestep on the top level. */
 
    float dtProc   = huge_number;
    Temp = LevelArray[0];
 
    while (Temp != NULL) {
      dtProc = min(dtProc, Temp->GridData->ComputeTimeStep());
      Temp = Temp->NextGridThisLevel;
    }
    dt = CommunicationMinValue(dtProc);
 
    if (Initialdt != 0) {
      dt = min(dt, Initialdt);
      Initialdt = 0;
    }
 
    /* Make sure timestep doesn't go past an output. */
 
    if (ComovingCoordinates)
      for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
	if (CosmologyOutputRedshift[i] != -1)
	  dt = min(1.0001*(CosmologyOutputRedshiftTime[i]-MetaData.Time), dt);
    for (i = 0; i < MAX_TIME_ACTIONS; i++)
      if (TimeActionTime[i] > 0 && TimeActionType[i] > 0)
	dt = min(1.0001*(TimeActionTime[i] - MetaData.Time), dt);
    if (MetaData.dtDataDump > 0.0) {
      while (MetaData.TimeLastDataDump+MetaData.dtDataDump < MetaData.Time)
	MetaData.TimeLastDataDump += MetaData.dtDataDump;
     dt = min(1.0001*(MetaData.TimeLastDataDump + MetaData.dtDataDump -
		       MetaData.Time), dt);
    }
 
    /* Set the time step.  If it will cause Time += dt > StopTime, then
       set dt = StopTime - Time */
 
    dt = min(MetaData.StopTime - MetaData.Time, dt);
    Temp = LevelArray[0];
 
    while (Temp != NULL) {
      Temp->GridData->SetTimeStep(dt);
      Temp = Temp->NextGridThisLevel;
    }
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      printf("enzo-test: MetaData.CycleNumber %d\n", MetaData.CycleNumber);
      printf("enzo-test: dt %.14g\n",dt);
      printf("enzo-test: MetaData.Time %"GOUTSYM"\n", MetaData.Time);
      fflush(stdout);
    }

    if (debug) {
      printf("TopGrid dt = %"FSYM"     time = %"GOUTSYM"    cycle = %"ISYM,
	     dt, MetaData.Time, MetaData.CycleNumber);
      if (ComovingCoordinates) {
	FLOAT a, dadt;
	CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
	printf("    z = %"GOUTSYM, (1 + InitialRedshift)/a - 1);
      }
      printf("\n");
    }
 
    /* Evolve the top grid (and hence the entire hierarchy). */
 
    MPI_Barrier(MPI_COMM_WORLD);
    tev0 = MPI_Wtime();
 
#ifdef ISO_GRAV
#ifdef RAD_HYDRO
    if (EvolveLevel(&MetaData, LevelArray, 0, dt, 
		    Exterior, PotentialBdry) == FAIL) {
#else
    if (EvolveLevel(&MetaData, LevelArray, 0, dt, 
		    Exterior, PotentialBdry, &FLDsolver) == FAIL) {
#endif
#else
#ifdef RAD_HYDRO
    if (EvolveLevel(&MetaData, LevelArray, 0, dt, Exterior) == FAIL) {
#else
    if (EvolveLevel(&MetaData, LevelArray, 0, dt, 
		    Exterior, &FLDsolver) == FAIL) {
#endif
#endif
      if (MyProcessorNumber == ROOT_PROCESSOR) {
	fprintf(stderr, "Error in EvolveLevel.\n");
	fprintf(stderr, "--> Dumping data (output number %"ISYM").\n",
		MetaData.DataDumpNumber);
      }
      if (NumberOfProcessors == 1)
#ifdef ISO_GRAV
	WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
		     &TopGrid, MetaData, Exterior, PotentialBdry);
#else
	WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
		     &TopGrid, MetaData, Exterior);
#endif
      return FAIL;
    }
 
    MPI_Barrier(MPI_COMM_WORLD);
    tev1 = MPI_Wtime();
 
    FILE *evlog;
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      evlog = fopen("Evtime", "a");
      fprintf(evlog, "%8"ISYM"  %16.9e\n", MetaData.CycleNumber, tev1-tev0);
      fclose(evlog);
    }
 
    /* Rebuild the grids from level 0. */
 
    if (ProblemType != 25)
      if (RebuildHierarchy(&MetaData, LevelArray, 0) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
 
    /* Add time and check stopping criteria (steps #21 & #22)
       (note the topgrid is also keeping its own time but this statement will
       keep the two in synch). */
 
    MetaData.Time += dt;
    MetaData.CycleNumber++;
	
    if (MetaData.Time >= MetaData.StopTime) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on top grid time limit.\n");
      Stop = TRUE;
    }
    if (MetaData.CycleNumber >= MetaData.StopCycle) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on top grid cycle limit.\n");
      Stop = TRUE;
    }
 
    /* Check for time-actions. */
 
    if (CheckForTimeAction(LevelArray, MetaData) == FAIL) {
      fprintf(stderr, "Error in CheckForTimeActions.\n");
      return FAIL;
    }
 
    /* Check for output. */
 
#ifdef ISO_GRAV
    if (CheckForOutput(&TopGrid, MetaData, Exterior, PotentialBdry, WroteData) == FAIL) {
#else
    if (CheckForOutput(&TopGrid, MetaData, Exterior, WroteData) == FAIL) {
#endif
      fprintf(stderr, "Error in CheckForOutput.\n");
      return FAIL;
    }
 
    /* Try to cut down on memory fragmentation. */
 
#ifdef REDUCE_FRAGMENTATION
 
    if (WroteData && !Stop)
      if (ReduceFragmentation(TopGrid, MetaData, Exterior, LevelArray)
	  == FAIL) {
	fprintf(stderr, "Error in ReduceFragmentation.\n");
	return FAIL;
      }
 
#endif /* REDUCE_FRAGMENTATION */

#ifdef USE_JBPERF
    jbPerf.stop("EL");
    if (((jb_iter+1) % JB_ITER_PER_SEGMENT)==0) jbPerf.end("EL");
#endif

#ifdef MEM_TRACE
    Eint64 MemInUse;
    MemInUse = mused();
    fprintf(memtracePtr, "%8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif
 
  } // ===== end of main loop ====
 
#ifdef USE_JBPERF
  if (((jb_iter+1) % JB_ITER_PER_SEGMENT)!=0) jbPerf.end("EL");
  jbPerf.attribute ("timestep",0, JB_NULL);
#endif

#ifdef USE_MPI
  MetaData.CPUTime = MPI_Wtime() - MetaData.CPUTime;
#endif
 
  /* Done, so report on current time, etc. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("Time     = %9"FSYM"   CycleNumber = %6"ISYM"    Wallclock   = %9"FSYM"\n",
	   MetaData.Time, MetaData.CycleNumber, MetaData.CPUTime);
    printf("StopTime = %9"FSYM"   StopCycle   = %6"ISYM"\n",
	   MetaData.StopTime, MetaData.StopCycle);
  }
 
  /* If we are running problem 23, TestGravity, then check the results. */
 
  if (ProblemType == 23)
    TestGravityCheckResults(LevelArray);
  if (ProblemType == 25 && NumberOfProcessors == 0)
    TestGravitySphereCheckResults(LevelArray);
 
  /* if we are doing data dumps, then do one last one */
 
  if (MetaData.dtDataDump != 0.0 || MetaData.CycleSkipDataDump != 0)
#ifdef ISO_GRAV
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
		     &TopGrid, MetaData, Exterior, PotentialBdry, -666) == FAIL) {
#else
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
		     &TopGrid, MetaData, Exterior,  -666) == FAIL) {
#endif
      fprintf(stderr, "Error in WriteAllData.\n");
      return FAIL;
    }
 
  if (NumberOfProcessors > 1)
    printf("Communication: processor %"ISYM" CommunicationTime = %"FSYM"\n",
	   MyProcessorNumber, CommunicationTime);
 
  /* done */
 
  return SUCCESS;
}

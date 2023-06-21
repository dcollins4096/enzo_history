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
#include "CosmologyParameters.h"

/* function prototypes */

int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior);
int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData, 
		 ExternalBoundary *Exterior, FLOAT WriteTime = -1);
int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int TestGravityCheckResults(LevelHierarchyEntry *LevelArray[]);
int TestGravitySphereCheckResults(LevelHierarchyEntry *LevelArray[]);
int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData, 
		   ExternalBoundary *Exterior, int &WroteData);
int CheckForTimeAction(LevelHierarchyEntry *LevelArray[], 
		       TopGridData &MetaData);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int OutputLevelInformation(FILE *fptr, TopGridData &MetaData,
			   LevelHierarchyEntry *LevelArray[]);
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids);
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
				LevelHierarchyEntry *LevelArray[], int level);
float CommunicationMinValue(float Value);
float ReturnCPUTime();
int ReduceFragmentation(HierarchyEntry &TopGrid, TopGridData &MetaData,
			ExternalBoundary *Exterior, 
			LevelHierarchyEntry *LevelArray[]);

#define NO_REDUCE_FRAGMENTATION

/* EvolveHierarchy function */

int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
                    ExternalBoundary *Exterior, 
		    LevelHierarchyEntry *LevelArray[], float Initialdt)
{

  /* declarations */

  float dt;

  /* initialize */

  int i, Stop = FALSE, WroteData;
  if (MetaData.Time        >= MetaData.StopTime ) Stop = TRUE;
  if (MetaData.CycleNumber >= MetaData.StopCycle) Stop = TRUE;
  MetaData.CPUTime = ReturnCPUTime();

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

  /* Check for output. */

  if (CheckForOutput(&TopGrid, MetaData, Exterior, WroteData) == FAIL) {
    fprintf(stderr, "Error in CheckForOutput.\n");
    return FAIL;
  }

  /* Compute the acceleration field so ComputeTimeStep can find dtAccel.
     (Actually, this is a huge pain-in-the-ass, so only do it if the
      problem really requires it). */
/*
  if (ProblemType == 21) {
    PrepareGravitatingMassField(&TopGrid, &MetaData, LevelArray, 0);
    ComputePotentialFieldLevelZero(&MetaData, Grids, NumberOfGrids);
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

  /* ====== MAIN LOOP ===== */

  while (!Stop) {

    /* Output level information to log file. */

    if (MyProcessorNumber == ROOT_PROCESSOR) {
      LevelInfofptr = fopen("OutputLevelInformation.out", "a");
      OutputLevelInformation(LevelInfofptr, MetaData, LevelArray);
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

    if (debug) {
      printf("TopGrid dt = %f     time = %"GOUTSYM"    cycle = %d", 
	     dt, MetaData.Time, MetaData.CycleNumber);
      if (ComovingCoordinates) {
	FLOAT a, dadt;
	CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
	printf("    z = %"GOUTSYM, (1 + InitialRedshift)/a - 1);
      }
      printf("\n");
    }

    /* Evolve the top grid (and hence the entire hierarchy). */

    //=================================================================
    JBPERF_START("EL");
    //=================================================================

    if (EvolveLevel(&MetaData, LevelArray, 0, dt, Exterior) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR) {
	fprintf(stderr, "Error in EvolveLevel.\n");
	fprintf(stderr, "--> Dumping data (output number %d).\n",
		MetaData.DataDumpNumber);
      }
      if (NumberOfProcessors == 1)
	WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
		     &TopGrid, MetaData, Exterior);
      return FAIL;
    }

    //=================================================================
    JBPERF_STOP("EL");
    JBPERF_ADVANCE;
    JBPERF_WRITE;
    JBMEM_MESSAGE(MyProcessorNumber);
    //=================================================================

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

    if (CheckForOutput(&TopGrid, MetaData, Exterior, WroteData) == FAIL) {
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

  } // ===== end of main loop ====

  MetaData.CPUTime = ReturnCPUTime() - MetaData.CPUTime;
  
  /* Done, so report on current time, etc. */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("Time     = %9f   CycleNumber = %6d    Wallclock   = %9f\n", 
	   MetaData.Time, MetaData.CycleNumber, MetaData.CPUTime);
    printf("StopTime = %9f   StopCycle   = %6d\n",
	   MetaData.StopTime, MetaData.StopCycle);
  }

  /* If we are running problem 23, TestGravity, then check the results. */

  if (ProblemType == 23)
    TestGravityCheckResults(LevelArray);
  if (ProblemType == 25 && NumberOfProcessors == 0)
    TestGravitySphereCheckResults(LevelArray);

  /* if we are doing data dumps, then do one last one */

  if (MetaData.dtDataDump != 0.0 || MetaData.CycleSkipDataDump != 0)
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
		     &TopGrid, MetaData, Exterior) == FAIL) {
      fprintf(stderr, "Error in WriteAllData.\n");
      return FAIL;
    }

  if (NumberOfProcessors > 1)
    printf("Communication: processor %d CommunicationTime = %f\n",
	   MyProcessorNumber, CommunicationTime);

  /* done */

  return SUCCESS;
}

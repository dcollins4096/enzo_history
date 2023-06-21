/***********************************************************************
/
/  COMPUTE AND OUTPUT SOME SUMMARY INFORMATION ABOUT PERFORMANCE
/
/  written by: Greg Bryan
/  date:       August, 2003
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <df.h>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"

double ReturnWallTime();

extern double LevelWallTime[MAX_DEPTH_OF_HIERARCHY];
extern double LevelZoneCycleCount[MAX_DEPTH_OF_HIERARCHY];
extern double LevelZoneCycleCountPerProc[MAX_DEPTH_OF_HIERARCHY];
extern double ProgramStartTime;

static double CumulativeZoneCycleCount = 0;
static double LastCallTime = -1;

int OutputPerformanceInformation(FILE *fptr, TopGridData &MetaData)
{

  /* Sum up the info over processors. */

  int i, level, proc, NumberOfLevels = MaximumRefinementLevel+1;
  double AverageLevelTime[MAX_DEPTH_OF_HIERARCHY], 
    MaxLevelTime[MAX_DEPTH_OF_HIERARCHY], AverageTime[MAX_PERFORMANCE_TIMERS],
    MaxTime[MAX_PERFORMANCE_TIMERS];
  float *CommTime = NULL;
  double *TimerData = NULL;
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    CommTime = new float[NumberOfProcessors];
    TimerData = new double[NumberOfProcessors*MAX_PERFORMANCE_TIMERS];
  }

  /* Collect performance information from other processors. */

#ifdef USE_MPI

  MPI_Reduce(LevelWallTime, AverageLevelTime, NumberOfLevels, MPI_DOUBLE, 
	     MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Reduce(LevelWallTime, MaxLevelTime, NumberOfLevels, MPI_DOUBLE, 
	     MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);

  MPI_Gather(&CommunicationTime, 1, MPI_FLOAT, CommTime, 1, MPI_FLOAT,
	     ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Gather(PerformanceTimers, MAX_PERFORMANCE_TIMERS, MPI_DOUBLE, 
	     TimerData, MAX_PERFORMANCE_TIMERS, MPI_DOUBLE,
	     ROOT_PROCESSOR, MPI_COMM_WORLD);

#else /* USE_MPI */

  /* If there is no MPI, assume there is only one processor. */

  for (level = 0; level < NumberOfLevels; level++) {
    AverageLevelTime[level] = LevelWallTime[level];
    MaxLevelTime[level]     = LevelWallTime[level];
  }
  CommTime[0] = CommunicationTime;

  for (i = 0; i < MAX_PERFORMANCE_TIMERS; i++)
    TimerData[i] = PerformanceTimers[i];

#endif /* USE_MPI */

  /* Only the root processor does I/O. */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

  /* Sum up over levels. */

  double TotalZoneCycle = 0;
  double TotalLevelTime = 0;
  for (level = 0; level < NumberOfLevels; level++) {
    TotalZoneCycle     += LevelZoneCycleCount[level];
    TotalLevelTime     += LevelWallTime[level];
    AverageLevelTime[level] /= double(NumberOfProcessors);
  }
  CumulativeZoneCycleCount += TotalZoneCycle;

  for (i = 0; i < MAX_PERFORMANCE_TIMERS; i++) {
    AverageTime[i] = 0;
    MaxTime[i] = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      AverageTime[i] += TimerData[proc*MAX_PERFORMANCE_TIMERS+i];
      MaxTime[i] = max(MaxTime[i], TimerData[proc*MAX_PERFORMANCE_TIMERS+i]);
    }
    AverageTime[i] /= double(NumberOfProcessors);
  }

  /* Calculate the total wall time since last called and cumulative value. */

  double CycleTime = ReturnWallTime() - 
    ((LastCallTime < 0) ? ProgramStartTime : LastCallTime);
  LastCallTime = ReturnWallTime();
  double CumulativeCycleTime = ReturnWallTime() - ProgramStartTime;

  /* Write global info. */

  CycleTime = max(CycleTime, 1.0e-20);
  TotalLevelTime = max(TotalLevelTime, 1.0e-20);
  CumulativeCycleTime = max(CumulativeCycleTime, 1.0e-20);
  fprintf(fptr, "global: %"GOUTSYM"  %g %g %g  (EvolveLevel: %g %g)   cumulative: %g %g %g\n", MetaData.Time, 
	  TotalZoneCycle, CycleTime, TotalZoneCycle/CycleTime, 
	  TotalLevelTime, TotalZoneCycle/TotalLevelTime,
	  CumulativeZoneCycleCount, CumulativeCycleTime, 
	         CumulativeZoneCycleCount/CumulativeCycleTime);

  /* Write per level info. */

  fprintf(fptr, "level: ");
  for (level = 0; level < NumberOfLevels; level++)
    fprintf(fptr, "%d %g %g %g %g %g   ", level,
	    LevelZoneCycleCount[level], LevelWallTime[level],
	    LevelZoneCycleCount[level]/max(LevelWallTime[level], 1.0e-20),
	    MaxLevelTime[level], AverageLevelTime[level]);
  fprintf(fptr, "\n");

  /* Write per proc info. */

  fprintf(fptr, "proc: ");
  for (proc = 0; proc < NumberOfProcessors; proc++)
    fprintf(fptr, "%g ", CommTime[proc]);
  fprintf(fptr, "\n");

  /* Write timer info.
       0 - RebuildHierarchy
       1 - SolveHydroEquations
       2 - SolveRateEquations
       3 - SolveRadiativeCooling
       4 - SetBoundaryConditions (interp from parent)
       5 - PrepareDensityField
       6 - UpdateFromFinerGrids
       7 - WriteAllData
       8 - CommunicationLoadBalanceGrids
       9 - Grid_CommunicationSendRegion
       10 - Grid_CommunicationReceiveRegion
       11 - Grid_CommunicationSendRegion wait only
       12 - Grid_CommunicationSendRegion data sent (MB)
       13 - Grid_CommunicationSendRegion data received (MB)
       14 - Grid_CommunicationSendRegion number of packages send
       15 - Grid_CommunicationSendRegion number of packages received
       16 - CommunicationReceiveHandler
       17 - SetBoundaryConditions (copy from siblings)
       18 - DepositParticleMassField
       19 - PrepareGravitatingMassField1
       20 - PrepareGravitatingMassField2
       21 - CopyOverlappingMassField (within EvolveLevel)
       22 - Compute potential within EvolveLevelRoutineOptimized
       23 - Generate list of sibling grids (within EvolveLevel)
       24 - RebuildHierarchy: move particles to top
       25 - CommunicationTransferParticles
       26 - FindSubgrids (in RebuildHierarchy)
       27 - MoveSubgridParticlesFast (in RebuildHierarchy)
       28 - CommunicationShareGrids (in RebuildHierarchy)
       29 - CopyZonesFromGridCountOnly (in RebuildHierarchy)
       30 - InterpolateFieldValues (in RebuildHierarchy)
       31 - CopyZonesFromGrid (in RebuildHierarchy)
   */

  for (i = 0; i < MAX_PERFORMANCE_TIMERS; i++) {
    fprintf(fptr, "timer %d: %g %g     ", i, AverageTime[i], MaxTime[i]);
    for (proc = 0; proc < NumberOfProcessors; proc++)
      fprintf(fptr, "%g ", TimerData[proc*MAX_PERFORMANCE_TIMERS+i]);
    fprintf(fptr, "\n");
  }

  delete [] CommTime;
  delete [] TimerData;

  } // end: if (MyProcessorNumber == ROOT_PROCESSOR)

  /* Zero all the timers and counters. */

  for (level = 0; level < NumberOfLevels; level++) {
    LevelWallTime[level] = 0;
    LevelZoneCycleCount[level] = 0;
    LevelZoneCycleCountPerProc[level] = 0;
  }

  for (i = 0; i < MAX_PERFORMANCE_TIMERS; i++)
    PerformanceTimers[i] = 0;

  return SUCCESS;
}

/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
double ReturnWallTime();

#define LOAD_BALANCE_RATIO 1.05

int CommunicationLoadBalanceGrids(HierarchyEntry *GridHierarchyPointer[],
				  int NumberOfGrids)
{

  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;

  /* Initialize */

  int i, GridMemory, NumberOfCells;
  float AxialRatio, GridVolume;
  float *ComputeTime = new float[NumberOfGrids];
  float *ProcessorComputeTime = new float[NumberOfProcessors];
  double time1 = ReturnWallTime();

  /*Zhiling Lan's instrumented part */
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
  moving_count ++;
  out_count ++;
#endif

  for (i = 0; i < NumberOfProcessors; i++)
    ProcessorComputeTime[i] = 0;
  
  /* Compute work for each grid. */

  for (i = 0; i < NumberOfGrids; i++) {
    GridHierarchyPointer[i]->GridData->CollectGridInformation(GridMemory, 
				   GridVolume, NumberOfCells, AxialRatio);
    //    ComputeTime[i] = GridMemory; // roughly speaking
    ComputeTime[i] = float(NumberOfCells);
    ProcessorComputeTime[GridHierarchyPointer[i]->GridData->ReturnProcessorNumber()] += ComputeTime[i];
  }

  /* Transfer grids from heavily-loaded processors. */
  
  int Done = FALSE, MinProc = 0, MaxProc = 0;
  while (!Done) {
    
    /* Find min and max */
    
    float MaxVal = 0, MinVal = huge_number;
    for (i = 0; i < NumberOfProcessors; i++) {
      if (ProcessorComputeTime[i] > MaxVal) {
	MaxVal = ProcessorComputeTime[i];
	MaxProc = i;
      }
      if (ProcessorComputeTime[i] < MinVal) {
	MinVal = ProcessorComputeTime[i];
	MinProc = i;
      }
    }
	 
    /* Transfer a grid if the ratio is large enough. */

    if (MaxVal > LOAD_BALANCE_RATIO*MinVal) {

      /* Find a grid to transfer. */

      for (i = 0; i < NumberOfGrids; i++) {
	int proc = GridHierarchyPointer[i]->GridData->ReturnProcessorNumber();
	if (proc == MaxProc && ComputeTime[i] < 0.5*(MaxVal-MinVal)) {

	  /* Transfer. */

//	  printf("%d: moving grid %d from %d -> %d\n", MyProcessorNumber, 
//		 i, proc, MinProc);
	  GridHierarchyPointer[i]->GridData->CommunicationMoveGrid(MinProc);
//	  printf("%d: done moving grid %d\n", MyProcessorNumber, i);

	  /* Update processor compute times. */

	  ProcessorComputeTime[MaxProc] -= ComputeTime[i];
	  ProcessorComputeTime[MinProc] += ComputeTime[i];

	  break;
	}
      }

      /* If we didn't find an appropriate transfer then quit. */

      if (i == NumberOfGrids) {
	Done = TRUE;
	if (MyProcessorNumber == ROOT_PROCESSOR) {
	  printf("LB ends with MaxVal = %g (%d) and MinVal = %g (%d)\n", MaxVal,MaxProc,MinVal,MinProc);
	  for (i=0; i < NumberOfGrids; i++)
	    if (GridHierarchyPointer[i]->GridData->ReturnProcessorNumber() == MaxProc)
	      printf("%g ", ComputeTime[i]);
	  printf("\n");
	}
#ifdef MPI_INSTRUMENTATION
	/* Zhiling Lan's instrumented part */
	if (MinVal == 0)
	  timer[3] = MaxVal;
	else
	  timer[3] = MaxVal/MinVal;
#endif /* MPI_INSTRUMENTATION */
      }
    }
    else {
      Done = TRUE;
#ifdef MPI_INSTRUMENTATION
      counter[3]++;
      if (MinVal == 0)
	timer[3] = MaxVal;
      else
	timer[3] = MaxVal/MinVal;
#endif /* MPI_INSTRUMENTATION */
    }
  }

#ifdef MPI_INSTRUMENTATION
  moving_pct += float(out_count)/NumberOfGrids;
#endif /* MPI_INSTRUMENTATION */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("LoadBalance (grids=%d): ", NumberOfGrids);
    float norm = ProcessorComputeTime[0];
    for (i = 1; i < NumberOfProcessors; i++)
      norm = max(norm, ProcessorComputeTime[i]);
    //    for (i = 0; i < NumberOfProcessors; i++)
    //      ProcessorComputeTime[i] /= max(norm, 1.0e-10);
    WriteListOfFloats(stdout, NumberOfProcessors, ProcessorComputeTime);
  }

  delete [] ComputeTime;
  delete [] ProcessorComputeTime;

#ifdef MPI_INSTRUMENTATION
  /* Zhiling Lan's instrumented part */
  endtime = MPI_Wtime();
  timer[2] += endtime - starttime;
  timer[22] += (endtime - starttime)*(endtime - starttime);
  counter[2] ++;
#endif /* MPI_INSTRUMENTATION */

  PerformanceTimers[8] += ReturnWallTime() - time1;
  return SUCCESS;
}

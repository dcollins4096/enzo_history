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
/  COMMUNICATION ROUTINE: LOAD BALANCE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  David Collins, 9/22/05
/              Updated the MinValue loop to avoid problems with huge_number
/              not being huge enough.
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */
#include "performance.h"
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

//#define LOAD_BALANCE_RATIO 1.05

int CommunicationLoadBalanceGrids(HierarchyEntry *GridHierarchyPointer[],
				  int NumberOfGrids)
{

  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;

  /* Initialize */

  int i, GridMemory, CellsActive, CellsTotal;
  float AxialRatio, GridVolume;
  float *ComputeTime = new float[NumberOfGrids];
  float *ProcessorComputeTime = new float[NumberOfProcessors];

  /*Zhiling Lan's instrumented part */


  ZLAN_START;
#ifdef MPI_INSTRUMENTATION
  moving_count ++;
  out_count ++;
#endif

  for (i = 0; i < NumberOfProcessors; i++)
    ProcessorComputeTime[i] = 0;
  
  /* Compute work for each grid. */

  for (i = 0; i < NumberOfGrids; i++) {
    GridHierarchyPointer[i]->GridData->CollectGridInformation
      (GridMemory, GridVolume, CellsActive, AxialRatio, CellsTotal);
    ComputeTime[i] = GridMemory; // roughly speaking
    ProcessorComputeTime[GridHierarchyPointer[i]->GridData->ReturnProcessorNumber()] += GridMemory;

  }

  //<dbg>
  /*
  if( MyProcessorNumber == ROOT_PROCESSOR && 0 == 1){
    for (i = 0; i < NumberOfProcessors; i++) {
      fprintf(stderr,"\nLoadBalancing Initial Load, proc %d = %f\n",
	      i,ProcessorComputeTime[i]);
    }
  }
  */
  //</dbg>

  if( LoadBalancing == 0 ){
    if( MyProcessorNumber == ROOT_PROCESSOR )
      fprintf(stderr,"Warning: Load Balancing Off.\n");
    return SUCCESS;
  }

  /* Transfer grids from heavily-loaded processors. */
  
  int Done = FALSE, MinProc = 0, MaxProc = 0;
  while (!Done) {
    
    
    
    /* Find min and max */
    
    float MaxVal = 0, MinVal = huge_number;

    MaxProc = -1;
    MinProc = -1;

    //dcc 09/22/05 updated this loop to avoid huge_number being too small.
    for (i = 0; i < NumberOfProcessors; i++) {
      if (ProcessorComputeTime[i] > MaxVal) {
	MaxVal = ProcessorComputeTime[i];
	MaxProc = i;
      }
    }      
    MinVal = MaxVal;
    for (i = 0; i < NumberOfProcessors; i++) {
      if (ProcessorComputeTime[i] <= MinVal) {
	MinVal = ProcessorComputeTime[i];
	MinProc = i;
      }
    }

    if(MinProc == -1 || MaxProc == -1 )
      fprintf(stderr, "TERRIBLE ERROR: CommunicationLoadBalance unable to find processors.\n");

    /* Transfer a grid if the ratio is large enough. */

    if (MaxVal > LoadBalancingRatio*MinVal) {
      //fprintf(stderr,"CLBG Max %f Min %f MyProc %d\n", MaxVal, MinVal, MyProcessorNumber);
    /* Find a grid to transfer. */

      for (i = 0; i < NumberOfGrids; i++) {
	//fprintf(stderr,"CLBG grid %d ComputeTime = %f MyProc %d\n", i, ComputeTime[i], MyProcessorNumber);
	int proc = GridHierarchyPointer[i]->GridData->ReturnProcessorNumber();
	if (proc == MaxProc && ComputeTime[i] < 0.5*(MaxVal-MinVal)) {

	  /* Transfer. */

	  if (RandomForcing)  //AK
            GridHierarchyPointer[i]->GridData->AppendForcingToBaryonFields();


	  GridHierarchyPointer[i]->GridData->CommunicationMoveGrid(MinProc);
          if (RandomForcing)  //AK
            GridHierarchyPointer[i]->GridData->RemoveForcingFromBaryonFields();
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
#ifdef MPI_INSTRUMENTATION
	/* Zhiling Lan's instrumented part */
	timer[3] = (MinVal != 0) ? MaxVal/MinVal : MaxVal;
#endif /* MPI_INSTRUMENTATION */
      }
    }//Ratio
    else {
      Done = TRUE;
#ifdef MPI_INSTRUMENTATION
      counter[3]++;
      timer[3] = (MinVal != 0) ? MaxVal/MinVal : MaxVal;
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
    for (i = 0; i < NumberOfProcessors; i++)
      ProcessorComputeTime[i] /= max(norm, 1.0e-10);
    WriteListOfFloats(stdout, NumberOfProcessors, ProcessorComputeTime);
  }

  delete [] ComputeTime;
  delete [] ProcessorComputeTime;

  ZLAN_STOP(2);

  return SUCCESS;
}

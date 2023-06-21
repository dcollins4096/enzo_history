/***********************************************************************
/
/  COMMUNICATION ROUTINE: FIND MINIMUM VALUE AMOUNG PROCESSORS
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
#endif /* USE_MPI */
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

float ReturnCPUTime();

float CommunicationMinValue(float Value)
{

  if (NumberOfProcessors == 1)
    return Value;

  float ReturnValue = Value;

#ifdef USE_MPI

  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;

//  MPI_Datatype DataType = MPI_FLOAT;
//  if (sizeof(float) == 8)
//    DataType = MPI_DOUBLE;

//  printf("min: %d sending %f\n", MyProcessorNumber, Value);
  float time1 = ReturnCPUTime();
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif
  MPI_Allreduce(&Value, &ReturnValue, 1, DataType, MPI_MIN, MPI_COMM_WORLD);
  
#ifdef MPI_INSTRUMENTATION
  /* Zhiling Lan's instrumented part */
  endtime = MPI_Wtime();
  timer[16]+= endtime-starttime;
  counter[16] ++;
  timer[36] += (endtime-starttime)*(endtime-starttime);
  GlobalCommunication += ReturnCPUTime() - time1;
#endif /* MPI_INSTRUMENTATION */

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return ReturnValue;
}

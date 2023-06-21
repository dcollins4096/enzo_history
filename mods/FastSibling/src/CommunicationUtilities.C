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

#ifdef USE_MPI
  
int CommunicationReduceValues(float *Values, int Number, 
			      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  float time1 = ReturnCPUTime();
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  
  int i;
  float *buffer = new float[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Reduce(buffer, Values, Number, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);

  delete [] buffer;

  CommunicationTime += ReturnCPUTime() - time1;

  return SUCCESS;
}

#endif /* USE_MPI */
  
int CommunicationSumValues(float *Values, int Number)
{
#ifdef USE_MPI
  return CommunicationReduceValues(Values, Number, MPI_SUM);
#else /* USE_MPI */
  return SUCCESS;
#endif /* USE_MPI */
}

  
int CommunicationAllSumValues(float *Values, int Number)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

#ifdef USE_MPI

  float time1 = ReturnCPUTime();
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  
  int i;
  float *buffer = new float[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Allreduce(buffer, Values, Number, DataType, MPI_SUM, MPI_COMM_WORLD);

  delete [] buffer;

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
}


int CommunicationBarrier()
{
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE_MPI */
  return SUCCESS;
}

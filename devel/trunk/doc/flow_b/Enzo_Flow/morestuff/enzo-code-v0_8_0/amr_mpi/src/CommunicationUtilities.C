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
#include "performance.h"
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
#include "error.h"

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

  ZLAN_START;

  CHECK_MPI_ERROR(MPI_Allreduce(&Value, &ReturnValue, 1, DataType, MPI_MIN, 
				MPI_COMM_WORLD));

  ZLAN_STOP_GLOBAL(16);

  CommunicationTime += ReturnCPUTime() - time1;
  
#endif /* USE_MPI */

  return ReturnValue;
}
  

  
int CommunicationSumValues(float *Values, int Number)
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

  JBPERF_START("MPI_Reduce");

  CHECK_MPI_ERROR(MPI_Reduce(buffer, Values, Number, DataType, MPI_SUM, 
			     ROOT_PROCESSOR, MPI_COMM_WORLD));

  JBPERF_STOP_BYTES("MPI_Reduce",Number,DataType);

  delete [] buffer;

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
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

  JBPERF_START("MPI_Allreduce");

  CHECK_MPI_ERROR(MPI_Allreduce(buffer, Values, Number, DataType, MPI_SUM, 
				MPI_COMM_WORLD));

  JBPERF_STOP_BYTES("MPI_Allreduce",Number,DataType);

  delete [] buffer;

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
}

int CommunicationAllSumIntegerValues(int *Values, int Number)
{
  if (NumberOfProcessors == 1)
    return SUCCESS;

#ifdef USE_MPI

  float time1 = ReturnCPUTime();
  MPI_Datatype DataType = MPI_INT;

  int i;
  int *buffer = new int[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  JBPERF_START("MPI_Allreduce");

  CHECK_MPI_ERROR(MPI_Allreduce(buffer, Values, Number, DataType, MPI_SUM, 
				MPI_COMM_WORLD));

  JBPERF_STOP_BYTES("MPI_Allreduce",Number,DataType);

  delete [] buffer;

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
}


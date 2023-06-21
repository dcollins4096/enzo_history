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
#include "error.h"

/* function prototypes */

float ReturnCPUTime();

//----------------------------------------------------------------------

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
  
//----------------------------------------------------------------------
  
float CommunicationMaxValue(float Value)
{

  if (NumberOfProcessors == 1)
    return Value;

  float ReturnValue = Value;

#ifdef USE_MPI

  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;

//  printf("min: %d sending %f\n", MyProcessorNumber, Value);

#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Allreduce(&Value, &ReturnValue, 1, DataType, MPI_MAX, MPI_COMM_WORLD);
  
#ifdef MPI_INSTRUMENTATION
  /* Zhiling Lan's instrumented part */
  endtime = MPI_Wtime();
  timer[16]+= endtime-starttime;
  counter[16] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */

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

  JBPERF_START_MPI_REDUCE("MPI_Reduce",Number,DataType);
  CHECK_MPI_ERROR(MPI_Reduce(buffer, Values, Number, DataType, MPI_SUM, 
			     ROOT_PROCESSOR, MPI_COMM_WORLD));
  JBPERF_STOP_MPI_REDUCE("MPI_Reduce",Number,DataType);

  delete [] buffer;

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
}

//----------------------------------------------------------------------
  
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

  JBPERF_START_MPI_REDUCE("MPI_Allreduce",Number,DataType);
  CHECK_MPI_ERROR(MPI_Allreduce(buffer, Values, Number, DataType, MPI_SUM, 
				MPI_COMM_WORLD));
  JBPERF_STOP_MPI_REDUCE("MPI_Allreduce",Number,DataType);

  delete [] buffer;

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
}

//----------------------------------------------------------------------

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

  JBPERF_START_MPI_REDUCE("MPI_Allreduce",Number,DataType);
  CHECK_MPI_ERROR(MPI_Allreduce(buffer, Values, Number, DataType, MPI_SUM, 
				MPI_COMM_WORLD));
  JBPERF_STOP_MPI_REDUCE("MPI_Allreduce",Number,DataType);

  delete [] buffer;

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
}

//----------------------------------------------------------------------

extern "C" {
  void FORTRAN_NAME(fc_mpi_comm_rank) (int *id, int *ierr) 
  {
#ifdef USE_MPI    
    *ierr = MPI_Comm_rank (MPI_COMM_WORLD, id);
#else
    *id = 0;
    *ierr = 0;
#endif
  }
}

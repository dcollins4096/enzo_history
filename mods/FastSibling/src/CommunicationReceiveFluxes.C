/***********************************************************************
/
/  COMMUNICATION ROUTINE: RECEIVE FLUXES FROM ANOTHER PROCESSOR
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
#include "communication.h"

/* function prototypes */

float ReturnCPUTime();


int CommunicationReceiveFluxes(fluxes *Fluxes, int FromProc, 
			       int NumberOfFields, int Rank)
{

  /* Count space and allocate buffer. */

  int dim1, dim2, field, i, TotalSize = 0, Sizes[MAX_DIMENSION], TempDim;
  float *buffer = NULL;
  for (dim1 = 0; dim1 < Rank; dim1++) {
    int size = 1;
    for (dim2 = 0; dim2 < Rank; dim2++) {
      TempDim = (Fluxes->LeftFluxEndGlobalIndex[dim1][dim2] -
	         Fluxes->LeftFluxStartGlobalIndex[dim1][dim2]) + 1;
      if (dim2 == dim1)
	TempDim = 1;
      size *= TempDim;
    }
    Sizes[dim1] = size;
    TotalSize += 2*size;
  }

#ifdef USE_MPI

  TotalSize *= NumberOfFields;
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
    buffer = new float[TotalSize];

  /* receive into buffer. */

  MPI_Status status;
  MPI_Datatype DataType = MPI_FLOAT;
  if (sizeof(float) == 8)
    DataType = MPI_DOUBLE;

  float time1 = ReturnCPUTime();
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  /* If in post-receive mode, then register the receive buffer to be filled
     when the data actually arrives. */

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
    MPI_Irecv(buffer, TotalSize, DataType, FromProc, 
	     MPI_FLUX_TAG, MPI_COMM_WORLD,
	     CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
    CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
    CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
      CommunicationReceiveCurrentDependsOn;
    CommunicationReceiveIndex++;
    return SUCCESS;
  }

  /* If in send-receive mode, then wait for the data to arrive now. */

  if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
    if (MPI_Recv(buffer, TotalSize, DataType, FromProc, MPI_FLUX_TAG, 
		 MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
      fprintf(stderr, "Proc %d MPI_Recv error %d\n", MyProcessorNumber,
	      status.MPI_ERROR);
      return FAIL;
    }

#ifdef MPI_INSTRUMENTATION
  /* Zhiling lan's instrumented part */
  endtime = MPI_Wtime();
  timer[12] += endtime-starttime;
  counter[12] ++;
  timer[32] += (endtime-starttime)*(endtime-starttime);
  RecvComm += ReturnCPUTime() - time1;
#endif /* MPI_INSTRUMENTATION */

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  /* Unpack buffer */

  int index = 0;
  for (dim1 = 0; dim1 < Rank; dim1++)
    for (field = 0; field < NumberOfFields; field++) {
      for (i = 0; i < Sizes[dim1]; i++)
	Fluxes->LeftFluxes[field][dim1][i] = buffer[index++];
      for (i = 0; i < Sizes[dim1]; i++)
	Fluxes->RightFluxes[field][dim1][i] = buffer[index++];
    }

  delete [] buffer;

  return SUCCESS;
}

/***********************************************************************
/
/  GRID CLASS (RECEIVES FROM 'FAKE' GRID TO REAL GRID)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "communication.h"

/* function prototypes */

extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest, 
                                   int *sdim1, int *sdim2, int *sdim3, 
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3, 
                                   int *dstart1, int *dstart2, int *dststart3);
float ReturnCPUTime();
double ReturnWallTime();
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
                              int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);
#endif /* USE_MPI */

int grid::CommunicationReceiveRegion(grid *FromGrid, int FromProcessor, 
				     int SendField, int NewOrOld, 
				     int RegionStart[], int RegionDim[],
				     int IncludeBoundary)
{

  int i, index, field, dim, Zero[] = {0, 0, 0};
  double start_time = ReturnWallTime();

  /* Return if not on this processor. */

  if (MyProcessorNumber != ProcessorNumber &&
      MyProcessorNumber != FromProcessor)
    return SUCCESS;

  /* Compute size of region to transfer. */

  int NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
                       ((NewOrOld == NEW_AND_OLD)? 2 : 1);
  int RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];

  /* Allocate buffer. */

  int TransferSize = RegionSize*NumberOfFields;
  float *buffer = NULL;
#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
#endif /* USE_MPI */
    buffer = new float[TransferSize];

  /* If this is the from processor, pack fields. */

  int FromDim[MAX_DIMENSION], FromOffset[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    FromOffset[dim] = (dim < GridRank && IncludeBoundary == FALSE)? 
      DEFAULT_GHOST_ZONES : 0;
    FromDim[dim] = RegionDim[dim] + 2*FromOffset[dim];
  }
    
  if (MyProcessorNumber == FromProcessor) {
    
    index = 0;

    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < FromGrid->NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(FromGrid->BaryonField[field], &buffer[index],
			       FromDim, FromDim+1, FromDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       FromOffset, FromOffset+1, FromOffset+2);
	  index += RegionSize;
	}
    
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < FromGrid->NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(FromGrid->OldBaryonField[field], &buffer[index],
			       FromDim, FromDim+1, FromDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       FromOffset, FromOffset+1, FromOffset+2);
	  index += RegionSize;
      }

  }

  /* Send buffer. */

#ifdef USE_MPI

  /* only send if processor numbers are not identical */

  if (ProcessorNumber != FromProcessor) {

    MPI_Status status;

    float time1 = ReturnCPUTime();
#ifdef MPI_INSTRUMENTATION
    starttime=MPI_Wtime();
#endif /* MPI_INSTRUMENTATION */

    if (MyProcessorNumber == FromProcessor) {
//      fprintf(stderr, "RF: Sending %d floats from %d to %d\n", TransferSize, 
//	      FromProcessor, ProcessorNumber);
      CommunicationBufferedSend(buffer, TransferSize, MPI_FLOAT, 
		       ProcessorNumber, 0, MPI_COMM_WORLD, BUFFER_IN_PLACE);
    }

    if (MyProcessorNumber == ProcessorNumber) {

//      fprintf(stderr, "RF: Waiting for %d floats at %d from %d\n", 
//	      TransferSize, MyProcessorNumber, FromProcessor);

      /* Post the receive message without waiting for the message to
	 be received.  When the data arrives, this will be called again
	 in (the real) receive mode. */

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	MPI_Irecv(buffer, TransferSize, MPI_FLOAT, FromProcessor, 0, 
		  MPI_COMM_WORLD, 
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
	CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;      }

      /* If in send-receive mode, then wait for the message now. */

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
	MPI_Recv(buffer, TransferSize, MPI_FLOAT, FromProcessor, 0, 
		 MPI_COMM_WORLD, &status);
      }

    }

#ifdef MPI_INSTRUMENTATION
    /* Zhiling Lan's instrumented part */
    endtime=MPI_Wtime();
    timer[5] += endtime-starttime;
    counter[5] ++;
    timer[6] += double(TransferSize);
    timer[25] += (endtime-starttime)*(endtime-starttime);
    timer[26] += double(TransferSize*TransferSize);
    RecvComm += ReturnCPUTime() - time1;
#endif /* MPI_INSTRUMENTATION */

    CommunicationTime += ReturnCPUTime() - time1;

  }

#endif /* USE_MPI */

  /* If this is the to processor, unpack fields. */

  int GridSize = GridDimension[0]*GridDimension[1]*GridDimension[2];

  if (MyProcessorNumber == ProcessorNumber &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {
    
    index = 0;
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  if (BaryonField[field] == NULL) {
	    BaryonField[field] = new float[GridSize];
	    for (i = 0; i < GridSize; i++)
	      BaryonField[field][i] = 0;
          }
	  FORTRAN_NAME(copy3d)(&buffer[index], BaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionStart, RegionStart+1, RegionStart+2,
			       Zero, Zero+1, Zero+2);

	  index += RegionSize;
	}
    
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  if (OldBaryonField[field] == NULL) {
	    OldBaryonField[field] = new float[GridSize];
	    for (i = 0; i < GridSize; i++)
	      BaryonField[field][i] = 0;
          }
	  FORTRAN_NAME(copy3d)(&buffer[index], OldBaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionStart, RegionStart+1, RegionStart+2,
			       Zero, Zero+1, Zero+2);

	  index += RegionSize;
	}

    /* Clean up. */

    delete [] buffer;

  }

  PerformanceTimers[10] += ReturnWallTime() - start_time;

  return SUCCESS;
}


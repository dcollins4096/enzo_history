/***********************************************************************
/
/  GRID CLASS (SEND FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
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

int grid::CommunicationSendRegion(grid *ToGrid, int ToProcessor,int SendField, 
			      int NewOrOld, int RegionStart[], int RegionDim[])
{

  int index, field, dim, Zero[] = {0, 0, 0};
  double start_time = ReturnWallTime();

  /* Return if not on processor. */

  if (MyProcessorNumber != ProcessorNumber && 
      MyProcessorNumber != ToProcessor)
    return SUCCESS;

  /* Compute size of region to transfer. */

  int NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
                       ((NewOrOld == NEW_AND_OLD)? 2 : 1);
  if (SendField == ACCELERATION_FIELDS)
    NumberOfFields = GridRank;
  int RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];
  int TransferSize = RegionSize * NumberOfFields;

  /* Allocate buffer. */

  float *buffer = NULL;
#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
#endif /* USE_MPI */
    buffer = new float[TransferSize];

  /* If this is the from processor, pack fields. */

  if (MyProcessorNumber == ProcessorNumber) {

/*    printf("SendRegion: RegionStart = %d %d %d\n", RegionStart[0], 
	   RegionStart[1], RegionStart[2]); */
    
    index = 0;

    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(BaryonField[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
	}
    
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(OldBaryonField[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
      }

    if (SendField == GRAVITATING_MASS_FIELD_PARTICLES)
      FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles, buffer,
			   GravitatingMassFieldParticlesDimension,
			   GravitatingMassFieldParticlesDimension+1,
			   GravitatingMassFieldParticlesDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);

    if (SendField == GRAVITATING_MASS_FIELD)
      FORTRAN_NAME(copy3d)(GravitatingMassField, buffer,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);

    if (SendField == POTENTIAL_FIELD)
      FORTRAN_NAME(copy3d)(PotentialField, buffer,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);

    if (SendField == ACCELERATION_FIELDS)
      for (dim = 0; dim < GridRank; dim++) {
	FORTRAN_NAME(copy3d)(AccelerationField[dim], &buffer[index],
			     GridDimension, GridDimension+1, GridDimension+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     RegionStart, RegionStart+1, RegionStart+2);
	index += RegionSize;
      }
  }

  /* Send buffer. */

#ifdef USE_MPI

  /* only send if processor numbers are not identical */

  if (ProcessorNumber != ToProcessor) {
    //    fprintf(stderr, "P(%d) send %d from %d to %d (phase %d)\n",
    //	    MyProcessorNumber, TransferSize, ProcessorNumber,
    //	    ToProcessor, CommunicationDirection);

    MPI_Status status;

    float time1 = ReturnCPUTime();
#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif

    /* Send the data if on send processor, but leave buffer until the data
       has been transfered out. */

    if (MyProcessorNumber == ProcessorNumber) {
//      fprintf(stderr, "Sending %d floats from %d to %d\n", TransferSize, 
//	      MyProcessorNumber, ToProcessor);
      CommunicationBufferedSend(buffer, TransferSize, MPI_FLOAT, ToProcessor, 
	       MPI_SENDREGION_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
      PerformanceTimers[12] += double(TransferSize*sizeof(float))/1.0e6;
      PerformanceTimers[14] += 1.0;
    }

    if (MyProcessorNumber == ToProcessor) {
//      fprintf(stderr, "Waiting for %d floats at %d from %d\n", TransferSize, 
//	      MyProcessorNumber, ProcessorNumber);
      double time2 = ReturnWallTime();

      /* Post the receive message without waiting for the message to
	 be received.  When the data arrives, this will be called again
	 in (the real) receive mode. */

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	MPI_Irecv(buffer, TransferSize, MPI_FLOAT, ProcessorNumber, 
		  MPI_SENDREGION_TAG, MPI_COMM_WORLD, 
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
	CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;
      }

      /* If in send-receive mode, then wait for the message now. */

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
	MPI_Recv(buffer, TransferSize, MPI_FLOAT, ProcessorNumber, 
		 MPI_SENDREGION_TAG, MPI_COMM_WORLD, &status);
      PerformanceTimers[11] += ReturnWallTime() - time2;
      PerformanceTimers[13] += double(TransferSize*sizeof(float))/1.0e6;
      PerformanceTimers[15] += 1.0;
    }

#ifdef MPI_INSTRUMENTATION
    /* Zhiling Lan's instrumented part */
    endtime = MPI_Wtime();
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

  /* If this is the to processor, and we're either in send-receive mode
     or receive mode, then unpack the data. */

  if (MyProcessorNumber == ToProcessor &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {
    
    index = 0;
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  delete ToGrid->BaryonField[field];
	  ToGrid->BaryonField[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->BaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	}
    
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  delete ToGrid->OldBaryonField[field];
	  ToGrid->OldBaryonField[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->OldBaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	}

    if (SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
      delete ToGrid->GravitatingMassFieldParticles;
      ToGrid->GravitatingMassFieldParticles = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassFieldParticles,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }

    if (SendField == GRAVITATING_MASS_FIELD) {
      delete ToGrid->GravitatingMassField;
      ToGrid->GravitatingMassField = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassField,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }

    if (SendField == POTENTIAL_FIELD) {
      delete ToGrid->PotentialField;
      ToGrid->PotentialField = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->PotentialField,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
    
    if (SendField == ACCELERATION_FIELDS)
      for (dim = 0; dim < GridRank; dim++) {
	delete ToGrid->AccelerationField[dim];
	ToGrid->AccelerationField[dim] = new float[RegionSize];
	FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->AccelerationField[dim],
			     RegionDim, RegionDim+1, RegionDim+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     Zero, Zero+1, Zero+2);
	index += RegionSize;
      }

    /* Only delete the buffer if we're in receive mode (in send mode it will
       be deleted by CommunicationBufferedSend and if we're in post-receive
       mode then it will be deleted when we get to receive-mode). */

    delete [] buffer;
			  
  } // end: if unpack buffer

  PerformanceTimers[9] += ReturnWallTime() - start_time;

  return SUCCESS;
}


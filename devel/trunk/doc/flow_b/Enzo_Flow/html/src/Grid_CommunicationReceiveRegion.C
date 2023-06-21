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
#include "performance.h"
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
#include "error.h"

/* function prototypes */

extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest, 
                                   int *sdim1, int *sdim2, int *sdim3, 
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3, 
                                   int *dstart1, int *dstart2, int *dststart3);
float ReturnCPUTime();

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */


int grid::CommunicationReceiveRegion(grid *FromGrid, int FromProcessor, 
				     int SendField, int NewOrOld, 
				     int RegionStart[], int RegionDim[],
				     int IncludeBoundary)
{

  int i, index, field, dim, Zero[] = {0, 0, 0};

  /* Compute size of region to transfer. */

  int NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
                       ((NewOrOld == NEW_AND_OLD)? 2 : 1);
  int RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];

  /* Allocate buffer. */

  int TransferSize = RegionSize*NumberOfFields;
  float *buffer = NULL;
  if (MyProcessorNumber == FromProcessor || 
      MyProcessorNumber == ProcessorNumber)
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
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;

//    MPI_Datatype DataType = MPI_FLOAT;
//    if (sizeof(float) == 8)
//      DataType = MPI_DOUBLE;

    float time1 = ReturnCPUTime();
    ZLAN_START;

    if (MyProcessorNumber == FromProcessor) {
//      fprintf(stderr, "RF: Sending %d floats from %d to %d\n", TransferSize, 
//	      FromProcessor, ProcessorNumber);
      CommunicationBufferedSend(buffer, TransferSize, DataType, 
		       ProcessorNumber, 0, MPI_COMM_WORLD, BUFFER_IN_PLACE);
    }

    if (MyProcessorNumber == ProcessorNumber) {
//      fprintf(stderr, "RF: Waiting for %d floats at %d from %d\n", 
//	      TransferSize, MyProcessorNumber, FromProcessor);
      JBPERF_START("MPI_Recv");

      CHECK_MPI_ERROR(MPI_Recv(buffer, TransferSize, DataType, FromProcessor, 
			       0, MPI_COMM_WORLD, &status));

      JBPERF_STOP_BYTES("MPI_Recv",TransferSize,DataType);
    }

//    if (MyProcessorNumber == FromProcessor) {
//      fprintf(stderr, "RF: Sending %d floats from %d to %d\n", TransferSize,
//            FromProcessor, ProcessorNumber);
//      CHECK_MPI_ERROR(MPI_Bsend(buffer, TransferSize, DataType, 
//                      ProcessorNumber, 0, MPI_COMM_WORLD));
//    }

//    if (MyProcessorNumber == ProcessorNumber) {
//      fprintf(stderr, "RF: Waiting for %d floats at %d from %d\n",
//            TransferSize, MyProcessorNumber, FromProcessor);
//      CHECK_MPI_ERROR(MPI_Recv(buffer, TransferSize, DataType, FromProcessor,
//                      0,MPI_COMM_WORLD, &status));
//    }

    ZLAN_STOP_RECV(5);

    ZLAN_COUNT(6,TransferSize);

    CommunicationTime += ReturnCPUTime() - time1;

  }

#endif /* USE_MPI */

  /* If this is the to processor, unpack fields. */

  int GridSize = GridDimension[0]*GridDimension[1]*GridDimension[2];

  if (MyProcessorNumber == ProcessorNumber) {
    
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

  }

  /* Clean up. */

  if (MyProcessorNumber == ProcessorNumber)
    delete [] buffer;

  return SUCCESS;
}


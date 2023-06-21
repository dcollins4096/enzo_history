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
/  COPY A SET OF REGIONS BETWEEN PROCESSORS (TRANSPOSE)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "error.h"

extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest, 
                                   int *sdim1, int *sdim2, int *sdim3, 
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3, 
                                   int *dstart1, int *dstart2, int *dststart3);
extern "C" void FORTRAN_NAME(copy3dft)(float *source, float *dest, 
                                   int *sdim1, int *sdim2, int *sdim3, 
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3, 
                                   int *dstart1, int *dstart2, int *dststart3);
extern "C" void FORTRAN_NAME(copy3drt)(float *source, float *dest, 
                                   int *sdim1, int *sdim2, int *sdim3, 
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3, 
                                   int *dstart1, int *dstart2, int *dststart3);
float ReturnCPUTime();

int CommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder)
{

  /* Declarations. */

  int dim, n, i, j, size, index, Zero[] = {0,0,0};
  int LeftIndex[MAX_DIMENSION], RightIndex[MAX_DIMENSION];
  float *ReceiveBuffer, *SendBuffer;
  //  fprintf(stderr, "CT(%d): start From=%d  To=%d\n", MyProcessorNumber, 
  //	  NumberOfFromRegions, NumberOfToRegions);
							  
  region *Sends = new region[NumberOfFromRegions];
  region *Receives = new region[NumberOfToRegions];

  /* Loop over processor jumps (number of processors ahead to send). */

  for (n = 0; n < NumberOfProcessors; n++) {

    /* Copy regions into communication buffer (or just set buffer
       if there is only one FromRegion per processor). */

    int sends = 0, receives = 0;
    int SendSize = 0, ReceiveSize = 0;

    for (j = 0; j < NumberOfFromRegions; j++)
      for (i = 0; i < NumberOfToRegions; i++)
        if ((ToRegion[i].Processor - FromRegion[j].Processor +
             NumberOfProcessors) % NumberOfProcessors == n &&
	    (MyProcessorNumber == FromRegion[j].Processor ||
	     MyProcessorNumber ==   ToRegion[i].Processor)) {

	  /* Determine if there is an overlap. */

	  size = 1;
	  for (dim = 0; dim < MAX_DIMENSION; dim++) {
	    LeftIndex[dim] = max(FromRegion[j].StartIndex[dim],
                                   ToRegion[i].StartIndex[dim]);
	    RightIndex[dim] = min(
	     FromRegion[j].StartIndex[dim] + FromRegion[j].RegionDim[dim],
	       ToRegion[i].StartIndex[dim] +   ToRegion[i].RegionDim[dim])-1;
	    size *= max(RightIndex[dim] - LeftIndex[dim] + 1, 0);
	  }

	  /* If there is an overlap, add it to the list of sends/receives. */
    
	  if (size > 0) {

	    if (MyProcessorNumber == FromRegion[j].Processor) {
//	      fprintf(stderr, "CT(%d): from: i,j=%d %d  %d->%d\n", 
//		      MyProcessorNumber, j, i, FromRegion[j].Processor,
//		      ToRegion[i].Processor);
	      for (dim = 0; dim < MAX_DIMENSION; dim++) {
		Sends[sends].StartIndex[dim] = LeftIndex[dim] -
		  FromRegion[j].StartIndex[dim];
		Sends[sends].RegionDim[dim] = RightIndex[dim] - 
		  LeftIndex[dim] + 1;
	      }
	      SendSize += size;
	      Sends[sends++].Processor = j;
	    }

	    if (MyProcessorNumber == ToRegion[i].Processor) {
//	      fprintf(stderr, "CT(%d): to: %d->%d from proc %d\n", 
//		      MyProcessorNumber, j, i, FromRegion[j].Processor);
	      for (dim = 0; dim < MAX_DIMENSION; dim++) {
		Receives[receives].StartIndex[dim] = LeftIndex[dim] -
		  ToRegion[i].StartIndex[dim];
		Receives[receives].RegionDim[dim] = RightIndex[dim] - 
		  LeftIndex[dim] + 1;
	      }
	      ReceiveSize += size;
	      Receives[receives++].Processor = i;
	    }

	  } // end: if (size > 0)

        } // end: if (proc jump == n)

    /* Allocate buffer and copy data into buffer. */

    ReceiveBuffer = NULL;
    SendBuffer = new float[SendSize];

    index = 0;
//  fprintf(stderr, "CT(%d): sends = %d  SendSize = %d\n", MyProcessorNumber, 
//	    sends, SendSize);
    for (i = 0; i < sends; i++) {
      j = Sends[i].Processor;
      if (TransposeOrder == TRANSPOSE_REVERSE)
	FORTRAN_NAME(copy3drt)(FromRegion[j].Data, SendBuffer+index,
		       FromRegion[j].RegionDim, FromRegion[j].RegionDim+1,
			   FromRegion[j].RegionDim+2,
		       Sends[i].RegionDim, Sends[i].RegionDim+1, 
			   Sends[i].RegionDim+2,
		       Zero, Zero+1, Zero+2,
		       Sends[i].StartIndex, Sends[i].StartIndex+1, 
			   Sends[i].StartIndex+2);
      else
	FORTRAN_NAME(copy3d)(FromRegion[j].Data, SendBuffer+index,
		       FromRegion[j].RegionDim, FromRegion[j].RegionDim+1,
			   FromRegion[j].RegionDim+2,
		       Sends[i].RegionDim, Sends[i].RegionDim+1, 
			   Sends[i].RegionDim+2,
		       Zero, Zero+1, Zero+2,
		       Sends[i].StartIndex, Sends[i].StartIndex+1, 
			   Sends[i].StartIndex+2);
      index += Sends[i].RegionDim[0]*Sends[i].RegionDim[1]*
	       Sends[i].RegionDim[2];
    }

    /* shift buffer by n processors. */

    if (n > 0) {

      ReceiveBuffer = new float[ReceiveSize];

#ifdef USE_MPI

      MPI_Status status;
      MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;

      int ToProc = (MyProcessorNumber + n) % NumberOfProcessors;
      int FromProc = (MyProcessorNumber - n + NumberOfProcessors) %
	NumberOfProcessors;

      float time1 = ReturnCPUTime();
      ZLAN_START;
//      fprintf(stderr, "CT(%d): MPI SS/RS = %d/%d From/To = %d %d\n", 
//	      MyProcessorNumber, SendSize, ReceiveSize, FromProc, ToProc);

      JBPERF_START_MPI_SENDRECV("MPI_Sendrecv",SendSize,DataType,
				ReceiveSize,DataType);
      CHECK_MPI_ERROR(MPI_Sendrecv((void*) SendBuffer, SendSize, DataType, 
				   ToProc, MPI_TRANSPOSE_TAG, 
				   (void*) ReceiveBuffer, ReceiveSize, 
				   DataType, FromProc, MPI_TRANSPOSE_TAG, 
				   MPI_COMM_WORLD, &status));
      JBPERF_STOP_MPI_SENDRECV("MPI_Sendrecv",SendSize,DataType,
			       ReceiveSize,DataType);

      ZLAN_STOP(14);

      CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

    } else {
      ReceiveBuffer = SendBuffer;
      SendBuffer = NULL;
    }

    /* Copy from communication buffer back to regions. */

    index = 0;
//    fprintf(stderr, "CT(%d): receives = %d\n", MyProcessorNumber, receives);
    for (i = 0; i < receives; i++) {
      j = Receives[i].Processor;
      if (ToRegion[j].Data == NULL)
	ToRegion[j].Data = new float[ToRegion[j].RegionDim[0]*
	      ToRegion[j].RegionDim[1]*ToRegion[j].RegionDim[2]];
      if (TransposeOrder == TRANSPOSE_FORWARD)
	FORTRAN_NAME(copy3dft)(ReceiveBuffer+index, ToRegion[j].Data,
		       Receives[i].RegionDim, Receives[i].RegionDim+1, 
			   Receives[i].RegionDim+2,
		       ToRegion[j].RegionDim, ToRegion[j].RegionDim+1,
			   ToRegion[j].RegionDim+2,
		       Receives[i].StartIndex, Receives[i].StartIndex+1, 
			   Receives[i].StartIndex+2,
		       Zero, Zero+1, Zero+2);
      else
	FORTRAN_NAME(copy3d)(ReceiveBuffer+index, ToRegion[j].Data,
		       Receives[i].RegionDim, Receives[i].RegionDim+1, 
			   Receives[i].RegionDim+2,
		       ToRegion[j].RegionDim, ToRegion[j].RegionDim+1,
			   ToRegion[j].RegionDim+2,
		       Receives[i].StartIndex, Receives[i].StartIndex+1, 
			   Receives[i].StartIndex+2,
		       Zero, Zero+1, Zero+2);
      index += Receives[i].RegionDim[0]*Receives[i].RegionDim[1]*
	       Receives[i].RegionDim[2];
    }

    /* Clean up. */

//    fprintf(stderr, "CT(%d): end jump %d\n", MyProcessorNumber, n);
    delete [] SendBuffer;
    delete [] ReceiveBuffer;

  } // end: loop over processors jumps

  /* Delete FromRegion data. */

//  fprintf(stderr, "CT(%d): Deleting FromRegions\n", MyProcessorNumber);
  for (i = 0; i < NumberOfFromRegions; i++) {
    delete [] FromRegion[i].Data;
    FromRegion[i].Data = NULL;
  }

  /* Clean up. */

  delete [] Sends;
  delete [] Receives;

  return SUCCESS;
}

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
/  COMMUNICATION ROUTINE: BUFFERED SEND WITH SELF-BUFFERING
/
/  written by: Greg Bryan
/  date:       January, 2001
/  modified1:
/
/  PURPOSE:
/    A replacement for MPI_Bsend, this routine allocates a buffer if
/     the BufferSize is non-negative and copies the data into that
/     buffer (otherwise the inbuffer is used directly).
/
************************************************************************/

#ifdef USE_MPI

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "error.h"

/* Records the number of times we've been called. */

static int CallCount = 0;

/* Defines the number of calls to wait before scanning. */

#define NUMBER_OF_CALLS_BETWEEN_SCANS 30

/* The MPI Handle and buffer storage area. */

#define MAX_NUMBER_OF_MPI_BUFFERS 10000

static MPI_Request  RequestHandle[MAX_NUMBER_OF_MPI_BUFFERS];
static char        *RequestBuffer[MAX_NUMBER_OF_MPI_BUFFERS];
static int          LastActiveIndex = -1;


/* function prototypes */

int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize)
{

  int i, RequestDone;
  MPI_Status Status;
  void *buffer_send;

  /* First, check to see if we should do a scan. */

  if (++CallCount % NUMBER_OF_CALLS_BETWEEN_SCANS == 0) {

    int NewLastActiveIndex = -1;
    for (i = 0; i < LastActiveIndex+1; i++) {
      if (RequestBuffer[i] != NULL) {
	CHECK_MPI_ERROR(MPI_Test(RequestHandle+i, &RequestDone, &Status));
	if (RequestDone) {

	  /* If the request is done, deallocate associated buffer. */

	  delete [] RequestBuffer[i];
	  RequestBuffer[i] = NULL;

	} else
	  NewLastActiveIndex = max(i, NewLastActiveIndex);
      }
    } // end: loop over request handles
    
  }

  /* If necessary, allocate buffer. */

  if (BufferSize != BUFFER_IN_PLACE) {
    buffer_send = new char[BufferSize];
    memcpy(buffer_send, buffer, BufferSize);
  }
  else
    buffer_send = (void *) buffer;
    
  /* Find open spot. */

  int index = LastActiveIndex+1;
  for (i = 0; i < LastActiveIndex+1; i++)
    if (RequestBuffer[i] == NULL)
      index = i;

  /* Error check. */

  if (index >= MAX_NUMBER_OF_MPI_BUFFERS-1) {
    fprintf(stderr, "CommunicationBufferedSend: increase MAX_NUMBER_OF_MPI_BUFFERs\n");
    exit(EXIT_FAILURE);
  }

  /* call MPI send and store handle. */

  JBPERF_START_MPI_SEND("MPI_Isend",size,Type);

  CHECK_MPI_ERROR(MPI_Isend(buffer_send, size, Type, Target, Tag, CommWorld, 
			    RequestHandle+index));

  JBPERF_STOP_MPI_SEND("MPI_Isend",size,Type);

  /* Store buffer info. */

  RequestBuffer[index] = (char *) buffer_send;
  LastActiveIndex = max(LastActiveIndex, index);

  return SUCCESS;
}

#endif /* USE_MPI */

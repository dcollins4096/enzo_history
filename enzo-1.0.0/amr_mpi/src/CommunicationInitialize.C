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
/  COMMUNICATION ROUTINE: INITIALIZE
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
#include <stdlib.h>
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

int CommunicationInitialize(int *argc, char **argv[])
{

#ifdef USE_MPI

  /* Initialize MPI and get info. */

  CHECK_MPI_ERROR(MPI_Init(argc, argv));
  CHECK_MPI_ERROR(MPI_Comm_rank(MPI_COMM_WORLD, &MyProcessorNumber));
  CHECK_MPI_ERROR(MPI_Comm_size(MPI_COMM_WORLD, &NumberOfProcessors));

#ifdef USE_MPE
  MPE_Init_log();
  MPE_Start_log();

  jb::mpi_send_start = MPE_Log_get_event_number();
  jb::mpi_send_stop  = MPE_Log_get_event_number();
  jb::mpi_recv_start = MPE_Log_get_event_number();
  jb::mpi_recv_stop  = MPE_Log_get_event_number();
  jb::mpi_sendrecv_start = MPE_Log_get_event_number();
  jb::mpi_sendrecv_stop  = MPE_Log_get_event_number();
  jb::mpi_barrier_start = MPE_Log_get_event_number();
  jb::mpi_barrier_stop  = MPE_Log_get_event_number();
  jb::mpi_gather_start = MPE_Log_get_event_number();
  jb::mpi_gather_stop  = MPE_Log_get_event_number();
  jb::mpi_reduce_start = MPE_Log_get_event_number();
  jb::mpi_reduce_stop  = MPE_Log_get_event_number();
  if (MyProcessorNumber == 0) {
    MPE_Describe_state(jb::mpi_send_start,jb::mpi_send_stop,
		       "MPI Send", "red");
    MPE_Describe_state(jb::mpi_recv_start,jb::mpi_recv_stop,
		       "MPI Recv", "blue");
    MPE_Describe_state(jb::mpi_sendrecv_start,jb::mpi_sendrecv_stop,
		       "MPI Sendrecv", "magenta");
    MPE_Describe_state(jb::mpi_barrier_start,jb::mpi_barrier_stop,
		       "MPI Barrier", "yellow");
    MPE_Describe_state(jb::mpi_gather_start,jb::mpi_gather_stop,
		       "MPI Gather", "cyan");
    MPE_Describe_state(jb::mpi_reduce_start,jb::mpi_reduce_stop,
		       "MPI Reduce", "green");
  }
#endif

  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("MPI_Init: NumberOfProcessors = %d\n", NumberOfProcessors);

#else /* USE_MPI */

  MyProcessorNumber  = 0;
  NumberOfProcessors = 1;
  
#endif /* USE_MPI */

  CommunicationTime = 0;

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  return SUCCESS;
}

int CommunicationFinalize()
{

#ifdef USE_MPI
#ifdef USE_MPE
  CHECK_MPI_ERROR(MPE_Finish_log("mpe"));
#endif
  CHECK_MPI_ERROR(MPI_Finalize());
#endif /* USE_MPI */

  return SUCCESS;
}

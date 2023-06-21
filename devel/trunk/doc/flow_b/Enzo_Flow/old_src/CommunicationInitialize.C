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
#ifdef USE_MPI
#include "mpi.h"
#include <stdlib.h>
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

int CommunicationInitialize(int *argc, char **argv[])
{

#ifdef USE_MPI

  /* Initialize MPI and get info. */

  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyProcessorNumber);
  MPI_Comm_size(MPI_COMM_WORLD, &NumberOfProcessors);

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
  MPI_Finalize();
#endif /* USE_MPI */

  return SUCCESS;
}

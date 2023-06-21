/***********************************************************************
/
/  COMMUNICATION ROUTINE: BROADCAST A VALUE FROM ROOT TO OTHERS
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

/* function prototypes */

float ReturnCPUTime();

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor)
{

  if (NumberOfProcessors == 1)
    return SUCCESS;

#ifdef USE_MPI

  float time1 = ReturnCPUTime();
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif
  MPI_Bcast((void*) Value, 1, MPI_INT, BroadcastProcessor, MPI_COMM_WORLD);
  
#ifdef MPI_INSTRUMENTATION
  /* Zhiling Lan's instrumented part */
  endtime = MPI_Wtime();
  timer[15]+= endtime-starttime;
  counter[15] ++;
  timer[35] += (endtime-starttime)*(endtime-starttime);
  GlobalCommunication += ReturnCPUTime() - time1;
#endif /* MPI_INSTRUMENTATION */

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
}

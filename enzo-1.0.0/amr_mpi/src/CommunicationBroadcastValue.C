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

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor)
{

  if (NumberOfProcessors == 1)
    return SUCCESS;

#ifdef USE_MPI

  float time1 = ReturnCPUTime();

  ZLAN_START;

  CHECK_MPI_ERROR(MPI_Bcast((void*) Value, 1, MPI_INT, 
			    BroadcastProcessor, MPI_COMM_WORLD));
  
  ZLAN_STOP_GLOBAL(15);

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  return SUCCESS;
}

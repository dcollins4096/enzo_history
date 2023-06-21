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
/  COMMUNICATION ROUTINE: UPDATE STAR PARTICLE COUNT
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modified1:
/
/  PURPOSE:
/    This routines keeps track of the number of new star particles.
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
#include "StarParticleData.h"
#include "error.h"

/* function prototypes */

float ReturnCPUTime();


int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids)
{

  int grid, *TotalParticleCount = new int[NumberOfGrids],
          *PartialParticleCount = new int[NumberOfGrids];

  /* Set ParticleCount to zero and record number of particles for grids
     on this processor. */

  for (grid = 0; grid < NumberOfGrids; grid++) {
    TotalParticleCount[grid] = 0;
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber)
      PartialParticleCount[grid] = 
	Grids[grid]->GridData->ReturnNumberOfParticles();
    else
      PartialParticleCount[grid] = 0;
  }

#ifdef USE_MPI

  /* Get counts from each processor to get total list of new particles. */

  float time1 = ReturnCPUTime();
  ZLAN_START;

  JBPERF_START_MPI_REDUCE("MPI_Allreduce",NumberOfGrids, MPI_INT);

  CHECK_MPI_ERROR(MPI_Allreduce(PartialParticleCount, TotalParticleCount, 
				NumberOfGrids, MPI_INT, MPI_SUM, 
				MPI_COMM_WORLD));
  JBPERF_STOP_MPI_REDUCE("MPI_Allreduce",NumberOfGrids, MPI_INT);

  ZLAN_STOP_GLOBAL(11);

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  /* Set new particle count for each grid. */

  for (grid = 0; grid < NumberOfGrids; grid++) {

    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber)

      /* If this grid is on this processor, then call routine to set the
	 particle index numbers.  This also updates NumberOfStarParticles. */

      Grids[grid]->GridData->SetNewParticleIndex(NumberOfStarParticles,
						 MetaData->NumberOfParticles);

    else {

      /* If not on this processor, then keep track of the number of new
	 star particles (which is the difference between the number of
	 got from the communication and what is currently stored).
	 Finally, correct the number of particles in our record. */

      NumberOfStarParticles += TotalParticleCount[grid] - 
                       Grids[grid]->GridData->ReturnNumberOfParticles();
      Grids[grid]->GridData->SetNumberOfParticles(TotalParticleCount[grid]);


    }

  }

  //  printf("NumberOfStarParticles = %d\n", NumberOfStarParticles);

  /* Clean up. */

  delete [] TotalParticleCount;
  delete [] PartialParticleCount;

  return SUCCESS;
}

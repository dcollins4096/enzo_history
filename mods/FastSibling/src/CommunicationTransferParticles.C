/***********************************************************************
/
/  COMMUNICATION ROUTINE: TRANSFER PARTICLES
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

#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_ParticleMoveList;
#endif


int CommunicationTransferParticles(grid *GridPointer[], int NumberOfGrids)
{

  /* Declarations. */

  struct ParticleMoveList {
    int FromGrid;
    int ToGrid[6];
    int NumberToMove[6];
    float_int *Pointer[6];
  };

  /* This is a loop that keeps going until no particles were moved. */

  int NumberOfParticlesMoved, Done = FALSE;
  while (!Done) {
  
  NumberOfParticlesMoved = 0;

  int i, j, grid, GridsToSend = 0;
  ParticleMoveList *SendList = new ParticleMoveList[NumberOfGrids];
  for (grid = 0; grid < NumberOfGrids; grid++)
    for (i = 0; i < 6; i++) {
      SendList[grid].NumberToMove[i] = 0;
      SendList[grid].Pointer[i]      = NULL;
      SendList[grid].ToGrid[i]       = -1;
    }

  /*  for (j = 0; j < NumberOfGrids; j++)
      printf("Pa(%d) CTP grid[%d] = %d\n", MyProcessorNumber, j, GridPointer[j]->ReturnNumberOfParticles()); */

  /* -------------------------------------------------------------------- */
  /* Generate the list of particle moves. */

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (GridPointer[grid]->ReturnProcessorNumber() == MyProcessorNumber) {
      SendList[GridsToSend].FromGrid = grid;
      if (GridPointer[grid]->CommunicationTransferParticles(GridPointer,
	     NumberOfGrids, 
             SendList[GridsToSend].ToGrid, 
             SendList[GridsToSend].NumberToMove, 
             SendList[GridsToSend].Pointer, 
	     COPY_OUT) == FAIL) {
	fprintf(stderr, "Error in grid->CommunicationTransferParticless\n");
	return FAIL;
      }
      GridsToSend++;
    }


  /* -------------------------------------------------------------------- */
  /* Share the Particle moves. */

  /* Allocate the array to receive subgrids. */

  ParticleMoveList *SharedList = NULL;
#ifdef USE_MPI
  int NumberOfSharedGrids = 0;
#endif /* USE_MPI */

  if (NumberOfProcessors > 1) {

#ifdef USE_MPI

    /* Generate a new MPI type corresponding to the ParticleMoveList struct. */

    MPI_Status status;

    if (FirstTimeCalled) {
      MPI_Type_contiguous(sizeof(ParticleMoveList), MPI_BYTE, 
			  &MPI_ParticleMoveList);
      MPI_Type_commit(&MPI_ParticleMoveList);
      FirstTimeCalled = FALSE;
    }
    int *SendListCount = new int[NumberOfProcessors],
        *SendListDisplacements = new int[NumberOfProcessors];

    /* Get counts from each processor to allocate buffers. */

    float time1 = ReturnCPUTime();
#ifdef MPI_INSTRUMENTATION
    starttime=MPI_Wtime();
#endif /* MPI_INSTRUMENTATION */
    MPI_Allgather(&GridsToSend, 1, MPI_INT, SendListCount, 1, MPI_INT,
		  MPI_COMM_WORLD);

    /* Allocate buffers and generated displacement list. */

    for (i = 0; i < NumberOfProcessors; i++) {
      SendListDisplacements[i] = NumberOfSharedGrids;
      NumberOfSharedGrids += SendListCount[i];
    }
    if (NumberOfSharedGrids != NumberOfGrids) {
      fprintf(stderr, "CTP error\n");
      return FAIL;
    }
    SharedList = new ParticleMoveList[NumberOfSharedGrids];

    /* Perform sharing operation. */

    MPI_Allgatherv(SendList, GridsToSend, MPI_ParticleMoveList, SharedList,
		   SendListCount, SendListDisplacements, 
		   MPI_ParticleMoveList, MPI_COMM_WORLD);

#ifdef MPI_INSTRUMENTATION
    /* Zhiling Lan's instrumented part */
    endtime = MPI_Wtime();
    timer[9] += endtime-starttime;
    counter[9] ++;
    timer[10] += double(NumberOfSharedGrids);
    timer[30] += double(NumberOfSharedGrids*NumberOfSharedGrids);
    timer[29] += (endtime-starttime)*(endtime-starttime);
    GlobalCommunication += ReturnCPUTime() - time1;
#endif /* MPI_INSTRUMENTATION */

    CommunicationTime += ReturnCPUTime() - time1;
    
    delete [] SendListCount;
    delete [] SendListDisplacements;
  
#endif /* USE_MPI */

    /* -------------------------------------------------------------------- */
    /* Move particles. */

    for (j = 0; j < NumberOfGrids; j++) {

      int FromGrid = SharedList[j].FromGrid;

      /* Loop over transfer directions. */

      for (i = 0; i < 6; i++) {

	int ToGrid = max(SharedList[j].ToGrid[i], 0);
	int FromProcessor = GridPointer[FromGrid]->ReturnProcessorNumber();
	int ToProcessor = GridPointer[ToGrid]->ReturnProcessorNumber();
	int TransferSize = SharedList[j].NumberToMove[i]*
	                   (9+NumberOfParticleAttributes);

	NumberOfParticlesMoved += SharedList[j].NumberToMove[i];
	if (TransferSize > 0 && FromProcessor != ToProcessor) {

#ifdef USE_MPI	

	  /* Send particles (Note: this is not good -- the data transfer will not work
                across heterogeneous machines). */

	  time1 = ReturnCPUTime();
	  if (FromProcessor == MyProcessorNumber)
	    MPI_Send(SharedList[j].Pointer[i], TransferSize*sizeof(float_int), MPI_BYTE, 
		     ToProcessor, MPI_TRANSFERPARTICLE_TAG, MPI_COMM_WORLD);

	  /* Receive particles. */

	  if (ToProcessor == MyProcessorNumber) {
	    SharedList[j].Pointer[i] = new float_int[TransferSize];
	    MPI_Recv(SharedList[j].Pointer[i], TransferSize*sizeof(float_int), MPI_BYTE,
		     FromProcessor, MPI_TRANSFERPARTICLE_TAG, MPI_COMM_WORLD, 
		     &status);
	  }
#ifdef MPI_INSTRUMENTATION
	  RecvComm += ReturnCPUTime() - time1;
#endif  /* MPI_INSTRUMENTATION */  
	  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

	} // end: if (TransferSize > 0...

	/* If this is not on my processor, then clean up pointer, since it
	   doesn't make sense on this processor (or if nothing moved). */

	if ((MyProcessorNumber != FromProcessor && 
	     MyProcessorNumber != ToProcessor) || TransferSize == 0)
	  SharedList[j].Pointer[i] = NULL;

      } // end: loop over i (directions)
    } // end: loop over j (grids)
    
  } else
    SharedList = SendList;  // if there is only one processor

  /* -------------------------------------------------------------------- */
  /* Copy particles back to grids. */

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (GridPointer[grid]->ReturnProcessorNumber() == MyProcessorNumber) {

      int LocalNumberToMove[6], counter = 0;
      float_int *LocalPointer[6];
      for (i = 0; i < 6; i++)
	LocalNumberToMove[i] = 0;

      /* Extract grid lists to put into this grid. */

      for (j = 0; j < NumberOfGrids; j++)
	for (i = 0; i < 6; i++)
	  if (SharedList[j].ToGrid[i] == grid) {
	    LocalNumberToMove[counter] = SharedList[j].NumberToMove[i];
	    LocalPointer[counter++] = SharedList[j].Pointer[i];
	  }

      /* Copy particles back (SharedList[grid].ToGrid is a dummy, not used). */

      if (GridPointer[grid]->CommunicationTransferParticles(GridPointer,
	      NumberOfGrids, SharedList[grid].ToGrid, LocalNumberToMove, 
	      LocalPointer, COPY_IN) == FAIL) {
	fprintf(stderr, "Error in grid->CommunicationTransferParticless\n");
	return FAIL;
      }

    } // end: if grid is on my processor

  /* Set number of particles so everybody agrees. */

  if (NumberOfProcessors > 1) {
    int *Changes = new int[NumberOfGrids];
    for (j = 0; j < NumberOfGrids; j++)
      Changes[j] = 0;
    for (j = 0; j < NumberOfGrids; j++)
      for (i = 0; i < 6; i++) 
	if (SharedList[j].ToGrid[i] != -1) {
	  Changes[SharedList[j].FromGrid] -= SharedList[j].NumberToMove[i];
	  Changes[SharedList[j].ToGrid[i]] += SharedList[j].NumberToMove[i];
	}
    for (j = 0; j < NumberOfGrids; j++) {
      if (GridPointer[j]->ReturnProcessorNumber() != MyProcessorNumber)
	GridPointer[j]->SetNumberOfParticles(
		         GridPointer[j]->ReturnNumberOfParticles()+Changes[j]);
      //      printf("Pb(%d) CTP grid[%d] = %d\n", MyProcessorNumber, j, GridPointer[j]->ReturnNumberOfParticles());
    }
    delete [] Changes;
  }

  /* -------------------------------------------------------------------- */
  /* CleanUp. */

  for (j = 0; j < NumberOfGrids; j++)
    for (i = 0; i < 6; i++) {
      delete [] SharedList[j].Pointer[i];
    }

  if (SendList != SharedList)
    delete [] SendList;

  delete [] SharedList;

  /* Check for completion. */

  if (debug)
    printf("CommunicationTransferParticles: moved = %d\n",
	   NumberOfParticlesMoved);
  if (NumberOfParticlesMoved == 0)
    Done = TRUE;

  }

  return SUCCESS;
}

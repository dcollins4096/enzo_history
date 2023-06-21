/***********************************************************************
/
/  EVOLVE LEVEL ROUTINES (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: This routine processes the receives stored in the 
/           CommunicationReceive stack.  Each receive is tagged with a 
/           type which indicates which method to call 
/           (and a record of the arguments).
/
************************************************************************/

#include <stdio.h>
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
 
#ifdef USE_MPI
static int ListOfIndices[MAX_RECEIVE_BUFFERS];
static MPI_Status ListOfStatuses[MAX_RECEIVE_BUFFERS];
#endif /* USE_MPI */
float ReturnCPUTime();
double ReturnWallTime();

int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[],
				int NumberOfSubgrids[])
{

#ifdef USE_MPI

  /* Set the communication mode. */

  CommunicationDirection = COMMUNICATION_RECEIVE;

  //  printf("P(%d) in CRH with %d requests\n", MyProcessorNumber,
  //	 CommunicationReceiveIndex);

  double start_time = ReturnWallTime();
  int ReceivesCompletedToDate = 0, NumberOfCompleteRequests, index, errcode,
    TotalReceives, igrid, isubgrid, dim;
  FLOAT EdgeOffset[MAX_DIMENSION];
  grid *grid_one, *grid_two;
  TotalReceives = CommunicationReceiveIndex;

  /* Define a temporary flux holder for the refined fluxes. */

  fluxes SubgridFluxesRefined;

  while (ReceivesCompletedToDate < TotalReceives) {

    /* Call the MPI wait handler. */

    float time1 = ReturnCPUTime();

    MPI_Waitsome(TotalReceives, CommunicationReceiveMPI_Request,
		 &NumberOfCompleteRequests, ListOfIndices, ListOfStatuses);
    //    printf("MPI: %d %d %d\n", TotalReceives, ReceivesCompletedToDate,
    //    	   NumberOfCompleteRequests);

    CommunicationTime += ReturnCPUTime() - time1;

    /* Error check */

#if 0
    if (NumberOfCompleteRequests == MPI_UNDEFINED) {
      fprintf(stderr, "Error in MPI_Waitsome\n");
      return FAIL;
    }
#endif

    /* Should loop over newly received completions and check error msgs now. */
    for (index = 0; index < NumberOfCompleteRequests; index++)
      if (ListOfStatuses[index].MPI_ERROR != 0)
	fprintf(stderr, "P(%d) mpi error %d\n", MyProcessorNumber,
		ListOfStatuses[index].MPI_ERROR);

    /* Loop over the receive handles, looking for completed (i.e. null)
       requests associated with unprocessed (i.e. non-null) grids. 
       It's insufficient to just loop over newly completed receives because
       there may be some completed receives which were not processed due
       to dependence issues. */

    for (index = 0; index < TotalReceives; index++) {
      //      fprintf(stderr, "%d %d %d %d %d\n", index, 
      //	      CommunicationReceiveCallType[index],
      //      	      CommunicationReceiveGridOne[index],
      //      	      CommunicationReceiveMPI_Request[index],
      //      	      CommunicationReceiveDependsOn[index]);
      if (CommunicationReceiveGridOne[index] != NULL &&
	  CommunicationReceiveMPI_Request[index] == MPI_REQUEST_NULL) {

	/* If this depends on an un-processed receive, then skip it. */

	if (CommunicationReceiveDependsOn[index] != 
	    COMMUNICATION_NO_DEPENDENCE)
	  if (CommunicationReceiveGridOne[CommunicationReceiveDependsOn[index]]
	      != NULL)
	    continue;

	grid_one = CommunicationReceiveGridOne[index];
	grid_two = CommunicationReceiveGridTwo[index];
	CommunicationReceiveIndex = index;

	/* Copy out the argument if needed */

	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  EdgeOffset[dim] = CommunicationReceiveArgument[dim][index];
	  
	/* Handle the buffers received, calling the appropriate method. */

	switch (CommunicationReceiveCallType[index]) {

	case 1:
	  errcode = grid_one->InterpolateBoundaryFromParent(grid_two);
	  break;

	case 2:
	  errcode = grid_one->CopyZonesFromGrid(grid_two, EdgeOffset);
	  break;

	case 3:
	  errcode = grid_one->DepositParticlePositions(grid_two,
			CommunicationReceiveArgument[0][index],
			CommunicationReceiveArgumentInt[0][index]);
	  break;

	case 4:
	  errcode = grid_one->CopyParentToGravitatingFieldBoundary(grid_two);
	  break;

	case 5:
	  errcode = grid_one->DepositBaryons(grid_two,
			           CommunicationReceiveArgument[0][index]);
	  break;

	case 6:
	  errcode = grid_one->AddOverlappingParticleMassField(grid_two,
							      EdgeOffset);
	  break;

	case 7:
	  errcode = grid_one->PreparePotentialField(grid_two);
	  break;

	case 8:
	  errcode = grid_one->CopyOverlappingMassField(grid_two, EdgeOffset);
	  break;

	case 9:
	  errcode = grid_one->CopyPotentialField(grid_two, EdgeOffset);
	  break;

	case 10:
	  errcode = grid_one->InterpolateAccelerations(grid_two);
	  break;
      
	case 11:  /* Note this one involves two calls. */

	  /* Project subgrid's refined fluxes to the level of this grid. */

	  if (grid_one->GetProjectedBoundaryFluxes(grid_two, 
					       SubgridFluxesRefined) == FAIL) {
	    fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	    return FAIL;
	  }
	
	  /* Correct this grid for the refined fluxes (step #19)
	     (this also deletes the fields in SubgridFluxesRefined). */

	  igrid = CommunicationReceiveArgumentInt[0][index];
	  isubgrid = CommunicationReceiveArgumentInt[1][index];
	  if ((errcode = grid_two->CorrectForRefinedFluxes
	      (SubgridFluxesEstimate[igrid][isubgrid], &SubgridFluxesRefined, 
	       SubgridFluxesEstimate[igrid][NumberOfSubgrids[igrid] - 1]     ))
	      == FAIL) {
	    fprintf(stderr, "Error in grid->CorrectForRefinedFluxes.\n");
	    return FAIL;
	  }
	  break;

	case 12:
	  errcode = grid_one->ProjectSolutionToParentGrid(*grid_two);
	  break;
      

	default:
	  fprintf(stderr, "Unrecognized call type %d\n", 
		  CommunicationReceiveCallType[index]);
	  return FAIL;

	} // end: switch on call type

	/* Report error if there has been one in any of the above calls. */

	if (errcode == FAIL) {
	  fprintf(stderr, "Error in CommunicationReceiveHandler, method %d\n",
		  CommunicationReceiveCallType[index]);
	  return FAIL;
	}

	/* Mark this receive complete. */

	CommunicationReceiveGridOne[index] = NULL;
	ReceivesCompletedToDate++;

      } // end: if statement to check if receive should be processed

    } // end: loop over all receives

  } // end: while loop waiting for all receives to be processed

  CommunicationReceiveIndex = 0;

  /* Reset the communication mode. */

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
  //  printf("P(%d) out of CRH\n", MyProcessorNumber);
  PerformanceTimers[16] += ReturnWallTime() - start_time;

#endif /* USE_MPI */

  return SUCCESS;

}

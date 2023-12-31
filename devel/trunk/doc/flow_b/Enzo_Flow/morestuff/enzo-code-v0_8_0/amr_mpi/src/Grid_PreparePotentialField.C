/***********************************************************************
/
/  GRID CLASS (PREPARE FIRST GUESS TO GRAVITATIONAL POTENTIAL)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#define SPLINE

/* function prototypes */

#ifdef SPLINE

extern "C" void FORTRAN_NAME(int_spline)(
                                      float *source, float *dest, int *ndim, 
				      int *sdim1, int *sdim2, int *sdim3,
				      int *ddim1, int *ddim2, int *ddim3, 
				      int *start1, int *start2, int *start3,
				      int *refine1, int *refine2, int *refine3,
				      float *temp1, float *temp2);
#else /* SPLINE */

extern "C" void FORTRAN_NAME(prolong)(float *source, float *dest, int *ndim, 
				      int *sdim1, int *sdim2, int *sdim3,
				      int *ddim1, int *ddim2, int *ddim3, 
				      int *start1, int *start2, int *start3,
				      int *refine1, int *refine2,int *refine3);
#endif /* SPLINE */


int grid::PreparePotentialField(grid *ParentGrid)
{

  /* Return if this doesn't involve us. */

  if (MyProcessorNumber != ProcessorNumber && 
      MyProcessorNumber != ParentGrid->ProcessorNumber)
    return SUCCESS;

  if (CommunicationDirection == COMMUNICATION_SEND &&
      (MyProcessorNumber != ParentGrid->ProcessorNumber || 
       ProcessorNumber == ParentGrid->ProcessorNumber))
    return SUCCESS;

  if (CommunicationDirection == COMMUNICATION_RECEIVE &&
      MyProcessorNumber == ParentGrid->ProcessorNumber &&
      ProcessorNumber != ParentGrid->ProcessorNumber)
    return SUCCESS;

  if (ParentGrid == NULL)
    return FAIL;

  /* Return if potential field already exists. */

  //  if (PotentialField != NULL && NumberOfProcessors == 1)
  //    return SUCCESS;

  /* Allocate space for potential field. */

  int dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldDimension[dim];

  if (MyProcessorNumber == ProcessorNumber) {
    if (PotentialField != NULL) {
      //      fprintf(stderr, "What???\n"); 
      delete [] PotentialField;
    }
    PotentialField = new float[size];
  }
      
  /* Declarations. */

  int ParentOffset[MAX_DIMENSION], ParentStartIndex[MAX_DIMENSION],
      ParentTempDim[MAX_DIMENSION], Refinement[MAX_DIMENSION],
      ParentDim[MAX_DIMENSION];

  /* Compute refinement factors. */

  ParentGrid->ComputeRefinementFactors(this, Refinement);

  /* Set unused dims to zero. */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    ParentOffset[dim] = ParentStartIndex[dim] = 0;
    ParentTempDim[dim] = ParentDim[dim] = 1;
  }

  /* Compute the ParentOffset (in grid units) and ParentStartIndex and
     the region dim (in parent units). */

  for (dim = 0; dim < GridRank; dim++) {
    ParentOffset[dim] = nint((GravitatingMassFieldLeftEdge[dim] -
		  ParentGrid->GravitatingMassFieldLeftEdge[dim])/
			     GravitatingMassFieldCellSize);
    ParentStartIndex[dim] = ParentOffset[dim]/Refinement[dim] - 1;
    ParentTempDim[dim] = (ParentOffset[dim] + 
		       GravitatingMassFieldDimension[dim]-1)/Refinement[dim] -
			 ParentStartIndex[dim] + 3;
    ParentDim[dim] = ParentGrid->GravitatingMassFieldDimension[dim];
    if (ParentStartIndex[dim] < 0 || 
	ParentStartIndex[dim]+ParentTempDim[dim] > ParentDim[dim]) {
      fprintf(stderr, "ParentStartIndex[%d] = %d ParentTempDim = %d(%d).\n",
	      dim, ParentStartIndex[dim], ParentTempDim[dim], ParentDim[dim]);
      return FAIL;
    }
  }

  /* Copy data from other processor if needed (modify ParentDim and
     ParentStartIndex to reflect the fact that we are only coping part of
     the grid. */

  if (ProcessorNumber != ParentGrid->ProcessorNumber) {
    ParentGrid->CommunicationSendRegion(ParentGrid, ProcessorNumber, 
	       POTENTIAL_FIELD, NEW_ONLY, ParentStartIndex, ParentTempDim);
    for (dim = 0; dim < GridRank; dim++) {
      ParentOffset[dim] -= Refinement[dim]*ParentStartIndex[dim];
      ParentDim[dim] = ParentTempDim[dim];
    }
  }

  /* Return if this is not our concern. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Interpolate */

#ifdef SPLINE

  int size1 = GravitatingMassFieldDimension[0] *
              (GravitatingMassFieldDimension[1]/Refinement[1] + 2);
  int size2 = GravitatingMassFieldDimension[0] *
              GravitatingMassFieldDimension[1] *
              (GravitatingMassFieldDimension[2]/Refinement[2] + 2);
  float *Temp1 = new float[size1];
  float *Temp2 = new float[size2];

  FORTRAN_NAME(int_spline)(ParentGrid->PotentialField, 
			   PotentialField, &GridRank,
			   ParentDim, ParentDim+1, ParentDim+2,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   ParentOffset, ParentOffset+1, ParentOffset+2,
			   Refinement, Refinement+1, Refinement+2,
			   Temp1, Temp2);

  delete [] Temp1;
  delete [] Temp2;

#else /* SPLINE */

  FORTRAN_NAME(prolong)(ParentGrid->PotentialField, 
			    PotentialField, &GridRank,
			    ParentDim, ParentDim+1, ParentDim+2,
			    GravitatingMassFieldDimension,
			    GravitatingMassFieldDimension+1,
			    GravitatingMassFieldDimension+2,
			    ParentOffset, ParentOffset+1, ParentOffset+2,
			    Refinement, Refinement+1, Refinement+2);

#endif /* SPLINE */

  /* Clean up parent. */

  if (MyProcessorNumber != ParentGrid->ProcessorNumber) {
    delete ParentGrid->PotentialField;
    ParentGrid->PotentialField = NULL;
  }

  return SUCCESS;
}

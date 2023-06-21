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
/  GRID CLASS (DIFFERENCES THE POTENTIAL TO GET THE ACCELERATION FIELD)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
extern "C" void FORTRAN_NAME(comp_accel)(float *source, float *dest1, 
			    float *dest2, float *dest3, int *ndim, int *iflag,
                           int *sdim1, int *sdim2, int *sdim3, int *ddim1, 
			    int *ddim2, int *ddim3, 
                           int *start1, int *start2, int *start3, 
			    float *delx, float *dely, float *delz);
extern "C" void FORTRAN_NAME(smooth)(float *source1, float *source2,
				     float *source3, int *ndim, int *sdim1,
				     int *sdim2, int *sdim3, int *nsmooth);


int grid::ComputeAccelerationField(int DifferenceType, int level) 
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* declarations */

  int dim, size = 1, DiffFlag = 1, Offset[MAX_DIMENSION] = {0,0,0};
  float CellSize[MAX_DIMENSION] = {1,1,1};

  /* Set the number of zones to difference based on DifferenceType. */

  if (DifferenceType == PARTICLES || DifferenceType == ZEUS_GRIDS)
    DiffFlag = 0;

  /* Compute adot/a at time = t+1/2dt (time-centered). */

  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }

  /* Set cell size. */

  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    CellSize[dim] = a * GravitatingMassFieldCellSize;
    Offset[dim] = nint((CellLeftEdge[dim][0] - 
			GravitatingMassFieldLeftEdge[dim])/ CellWidth[dim][0]);
  }

  /* Loop over dimensions and difference acceleration. */

  for (dim = 0; dim < GridRank; dim++) {

    /* Allocate acceleration field. */

    if (AccelerationField[dim] != NULL) {
      //      fprintf(stderr, "AccelerationField allocated?\n");
      //      return FAIL;
      delete [] AccelerationField[dim];
    }

    AccelerationField[dim] = new float[size];

  }

  /* Difference potential. */

  DiffFlag = 1;
  FORTRAN_NAME(comp_accel)(PotentialField, AccelerationField[0], 
	      AccelerationField[1], AccelerationField[2], &GridRank, &DiffFlag,
	    GravitatingMassFieldDimension, GravitatingMassFieldDimension+1,
	      GravitatingMassFieldDimension+2,
	    GridDimension, GridDimension+1, GridDimension+2,
            Offset, Offset+1, Offset+2, CellSize, CellSize+1, CellSize+2);

  int i, j, k, index1, Skip = 1;
 
  if (DifferenceType == PARTICLES || DifferenceType == ZEUS_GRIDS )
    for (dim = 0; dim < GridRank; dim++) {
      index1 = size-1;
      for (k = GridDimension[2]-1; k >= 0; k--)
	for (j = GridDimension[1]-1; j >= 0; j--)
	  for (i = GridDimension[0]-1; i >= 0; i--, index1--)
	    if ((dim != 0 || i > 0) &&
		(dim != 1 || j > 0) && (dim != 2 || k > 0))
		AccelerationField[dim][index1] = 0.5*(
		 AccelerationField[dim][index1] + 
		 AccelerationField[dim][index1-Skip]);
      Skip *= GridDimension[dim];
    }

  /* Smooth if necessary. */

#define NO_SMOOTH_ACCEL
#ifdef SMOOTH_ACCEL
  int nsmooth = max(level - MaximumGravityRefinementLevel, 0);
  if (nsmooth > 0) {

    nsmooth = nint(0.5*POW(RefineBy, nsmooth-1));
    FORTRAN_NAME(smooth)(AccelerationField[0], AccelerationField[1], 
			 AccelerationField[2], &GridRank, 
			 GridDimension, GridDimension+1, GridDimension+2,
			 &nsmooth);
  }
#endif /* SMOOTH_ACCEL */

  return SUCCESS;
}


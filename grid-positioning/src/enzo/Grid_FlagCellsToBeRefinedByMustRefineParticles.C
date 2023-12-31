/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINE BY POSITION OF MUST-REFINE PARTS)
/
/  written by: Greg Bryan
/  date:       February, 1998
/  modified1:  May, 2009 by John Wise
/                Works with local particle storage in RebuildHierarchy,
/                which appropriately distributes memory usage.
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
extern "C" void PFORTRAN_NAME(cic_flag)(FLOAT *posx, FLOAT *posy,
			FLOAT *posz, int *ndim, int *npositions,
                        int *itype, int *ffield, FLOAT *leftedge,
                        int *dim1, int *dim2, int *dim3, FLOAT *cellsize,
			int *imatch);
 
 
int grid::FlagCellsToBeRefinedByMustRefineParticles()
{
  /* declarations */
 
  int i, dim, ParticleTypeToMatch, size = 1;
  FLOAT LeftEdge[MAX_DIMENSION], CellSize;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* Set Left edge of grid. */
 
  for (dim = 0; dim < GridRank; dim++) {
    LeftEdge[dim] = CellLeftEdge[dim][0];
    size *= GridDimension[dim];
  }
  ParticleTypeToMatch = PARTICLE_TYPE_MUST_REFINE;
  CellSize = float(CellWidth[0][0]);
 
  /* Loop over all the particles, using only particles marked as
     must-refine particles. */
 
  PFORTRAN_NAME(cic_flag)(
           ParticlePosition[0], ParticlePosition[1], ParticlePosition[2],
	   &GridRank, &NumberOfParticles, ParticleType, FlaggingField,
	   LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
	   &CellSize, &ParticleTypeToMatch);
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)
      NumberOfFlaggedCells++;

  return NumberOfFlaggedCells;
 
}

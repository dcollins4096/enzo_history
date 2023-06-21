/***********************************************************************
/
/  GRID CLASS (COMPUTE PARTICLE AND GRID ACCELERATIONS)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
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
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
 
/* function prototypes */
 
 
/* EvolveHierarchy function */
 
int grid::ComputeAccelerations(int level)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfParticles > 0) {
 
    /* Compute the acceleration field for particles from potential. */
 
    if (this->ComputeAccelerationField(PARTICLES, level) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeAccelerationField.\n");
      return FAIL;
    }
 
    /* Add any fixed (external) acceleration to field. */
/*
    if (this->AddExternalAcceleration() == FAIL) {
      fprintf(stderr, "Error in grid->AddFixedAcceleration.\n");
      return FAIL;
    }
*/
    /* Clear particle accelerations. */
 
    this->ClearParticleAccelerations();
 
    /* Move particles 1/2 step forward in preparation for interpolation. */
 
    this->UpdateParticlePosition(0.5*dtFixed);
 
    /* Interpolate the accelerations back to the grid and particles. */
 
    this->InterpolateParticlePositions(this, PARTICLES);
 
    /* Move particles 1/2 step backwards to return to original positions. */
 
    this->UpdateParticlePosition(-0.5*dtFixed);
 
    /* Clean up. */
 
    this->DeleteAccelerationField();
 
  } // end: if (NumberOfParticles > 0)
 
  /* Compute acceleration field for cells. */
 
  if (NumberOfBaryonFields > 0)
    if (this->ComputeAccelerationField(
          (HydroMethod == Zeus_Hydro) ? ZEUS_GRIDS : GRIDS, level) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeAccelerationField.\n");
      return FAIL;
    }
 
  return SUCCESS;
}

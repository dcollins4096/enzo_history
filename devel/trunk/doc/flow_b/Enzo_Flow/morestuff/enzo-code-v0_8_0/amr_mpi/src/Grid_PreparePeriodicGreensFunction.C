/***********************************************************************
/
/  GRID CLASS (PREPARES THE GREENS FUNCTION IN FOURIER SPACE)
/
/  written by: Greg Bryan
/  date:       November, 1998
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//  Compute derived quantites
//

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::PreparePeriodicGreensFunction(region *GreensRegion)
{

  /* Declarations. */

  int i, j, k, ii, jj, kk, dim, size;
  const float twopi = 2.0*3.14159;

  /* Error check. */

  if (GravityBoundaryType != TopGridPeriodic) {
    fprintf(stderr, "GravityBoundaryType %d not supported.\n", 
	    GravityBoundaryType);
    return FAIL;
  }

  /* -------------------------------------------------- */
  /* Generate a description of the initial data layout. */

  /* Set start index and dimension of region,
     and add 2 if we are on the right-hand side of the x-dimension,
     which is required for the real-to-complex FFT. */

  for (dim = 0, size = 1; dim < GridRank; dim++) {
    GreensRegion->StartIndex[dim] = nint((GridLeftEdge[dim] -
	       DomainLeftEdge[dim])/GravitatingMassFieldCellSize);
    GreensRegion->RegionDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    if (dim == 0 && GridRightEdge[dim]+0.5*CellWidth[dim][0] >
	            DomainRightEdge[dim]) 
      GreensRegion->RegionDim[dim] += 2;
    if (dim == 0)
      GreensRegion->RegionDim[dim] /= 2;
    size *= GreensRegion->RegionDim[dim];
  }

  /* set start index and region for unused dims. */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    GreensRegion->StartIndex[dim] = 0;
    GreensRegion->RegionDim[dim]  = 1;
  }

  /* Set processor number where data lives. */

  GreensRegion->Processor = ProcessorNumber;
  GreensRegion->Data      = NULL;

  /* If the data is on this processor then allocate and set it. */

  if (MyProcessorNumber == GreensRegion->Processor) {

    if (debug) printf("PrepareGreens: Start = %d %d %d  Dim = %d %d %d\n",
		      GreensRegion->StartIndex[0], GreensRegion->StartIndex[1],
		      GreensRegion->StartIndex[2],
		      GreensRegion->RegionDim[0], GreensRegion->RegionDim[1],
		      GreensRegion->RegionDim[2]);

    GreensRegion->Data = new float[size];

    /* Set Greens' function. */

    int n = 0, DomainDim[] = {1,1,1};
    for (dim = 0; dim < GridRank; dim++)
      DomainDim[dim] = nint((DomainRightEdge[dim] - DomainLeftEdge[dim])/
			    GravitatingMassFieldCellSize);
    float kx, ky, kz, ksqr, value;
    for (k = 0; k < GreensRegion->RegionDim[2]; k++) {
      kk = k + GreensRegion->StartIndex[2];
      kk = (kk > DomainDim[2]/2 + 1) ? kk - DomainDim[2] : kk;
      kz = kk * twopi;
      for (j = 0; j < GreensRegion->RegionDim[1]; j++) {
	jj = j + GreensRegion->StartIndex[1];
	jj = (jj > DomainDim[1]/2 + 1) ? jj - DomainDim[1] : jj;
	ky = jj * twopi;
	for (i = 0; i < GreensRegion->RegionDim[0]; i++, n++) {
	  ii = i + GreensRegion->StartIndex[0]/2;
	  kx = ii * twopi;
	  
	  ksqr = kx*kx + ky*ky + kz*kz;

#define INFLUENCE1

#ifdef INFLUENCE1
	  value = (ksqr == 0) ? 0.0 : -1.0/ksqr;
#endif /* INFLUENCE1 */

#ifdef INFLUENCE2
	  if (ksqr != 0)
	    value = -POW(0.5*GravitatingMassFieldCellSize, 2) /
                    (POW(sin(0.5*kx*GravitatingMassFieldCellSize), 2) +
                     POW(sin(0.5*ky*GravitatingMassFieldCellSize), 2) +
                     POW(sin(0.5*kz*GravitatingMassFieldCellSize), 2) );
	  else
	    value = 0;
#endif /* INFLUENCE2 */

	  GreensRegion->Data[n] = value;
	
	}
      }
    }

  } // end: if (MyProcessor)

  return SUCCESS;
}


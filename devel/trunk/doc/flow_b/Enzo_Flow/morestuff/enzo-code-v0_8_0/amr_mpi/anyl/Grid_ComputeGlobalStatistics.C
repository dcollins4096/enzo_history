/***********************************************************************
/
/  GRID CLASS (ADD THIS GRIDS INFO TO GLOBAL STATISTICS)
/
/  written by: Greg Bryan
/  date:       May, 2000
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "StarParticleData.h"

int grid::ComputeGlobalStatistics(float RegionLeftEdge[], 
				  float RegionRightEdge[],
	   int NumberOfTempBins, float *TempBinEdge, double *TempBinInfo,
	   int NumberOfDensityBins, float DensBinStart, float DensBinEnd,
			      double *TempDensInfo)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  
  int i, j, k, n, n1, dim, dmindex = 0, index, size = 1, DMOffset[3] = {0,0,0},
      start[] = {0,0,0}, stop[] = {0,0,0};

  /* Check To see if grid overlaps the projected field. */

  float dens, CellVolume = 1;
  for (dim = 0; dim < GridRank; dim++) {
    if (GridLeftEdge[dim] > RegionRightEdge[dim] ||
	GridRightEdge[dim] < RegionLeftEdge[dim])
      return SUCCESS;
    CellVolume *= CellWidth[dim][0];
    size *= GridDimension[dim];
  }

  /* Loop over grid, output density, temperature, etc. */

  if (NumberOfBaryonFields > 0) {

    /* Compute the temperature. */

    float *temperature = new float[size];
    if (temperature == NULL)
      printf("malloc failed\n");
    this->ComputeTemperatureField(temperature);

    /* Find the start and stop indicies for the region of interest. */

    for (dim = 0; dim < GridRank; dim++) {
      start[dim] = max(int((RegionLeftEdge[dim] - GridLeftEdge[dim]) /
			   CellWidth[dim][0]), 0) + GridStartIndex[dim];
      stop[dim]  = min(int((RegionRightEdge[dim] - GridLeftEdge[dim]) /
			   CellWidth[dim][0]),
		       GridEndIndex[dim] - GridStartIndex[dim]) + 
	GridStartIndex[dim];
    }

    /* Set dm field offset. */

    if (GravitatingMassFieldParticles != NULL)
      for (dim = 0; dim < GridRank; dim++)
	DMOffset[dim] = nint((CellLeftEdge[dim][0] - 
			      GravitatingMassFieldParticlesLeftEdge[dim])/
			     GravitatingMassFieldParticlesCellSize);

    /* Find fields. */

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);

    /* Loop over region of interest and sum stats. */

    for (k = start[2]; k <= stop[2]; k++)
      for (j = start[1]; j <= stop[1]; j++) {
	index = (j+k*GridDimension[1])*GridDimension[0] + start[0];
	if (GravitatingMassFieldParticles != NULL)
	  dmindex = ((j+DMOffset[1]) + 
		   (k+DMOffset[2])*GravitatingMassFieldParticlesDimension[1])*
	                           GravitatingMassFieldParticlesDimension[0] +
	    (start[0]+DMOffset[0]);
	for (i = start[0]; i <= stop[0]; i++, index++, dmindex++)

	  /* Only count the cell if it is not further refined. */

	  if (BaryonField[NumberOfBaryonFields][index] == 0) {

	    /* Sum up mass and volume in region and for cells within
	       each temperature range. */

	    dens = BaryonField[DensNum][index];

	    TempBinInfo[0] += dens*CellVolume;
	    TempBinInfo[1] += CellVolume;
	    TempBinInfo[2] += dens*dens*CellVolume; /* for clumping factor. */
                              
	    for (n = 0; n < NumberOfTempBins; n++)
	      if (temperature[index] > TempBinEdge[n] &&
		  temperature[index] <= TempBinEdge[n+1]) {

		/* Add mass in this cell into temperature bin. */

		TempBinInfo[n*3+3] += dens*CellVolume;
		TempBinInfo[n*3+4] += CellVolume;
		TempBinInfo[n*3+5] += dens*dens*CellVolume;

		/* Add mass in this cell into density/temperature bin. */

		n1 = int((log10(dens) - DensBinStart)*NumberOfDensityBins/
			 (DensBinEnd - DensBinStart));
		n1 = min(max(n1, 0), NumberOfDensityBins-1);
		TempDensInfo[n*NumberOfDensityBins + n1] += dens*CellVolume;
		TempDensInfo[(n+NumberOfTempBins)*NumberOfDensityBins 
			    + n1] += CellVolume;
                                      
	      }
		
#ifdef UNUSED
	    if (GravitatingMassFieldParticles != NULL)
	      fprintf(gridfptr, "%g ", GravitatingMassFieldParticles[dmindex]);
	    fprintf(gridfptr, "%g\n", CellWidth[0][0]);
#endif
	  }
      }

    delete [] temperature;

  }

  return SUCCESS;
}

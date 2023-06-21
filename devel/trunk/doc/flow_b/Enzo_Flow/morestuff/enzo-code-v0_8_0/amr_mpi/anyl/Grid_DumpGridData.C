/***********************************************************************
/
/  GRID CLASS (OUTPUT GRID/DM/STAR DATA STORED ON THIS GRID)
/
/  written by: Greg Bryan
/  date:       December, 1999
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

int grid::DumpGridData(FILE *gridfptr, FILE *dmfptr, FILE *starfptr)
{
  
  int i, j, k, dim, dmindex = 0, index, size, DMOffset[3] = {0,0,0};

  /* Loop over grid, output density, temperature, etc. */

  if (gridfptr != NULL && NumberOfBaryonFields > 0) {

    /* Compute the temperature. */

    size = GridDimension[0]*GridDimension[1]*GridDimension[2];
    float *temperature = new float[size];
    this->ComputeTemperatureField(temperature);

    /* Set dm field offset. */

    if (GravitatingMassFieldParticles != NULL)
      for (dim = 0; dim < GridRank; dim++)
	DMOffset[dim] = nint((CellLeftEdge[dim][0] - 
			      GravitatingMassFieldParticlesLeftEdge[dim]) /
			     GravitatingMassFieldParticlesCellSize);

    /* Find fields. */

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = (j+k*GridDimension[1])*GridDimension[0] + GridStartIndex[0];
	if (GravitatingMassFieldParticles != NULL)
	  dmindex = ((j+DMOffset[1]) + 
		   (k+DMOffset[2])*GravitatingMassFieldParticlesDimension[1])*
	                           GravitatingMassFieldParticlesDimension[0] +
	    (GridStartIndex[0]+DMOffset[0]);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0];
	     i++, index++, dmindex++)
	  if (BaryonField[NumberOfBaryonFields][index] == 0) {
	    fprintf(gridfptr, "%g %g ", BaryonField[DensNum][index],
		    temperature[index]);
	    if (GravitatingMassFieldParticles != NULL)
	      fprintf(gridfptr, "%g ", GravitatingMassFieldParticles[dmindex]);
	    fprintf(gridfptr, "%g\n", CellWidth[0][0]);
	  }
      }

    delete [] temperature;

  }

  /* Loop over particles, out dm and star particles */

  if (NumberOfParticleAttributes == 0) {
    if (dmfptr != NULL)
      for (i = 0; i < NumberOfParticles; i++)
	fprintf(dmfptr, "%f %f %f\n", ParticlePosition[0][i], 
		ParticlePosition[1][i], ParticlePosition[2][i]);
  } else {
    for (i = 0; i < NumberOfParticles; i++) 
      if (ParticleAttribute[0][i] < 0) {
	if (dmfptr != NULL)
	  fprintf(dmfptr, "%f %f %f\n", ParticlePosition[0][i], 
		  ParticlePosition[1][i], ParticlePosition[2][i]);
      } else {
	if (starfptr != NULL) {
	  fprintf(starfptr, "%f %f %f ", ParticlePosition[0][i], 
		  ParticlePosition[1][i], ParticlePosition[2][i]);
	  for (j = 0; j < NumberOfParticleAttributes; j++)
	    fprintf(starfptr, "%g ", ParticleAttribute[j][i]);
	  fprintf(starfptr, "\n");
	}
      }
  }

  return SUCCESS;
}

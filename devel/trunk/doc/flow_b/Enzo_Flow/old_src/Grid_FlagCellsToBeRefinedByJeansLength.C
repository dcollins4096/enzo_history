/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE JEAN'S CRITERION)
/
/  written by: Greg Bryan
/  date:       February, 1998
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
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

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);


int grid::FlagCellsToBeRefinedByJeansLength()
{
  /* declarations */

  int i, dim;

  /* error check */

  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  /* compute size */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Compute the temperature field. */

  float *temperature = new float[size];
  if (this->ComputeTemperatureField(temperature) == FAIL) {
    fprintf(stderr, "Error in grid->ComputeTemperature.\n");
    return -1;
  }
    
  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  /* Get density units. */

  float DensityUnits = 1, LengthUnits = 1, VelocityUnits, TimeUnits,
        TemperatureUnits;
  if (ComovingCoordinates)
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }

  /* Compute constant for Jean's length computation.
      l_j = sqrt((pi*k*T) / (G \rho m_p))  . */

  FLOAT JLSquared = (double(3.14159*1.38e-16/6.67e-8)/
                 (double(DensityUnits)*double(1.67e-24))) /
                (double(LengthUnits)*double(LengthUnits));

  /* This is the safety factor to decrease the Jean's length by. */

  JLSquared /= POW(RefineByJeansLengthSafetyFactor, 2);

/* printf("jl: JL, dx, t, d = %g %g %g %g\n", sqrt(JLSquared), CellWidth[0][3],
	 temperature[(3 + 3*GridDimension[1])*GridDimension[0]+3],
	 BaryonField[DensNum][(3 + 3*GridDimension[1])*GridDimension[0]+3]);*/

  /* Loop over grid. */

#ifdef UNUSED
  int j, k, index;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (j + k*GridDimension[1])*GridDimension[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	if ((CellWidth[0][i])*(CellWidth[0][i]) >
	    JLSquared*temperature[index+i]/BaryonField[DensNum][index+i])
	  FlaggingField[index+i]++;
    }
#endif /* UNUSED */

  FLOAT CellWidthSquared = CellWidth[0][0]*CellWidth[0][0];
  for (i = 0; i < size; i++)
    if (CellWidthSquared > JLSquared*temperature[i]/BaryonField[DensNum][i])
      FlaggingField[i]++;

  /* clean up */

  delete temperature;

  /* Count number of flagged Cells. */

  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  return NumberOfFlaggedCells;

}

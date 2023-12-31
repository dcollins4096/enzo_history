/***********************************************************************
/
/  GRID CLASS (APPLY A TIME-ACTION TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       September, 2000
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
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
#include "CosmologyParameters.h"
 
/* function prototypes */
 
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
 
 
 
int grid::ApplyTimeAction(int Type, float Parameter)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Get the cosmology units so we can convert temperature later. */
 
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
        VelocityUnits;
  if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		&TimeUnits, &VelocityUnits, InitialTimeInCodeUnits) == FAIL) {
    fprintf(stderr, "Error in CosmologyGetUnits.\n");
    return FAIL;
  }
 
  /* Determine the size of the grids. */
 
  int size = 1, i, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find the density, gas energy, velocities & total energy
     (where appropriate). */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
 
  /* ----------------------------------------------------------------
     1) - Add a specific amount of heat to gas  Parameter is temperature
          increase in degrees K (i.e. energy per unit particle). */
 
  if (Type == 1) {
 
    if (NumberOfBaryonFields > 0) {
 
      /* Set energy increase (per particle) in code units
	 (Note: assumes mu=0.59). */
 
      float EnergyIncrease = Parameter/TemperatureUnits/
	                     ((Gamma - 1.0)*0.59);
 
      /* Add to total energy field (if zeus -- thermal energy). */
 
      for (i = 0; i < size; i++)
	BaryonField[TENum][i] += EnergyIncrease;
 
      /* Add to thermal energy field, if required. */
 
      if (DualEnergyFormalism)
	for (i = 0; i < size; i++)
	  BaryonField[GENum][i] += EnergyIncrease;
 
    }
 
  } else {
 
    /* Type unknown. */
 
    fprintf(stderr, "TimeAction Type %"ISYM" unknown.\n", Type);
    return FAIL;
  }
 
  return SUCCESS;
}

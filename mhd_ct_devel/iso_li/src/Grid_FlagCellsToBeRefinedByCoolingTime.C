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
/  GRID CLASS (FLAG CELLS TO BE REFINED WHERE COOLING TIME < DX/CSOUND)
/
/  written by: Greg Bryan
/  date:       February, 2000
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

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


int grid::FlagCellsToBeRefinedByCoolingTime()
{
  /* declarations */

  int i, dim;

  /* error check */

  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */

  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);

  /* compute size */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Compute the cooling time. */

  float *cooling_time = new float[size];
  if (this->ComputeCoolingTime(cooling_time) == FAIL) {
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

  /* Loop over grid, looking for cells for which tcool/tsound < 1
     where tsound = dx/csound. */

  float Coef = Gamma*(Gamma - 1.0) / POW(a*CellWidth[0][0], 2);
  float gas_energy;

  int FieldIndex = TENum;
  if (DualEnergyFormalism)
    FieldIndex = GENum;
  if (HydroMethod == Zeus_Hydro || DualEnergyFormalism) {
    for (i = 0; i < size; i++)
      if (cooling_time[i]*cooling_time[i]*BaryonField[FieldIndex][i]*Coef
	  < 1.0)
	FlaggingField[i]++;
  } else {
    for (i = 0; i < size; i++) {
      gas_energy = BaryonField[TENum][i];
      for (dim = 0; dim < GridRank; dim++)
	gas_energy -= 0.5*BaryonField[Vel1Num+dim][i]*
	                  BaryonField[Vel1Num+dim][i];
      if (cooling_time[i]*cooling_time[i]*gas_energy*Coef < 1.0)
	FlaggingField[i]++;
    }
  }

  /* clean up */

  delete cooling_time;

  /* Count number of flagged Cells. */

  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  return NumberOfFlaggedCells;

}

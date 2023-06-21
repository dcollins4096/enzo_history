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
/  GRID CLASS (FIND PEAKS)
/
/  written by: Greg Bryan
/  date:       February, 1996
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

/* function prototypes */


int grid::FindDensityPeaks(int NumberOfPeaks, int &NumberOfPeaksFound,
			  float *PeakPosition[MAX_DIMENSION], float *PeakValue,
			  float PeakSeparation, float MinDensity,
			  int PeakField)
{
  
  int i, j, k, i0, j0, k0, i1, j1, k1, index, index1, index2, n, n1, dim;
  float Position[MAX_DIMENSION], delx, dely, delz, DomainWidth[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (NumberOfBaryonFields > 0)
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum);

  /* Compute temperature */

  float *temperature, CurrentPeakValue;
  if (PeakField == 3) {
    int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
    temperature = new float[size];
    if (this->ComputeTemperatureField(temperature) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeTemperatureField\n");
      return FAIL;
    }
  }

  /* Loop over grid. */

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	/* Check to see if it is above the minimum density requirement. */

	if (NumberOfBaryonFields > 0)
	  if (BaryonField[DensNum][index] < MinDensity)
	    goto NotAPeak;

	/* Check to see if it is a peak. */

	switch (PeakField) {

	case 0:
	  for (k1 = -1; k1 <= 1; k1++)
	    for (j1 = -1; j1 <= 1; j1++) {
	      index1 = ((k+k1)*GridDimension[1] + j+j1)*GridDimension[0] + 
		         i-1;
	      for (i1 = -1; i1 <= 1; i1++, index1++)
		if (BaryonField[DensNum][index] < 
		    BaryonField[DensNum][index1])
		  goto NotAPeak;
	    }
	  CurrentPeakValue = BaryonField[DensNum][index];
	  break;

	case 1:
	  i0 = nint((CellLeftEdge[0][0]-GravitatingMassFieldLeftEdge[0])/
                      GravitatingMassFieldCellSize);
	  j0 = nint((CellLeftEdge[1][0]-GravitatingMassFieldLeftEdge[1])/
                      GravitatingMassFieldCellSize);
	  k0 = nint((CellLeftEdge[2][0]-GravitatingMassFieldLeftEdge[2])/
                      GravitatingMassFieldCellSize);
	  index2 = ((k+k0)*GravitatingMassFieldDimension[1] + j+j0)*
	              GravitatingMassFieldDimension[0] + i+i0;
	  for (k1 = -1; k1 <= 1; k1++)
	    for (j1 = -1; j1 <= 1; j1++) {
	      index1 = ((k+k0+k1)*GravitatingMassFieldDimension[1] +j+j0+j1)*
		  GravitatingMassFieldDimension[0] + i+i0-1;
	      for (i1 = -1; i1 <= 1; i1++, index1++)
		if (GravitatingMassField[index2] < 
		    GravitatingMassField[index1])
		  goto NotAPeak;
	    }
	  CurrentPeakValue = GravitatingMassField[index2];
	  break;

	case 2:
	  i0 = nint((CellLeftEdge[0][0]-GravitatingMassFieldLeftEdge[0])/
                      GravitatingMassFieldCellSize);
	  j0 = nint((CellLeftEdge[1][0]-GravitatingMassFieldLeftEdge[1])/
                      GravitatingMassFieldCellSize);
	  k0 = nint((CellLeftEdge[2][0]-GravitatingMassFieldLeftEdge[2])/
                      GravitatingMassFieldCellSize);
	  index2 = ((k+k0)*GravitatingMassFieldDimension[1] + j+j0)*
	            GravitatingMassFieldDimension[0] + i+i0;
	  for (k1 = -1; k1 <= 1; k1++)
	    for (j1 = -1; j1 <= 1; j1++) {
	      index1 = ((k+k0+k1)*GravitatingMassFieldDimension[1] +j+j0+j1)*
		  GravitatingMassFieldDimension[0] + i+i0-1;
	      for (i1 = -1; i1 <= 1; i1++, index1++)
		if (AccelerationField[GridRank][index2] >
		      AccelerationField[GridRank][index1])
		  goto NotAPeak;
	    }
	  CurrentPeakValue = -GravitatingMassField[index2];
	  break;

	case 3:
	  for (k1 = -1; k1 <= 1; k1++)
	    for (j1 = -1; j1 <= 1; j1++) {
	      index1 = ((k+k1)*GridDimension[1] + j+j1)*GridDimension[0] + i-1;
	      for (i1 = -1; i1 <= 1; i1++, index1++)
		if (BaryonField[DensNum][index]*BaryonField[DensNum][index]*
		      temperature[index] < 
		    BaryonField[DensNum][index1]*BaryonField[DensNum][index1]*
		      temperature[index1])
		  goto NotAPeak;
	    }
	  CurrentPeakValue = BaryonField[DensNum][index]*
	    BaryonField[DensNum][index]*temperature[index];
	  break;

	}

	/* Determine Position. */

	Position[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 0)
	  Position[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 1)
	  Position[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	/* Check to see if we're too close to another peak. */
	  
	for (n = 0; n < NumberOfPeaksFound; n++) {

	  /* Compute radius */

	  if (GridRank > 1) {
	    delz = Position[2] - PeakPosition[2][n];
	    delz = min(delz, DomainWidth[2]-delz);
	  }
	  else
	    delz = 0;
	  if (GridRank > 0) {
	    dely = Position[1] - PeakPosition[1][n];
	    dely = min(dely, DomainWidth[1]-dely);
	  }
	  else
	    dely = 0;
	  delx = Position[0] - PeakPosition[0][n];
	  delx = min(delx, DomainWidth[0]-delx);

	  /* Check if we are overlapping this other peak. */

	  if (delx*delx+dely*dely+delz*delz < PeakSeparation*PeakSeparation)

	    /* If yes, throw out peak with lower density. */

	    if (CurrentPeakValue < PeakValue[n])
	      goto NotAPeak;
	    else {
	      NumberOfPeaksFound--;
	      for (n1 = n; n1 < NumberOfPeaksFound; n1++) {
		for (dim = 0; dim < GridRank; dim++)
		  PeakPosition[dim][n1] = PeakPosition[dim][n1+1];
		PeakValue[n1] = PeakValue[n1+1];
	      }
	    }
	    
	} // end loop over other peaks
	
	/* Rearrange the peak table. */
	/* Find position in list of peaks. */

	for (n = 0; n < NumberOfPeaksFound; n++)
	  if (CurrentPeakValue > PeakValue[n])
	    break;

	/* Move peak n and all below one down list. */
	
	for (n1 = min(NumberOfPeaksFound, NumberOfPeaks-1); n1 > n; n1--) {
	  for (dim = 0; dim < GridRank; dim++)
	    PeakPosition[dim][n1] = PeakPosition[dim][n1-1];
	  PeakValue[n1] = PeakValue[n1-1];
	}

	/* Overwrite peak position n. */

	if (n != NumberOfPeaks) {
	  if (NumberOfPeaksFound != NumberOfPeaks)
	    NumberOfPeaksFound++;
	  for (dim = 0; dim < GridRank; dim++)
	    PeakPosition[dim][n] = Position[dim];
	  PeakValue[n] = CurrentPeakValue;
	}

      NotAPeak:
	;
	
      } // end: loop over i
    } // end: loop over j

  if (PeakField == 3)
    delete temperature;

  return SUCCESS;
}

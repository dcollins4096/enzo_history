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
/  GRID CLASS (RESTORE CONSISTENCY BETWEEN TOTAL AND INTERNAL ENERGY)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/    If using the dual energy formalism, this routine will restore
/     consistency between the total and internal (gas) energy fields.
/     This is done to either the entire field or just the boundary zones.
/     This can result in the supression of small shocks, so the routine
/     should be used with caution (i.e. just after interpolation which
/     generates errors of order a few percent anyway).
/     
/
/  RETURNS: FAIL or SUCCESS
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


/* Interpolate Boundary From Parent function */

int grid::RestoreEnergyConsistency(int Region)
{

  int i, j, k, n;

  /* Error check */

  if (Region != ENTIRE_REGION && Region != ONLY_BOUNDARY) {
    fprintf(stderr, "Region type %d unknown.\n", Region);
    return FAIL;
  }

  /* If there is no work, we're done. */

  if (NumberOfBaryonFields <= 0 || DualEnergyFormalism == 0)
    return SUCCESS;

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  /* a) Correct the entire field. */

  if (Region == ENTIRE_REGION) {

    /* Compute size of field. */

    int size = 1;
    for (int dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

    /* Perform correction:  E = e + 1/2*v^2. */

    for (i = 0; i < size; i++)
      BaryonField[TENum][i] = BaryonField[GENum][i] +
	         0.5*(BaryonField[Vel1Num][i])*(BaryonField[Vel1Num][i]);

    if (GridRank > 1)
      for (i = 0; i < size; i++)
	BaryonField[TENum][i] += 
	         0.5*(BaryonField[Vel2Num][i])*(BaryonField[Vel2Num][i]);

    if (GridRank > 2)
      for (i = 0; i < size; i++)
	BaryonField[TENum][i] += 
	         0.5*(BaryonField[Vel3Num][i])*(BaryonField[Vel3Num][i]);

#ifdef HAOXU
     if(MHD_Used == TRUE)
            for (i = 0; i < size; i++){
                BaryonField[TENum][i] *=BaryonField[DensNum][i];
                BaryonField[TENum][i] +=0.5*(pow(CenteredB[0][i],2)+pow(CenteredB[1][i],2)+pow(CenteredB[2][i],2));
            }
#endif

  } // end: Region == ENTIRE_FIELD

  /* b) Correct just the boundary zones. */

  if (Region == ONLY_BOUNDARY) {

    for (k = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++)
	  if (i < GridStartIndex[0] || i > GridEndIndex[0] ||
	      j < GridStartIndex[1] || j > GridEndIndex[1] ||
	      k < GridStartIndex[2] || k > GridEndIndex[2]   ) {

	    n = k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i;

	    /* Perform correction:  E = e + 1/2*v^2. */

	    BaryonField[TENum][n] = BaryonField[GENum][n] +
	         0.5*(BaryonField[Vel1Num][n])*(BaryonField[Vel1Num][n]);

	    if (GridRank > 1)
	      BaryonField[TENum][n] += 
	         0.5*(BaryonField[Vel2Num][n])*(BaryonField[Vel2Num][n]);

	    if (GridRank > 2)
	      BaryonField[TENum][n] += 
	         0.5*(BaryonField[Vel3Num][n])*(BaryonField[Vel3Num][n]);

#ifdef HAOXU
     if(MHD_Used == TRUE){
                BaryonField[TENum][n] *=BaryonField[DensNum][n];
                BaryonField[TENum][n] +=0.5*(pow(CenteredB[0][n],2)+pow(CenteredB[1][n],2)+pow(CenteredB[2][n],2));                                                                                                                                                                                                                  
      }
#endif
	  }

  } // end: Region == BOUNDARY_ONLY
  
  return SUCCESS;

}

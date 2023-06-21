/*****************************************************************************
 *                                                                           *
 * Copyright 2007 Rick Wagner                                                *
 * Copyright 2007 Laboratory for Computational Astrophysics                  *
 * Copyright 2007 Board of Trustees of the University of Illinois            *
 * Copyright 2007 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRID CLASS (Bin PDFs for various values)
/
/  written by: Rick Wagner
/  date:       October, 2007
/  modified1:
/
/  PURPOSE: Simple example.
/
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
#include "ealInt.h"
#include "ealFloat.h"
#include "ProbabilityDistributionFunction.h"

/* function prototypes */
void grid::BinPDF( void *dens_func,
		   void *en_func ){

  ProbabilityDistributionFunction *dens_func_p =
    (ProbabilityDistributionFunction *)dens_func;

  ProbabilityDistributionFunction *en_func_p =
    (ProbabilityDistributionFunction *)en_func;

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber )
    return;

  /* declarations */
    
  int i, j, k, bin, dim, index;
  float cellvol = 1.0, dens, en;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    cellvol*= CellWidth[dim][0];
  }

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	if(!FlaggingField[index + i]){
	  dens = BaryonField[DensNum][index+i];

	  if(dens < dens_func_p->pdf_min){
	    dens_func_p->under->Array[0] += cellvol;
	  }else if(dens > dens_func_p->pdf_max){
	    dens_func_p->over->Array[0] += cellvol;
	  }else{
	    bin = (int)floor(log10((dens)/(dens_func_p->pdf_min))/(dens_func_p->log_bin_width));
	    dens_func_p->pdf->Array[bin] += cellvol;
	  }

	  en = BaryonField[TENum][index+i];
	  en = en - 0.5*(BaryonField[Vel1Num][index+i]*BaryonField[Vel1Num][index+i]+
		    BaryonField[Vel2Num][index+i]*BaryonField[Vel2Num][index+i]+
		    BaryonField[Vel3Num][index+i]*BaryonField[Vel3Num][index+i]);

	  if(en < en_func_p->pdf_min){
	    en_func_p->under->Array[0] += cellvol;
	  }else if(en > en_func_p->pdf_max){
	    en_func_p->over->Array[0] += cellvol;
	  }else{
	    bin = (int)floor(log10((en)/(en_func_p->pdf_min))/(en_func_p->log_bin_width));
	    en_func_p->pdf->Array[bin] += cellvol;
	  }
	}
    }

  return;
}

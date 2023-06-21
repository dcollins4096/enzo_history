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
/  GRID CLASS (Set the Miminum Speed for ES Switch )
/
/  written by: Hao Xu
/  date:       March, 2006
/  modified1:
/
/  PURPOSE:
/
********************************************************************/
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include <stdlib.h>

#ifdef HAOXU
float grid::SetESSpeed()
{

   /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return huge_number;


  /* Compute grid size. */

  int size = 1, field, i;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  float MHDTVDSpeed,temp;

 //
  //Figure out which BaryonField contains which hydro variable.
  //
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
  }

  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
    
  
    MHDTVDSpeed = huge_number;
   
    for (i = 0; i < size; i++) {
     temp=  BaryonField[GENum][i]*Gamma*(Gamma-1);
     for(field=0;field<3;field++)
      temp+= CenteredB[field][i]*CenteredB[field][i]/BaryonField[DensNum][i];
     MHDTVDSpeed = min(sqrt(max(0.0,temp)),MHDTVDSpeed); 
   }
   
     MHDTVDSpeed*=100;
   
     return MHDTVDSpeed;

}
#endif //HAOXU

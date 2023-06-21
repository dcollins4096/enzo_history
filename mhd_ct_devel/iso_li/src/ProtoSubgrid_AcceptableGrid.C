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
/  PROTOSUBGRID CLASS (CHECK IF A SUBGRID NEEDS TO BE SUBDIVIDED)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
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

#define MINIMUM_SIDE_LENGTH 5
//64^3 = 262144
#define MAXIMUM_SIZE 262144
//#define MAXIMUM_SIZE 288
int ProtoSubgrid::AcceptableSubgrid()
{
  //fprintf(stderr,"kludge: max_size = %d\n",MAXIMUM_SIZE);
  /* If NumberFlagged hasn't been computed yet, then compute it. */

  if (NumberFlagged == INT_UNDEFINED) {
    this->ComputeSignature(0);
    NumberFlagged = 0;
    for (int i = 0; i < GridDimension[0]; i++)
      NumberFlagged += Signature[0][i];
  }

  /* Compute size and efficiency. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  float efficiency = float(NumberFlagged)/float(size);

  /* Subgrid is acceptable if it is efficient enough or small enough,
     but if we are using multiple processors, then make sure it is
     not too big (for good load balancing). */

  if (size <= POW(float(MINIMUM_SIDE_LENGTH), GridRank))
    return TRUE;


  //if (size > MAXIMUM_SIZE && NumberOfProcessors > 1)
    if (size > MAXIMUM_SIZE )
    return FALSE;
   
  if (efficiency > MinimumEfficiency)
    return TRUE;

  /* Not acceptable yet -- keep refining. */
  //return TRUE;
  return FALSE;
}


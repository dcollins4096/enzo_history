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
/  GRID CLASS (CHECKS FOR NANS)
/
/  written by: Greg Bryan
/  date:       1999
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"

// Turn the first on to check the gas and dark matter data for nan's.
// Turn the second on just to print the message (i.e. which routine
//     has called it.

#define DEBUG_CHECK_OFF
#define TRACE_OFF

/* function prototypes */

int  ReportMemoryUsage(char *header = NULL);

int grid::DebugCheck(char *message)
{

#ifdef TRACE_ON

  if (ProcessorNumber == MyProcessorNumber)
    fprintf(stderr, "P(%d): %s\n", MyProcessorNumber, message);
    //    ReportMemoryUsage(message);

#endif /* TRACE_ON */

#ifdef DEBUG_CHECK_ON

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Set this to zero but don't the compiler so it doesn't optimize everything
     away. */

  int ThisIsZero = GridStartIndex[0] - DEFAULT_GHOST_ZONES, size = 1, 
      dim, k1, k2;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (k1 = 0; k1 < NumberOfBaryonFields; k1++)
    for (k2 = 0; k2 < size; k2++)
      if (BaryonField[k1][k2+ThisIsZero] != BaryonField[k1][k2]) {
	fprintf(stderr, "DebugCheck[%s](Proc %d): gas %d %d %g\n", message, 
		MyProcessorNumber, k1, k2, BaryonField[k1][k2]);
	exit(EXIT_FAILURE);
      }

  for (k1 = 0; k1 < GridRank; k1++)
    for (k2 = 0; k2 < NumberOfParticles; k2++)
      if (ParticlePosition[k1][k2] != ParticlePosition[k1][k2+ThisIsZero] ||
	  ParticleVelocity[k1][k2] != ParticleVelocity[k1][k2+ThisIsZero]  ) {
	fprintf(stderr, "DebugCheck[%s](Proc %d): dm %d (%d/%d) %g %g\n", 
		message, MyProcessorNumber, k1, k2, NumberOfParticles,
		ParticlePosition[k1][k2], ParticleVelocity[k1][k2]);
	exit(EXIT_FAILURE);
      }

#if 0		
  if (NumberOfBaryonFields > 0)
    for (k1 = 0; k1 < 2+DualEnergyFormalism; k1++)
      for (k2 = 0; k2 < size; k2++)
	if (BaryonField[k1][k2+ThisIsZero] <= 0) {
	  fprintf(stderr, "DebugCheck[%s](Proc %d): <0 %d %d %g\n", message, 
		  MyProcessorNumber, k1, k2, BaryonField[k1][k2]);
	  exit(EXIT_FAILURE);
	}
#endif

#endif /* DEBUG_CHECK_ON */
		
  return SUCCESS;
}

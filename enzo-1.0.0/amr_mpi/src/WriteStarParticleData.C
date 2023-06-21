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
/  WRITES A STAR PARTICLE DATA
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include "hdf4.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"


int WriteStarParticleData(FILE *fptr)
{

  if (StarParticleCreation == FALSE)
    return SUCCESS;

  /* Write out number data. */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "NumberOfStarParticles      = %d\n", NumberOfStarParticles);

  return SUCCESS;
}

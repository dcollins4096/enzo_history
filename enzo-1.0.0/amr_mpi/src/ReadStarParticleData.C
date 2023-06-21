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
/  READS THE STAR PARTICLE DATA
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


int ReadStarParticleData(FILE *fptr)
{

  if (StarParticleCreation == FALSE)
    return SUCCESS;

  /* read in number data. */

  if (fscanf(fptr, "NumberOfStarParticles = %d\n", 
	     &NumberOfStarParticles) != 1) {
    fprintf(stderr, "Error reading NumberOfStarParticles.\n");
    return FAIL;
  }

  return SUCCESS;
}

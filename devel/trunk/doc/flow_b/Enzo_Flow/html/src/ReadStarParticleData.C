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

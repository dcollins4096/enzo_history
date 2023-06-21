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
/  WRITES THE RADIATION FIELD DATA
/
/  written by: Greg Bryan
/  date:       October, 1999
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


int WriteRadiationData(FILE *fptr)
{
  int i;

  /* write scalar data. */

  fprintf(fptr, "TimeFieldLastUpdated = %"GOUTSYM"\n", 
	  RadiationData.TimeFieldLastUpdated);

  /* write field. */

  for (i = 0; i < RadiationData.NumberOfFrequencyBins; i++)
    fprintf(fptr, "%g %g %g %g\n",
	       RadiationData.Spectrum[0][i], RadiationData.Spectrum[1][i],
	       RadiationData.Spectrum[2][i], RadiationData.Spectrum[3][i]);

  return SUCCESS;
}

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
/  READS COSMOLOGY PARAMETERS FROM INPUT FILE
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "CosmologyParameters.h"

int CosmologyComputeTimeFromRedshift(float Redshift, float *TimeCodeUnits);

int CosmologyReadParameters(FILE *fptr)
{

  char line[MAX_LINE_LENGTH], *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  /* Set defaults. */

  HubbleConstantNow    = 0.5;
  OmegaMatterNow       = 1.0;
  OmegaLambdaNow       = 0.0;
  OmegaWDMNow          = 0.0;
  OmegaHDMNow          = 0.0;
  OmegaBaryonNow       = 0.06;
  ComovingBoxSize      = 64.0;
  InitialRedshift      = 20.0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    int ret = 0;

    /* read parameters */
    
    ret += sscanf(line, "CosmologyHubbleConstantNow = %f", &HubbleConstantNow);
    ret += sscanf(line, "CosmologyOmegaMatterNow = %f", &OmegaMatterNow);
    ret += sscanf(line, "CosmologyOmegaLambdaNow = %f", &OmegaLambdaNow);
    ret += sscanf(line, "CosmologyOmegaWDMNow = %f", &OmegaWDMNow);
    ret += sscanf(line, "CosmologyOmegaHDMNow = %f", &OmegaHDMNow);
    ret += sscanf(line, "CosmologyOmegaBaryonNow = %f", &OmegaBaryonNow);
    ret += sscanf(line, "CosmologyComovingBoxSize = %f", &ComovingBoxSize);
    ret += sscanf(line, "CosmologyInitialRedshift = %f", &InitialRedshift);

    /* If the dummy char space was used, then make another. */

    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#' && 
	strstr(line, "Cosmology"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }
  
  return SUCCESS;
}

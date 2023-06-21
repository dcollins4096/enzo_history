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
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

int CosmologyComputeTimeFromRedshift(FLOAT Redshift, FLOAT *TimeCodeUnits);

int CosmologyReadParameters(FILE *fptr, FLOAT *StopTime, FLOAT *InitTime)
{

  int i, OutputNumber;
  FLOAT FinalRedshift, CurrentRedshift;
  char line[MAX_LINE_LENGTH], *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  /* Set defaults. */

  HubbleConstantNow    = 0.5;
  OmegaMatterNow       = 1;
  OmegaLambdaNow       = 0;
  ComovingBoxSize      = 64;
  MaxExpansionRate     = 0.01;
  InitialRedshift      = 20;
  FinalRedshift        = 0;

#ifdef HAOXU
  MaxExpansion         = 1.0e20;
  MHD_Equation         = 2;
#endif

#ifndef HAOXU
  MHD_Equation         = 2;
#endif 

  for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++) {
    CosmologyOutputRedshift[i]     = -1;  // Never!!
    CosmologyOutputRedshiftName[i] = NULL;
  }

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    int ret = 0;

    /* read parameters */
    
    ret += sscanf(line, "CosmologyHubbleConstantNow = %"FSYM, 
		  &HubbleConstantNow);
    ret += sscanf(line, "CosmologyOmegaMatterNow = %"FSYM, &OmegaMatterNow);
    ret += sscanf(line, "CosmologyOmegaLambdaNow = %"FSYM, &OmegaLambdaNow);
    ret += sscanf(line, "CosmologyComovingBoxSize = %"FSYM, &ComovingBoxSize);
    ret += sscanf(line, "CosmologyMaxExpansionRate = %"FSYM, 
		  &MaxExpansionRate);
    ret += sscanf(line, "CosmologyInitialRedshift = %"PSYM, &InitialRedshift);
    ret += sscanf(line, "CosmologyFinalRedshift = %"PSYM, &FinalRedshift);
    ret += sscanf(line, "CosmologyCurrentRedshift = %"PSYM, &CurrentRedshift);

#ifdef HAOXU
    ret += sscanf(line, "CosmologyMaxExpansion = %"FSYM,&MaxExpansion);
    ret += sscanf(line, "MHD_Equation = %d",&MHD_Equation);
#endif

    if (sscanf(line, "CosmologyOutputRedshift[%d] =", &OutputNumber) == 1)
      ret += sscanf(line, "CosmologyOutputRedshift[%d] = %"PSYM,
		    &OutputNumber, &CosmologyOutputRedshift[OutputNumber]);
    if (sscanf(line, "CosmologyOutputRedshiftName[%d] = %s", 
	       &OutputNumber, dummy) == 2)
      CosmologyOutputRedshiftName[OutputNumber] = dummy;

    /* If the dummy char space was used, then make another. */

    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#' && 
	strstr(line, "Cosmology") && !strstr(line, "CosmologySimulation"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }

  /* Initialize by finding the time at the initial redshift. */

  if (CosmologyComputeTimeFromRedshift(InitialRedshift, 
				       &InitialTimeInCodeUnits) == FAIL) {
    fprintf(stderr, "Error in ComputeTimeFromRedshift.\n");
    return FAIL;
  }
  if (*InitTime == 0.0)
    *InitTime = InitialTimeInCodeUnits;

  /* Now find the time at the end of the simulation. */

  if (CosmologyComputeTimeFromRedshift(FinalRedshift, StopTime) == FAIL) {
    fprintf(stderr, "Error in ComputeTimeFromRedshift.\n");
    return FAIL;
  }  

  /* Convert the output redshifts into time, for later convenience. */

  for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
    if (CosmologyOutputRedshift[i] != -1)
      CosmologyComputeTimeFromRedshift(CosmologyOutputRedshift[i],
				       &CosmologyOutputRedshiftTime[i]);

  /* Convert the time action redshift into time. */

  for (i = 0; i < MAX_TIME_ACTIONS; i++)
    if (TimeActionRedshift[i] != -1)
      CosmologyComputeTimeFromRedshift(TimeActionRedshift[i],
				       &TimeActionTime[i]);
  
  return SUCCESS;
}

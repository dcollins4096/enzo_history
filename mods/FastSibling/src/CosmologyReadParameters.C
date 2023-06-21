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

  for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++) {
    CosmologyOutputRedshift[i]     = -1;  // Never!!
    CosmologyOutputRedshiftName[i] = NULL;
  }

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    int ret = 0;

    /* read parameters */
    
    ret += sscanf(line, "CosmologyHubbleConstantNow = %f", 
		  &HubbleConstantNow);
    ret += sscanf(line, "CosmologyOmegaMatterNow = %f", &OmegaMatterNow);
    ret += sscanf(line, "CosmologyOmegaLambdaNow = %f", &OmegaLambdaNow);
    ret += sscanf(line, "CosmologyComovingBoxSize = %f", &ComovingBoxSize);
    ret += sscanf(line, "CosmologyMaxExpansionRate = %f", 
		  &MaxExpansionRate);
    ret += sscanf(line, "CosmologyInitialRedshift = %"FSYM, &InitialRedshift);
    ret += sscanf(line, "CosmologyFinalRedshift = %"FSYM, &FinalRedshift);
    ret += sscanf(line, "CosmologyCurrentRedshift = %"FSYM, &CurrentRedshift);

    if (sscanf(line, "CosmologyOutputRedshift[%d] =", &OutputNumber) == 1) {
      if (OutputNumber > MAX_NUMBER_OF_OUTPUT_REDSHIFTS-1) {
	fprintf(stderr, "OutputRedshift number(%d) > MAX\n", OutputNumber);
	return FAIL;
      }
      ret += sscanf(line, "CosmologyOutputRedshift[%d] = %"FSYM,
		    &OutputNumber, &CosmologyOutputRedshift[OutputNumber]);
    }
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

/***********************************************************************
/
/  WRITES COSMOLOGY PARAMETERS TO AN OUTPUT FILE
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

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


int CosmologyWriteParameters(FILE *fptr, FLOAT StopTime, FLOAT CurrentTime)
{

  /* Compute the final redshift from StopTime. */

  FLOAT a, dadt, FinalRedshift, CurrentRedshift;
  if (CosmologyComputeExpansionFactor(StopTime, &a, &dadt) == FAIL) {
    fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
    return FAIL;
  }
  FinalRedshift = (1 + InitialRedshift)/a - 1;

  /* Compute the current redshift (for information only). */

  CosmologyComputeExpansionFactor(CurrentTime, &a, &dadt);
  CurrentRedshift = (1 + InitialRedshift)/a - 1;

  /* write output to file */

  for (int i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
    if (CosmologyOutputRedshift[i] != -1) {
      fprintf(fptr, "CosmologyOutputRedshift[%d] = %"GOUTSYM"\n", i, 
	      CosmologyOutputRedshift[i]);
      if (CosmologyOutputRedshiftName[i] != NULL)
	fprintf(fptr, "CosmologyOutputRedshiftName[%d] = %s\n", i,
		CosmologyOutputRedshiftName[i]);
    }
    
  fprintf(fptr, "CosmologyHubbleConstantNow = %g\n", HubbleConstantNow);
  fprintf(fptr, "CosmologyOmegaMatterNow    = %g\n", OmegaMatterNow);
  fprintf(fptr, "CosmologyOmegaLambdaNow    = %g\n", OmegaLambdaNow);
  fprintf(fptr, "CosmologyComovingBoxSize   = %g\n", ComovingBoxSize);
  fprintf(fptr, "CosmologyMaxExpansionRate  = %g\n", MaxExpansionRate);
  fprintf(fptr, "CosmologyInitialRedshift   = %"GOUTSYM"\n", InitialRedshift);
  fprintf(fptr, "CosmologyFinalRedshift     = %"GOUTSYM"\n", FinalRedshift);
  fprintf(fptr, "CosmologyCurrentRedshift   = %"GOUTSYM"\n\n", CurrentRedshift);
  
  return SUCCESS;
}

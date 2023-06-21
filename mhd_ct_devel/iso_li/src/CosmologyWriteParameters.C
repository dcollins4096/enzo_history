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
#ifdef HAOXU
  fprintf(fptr, "CosmologyMaxExpansion      = %g\n", MaxExpansion);
  fprintf(fptr, "MHD_Equation      = %d\n", MHD_Equation);
#endif
  fprintf(fptr, "CosmologyInitialRedshift   = %"GOUTSYM"\n", InitialRedshift);
  fprintf(fptr, "CosmologyFinalRedshift     = %"GOUTSYM"\n", FinalRedshift);
  fprintf(fptr, "CosmologyCurrentRedshift   = %"GOUTSYM"\n\n", CurrentRedshift);
  
  return SUCCESS;
}

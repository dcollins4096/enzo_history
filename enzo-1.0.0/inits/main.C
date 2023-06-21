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
/  INITS MAIN CODE
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/    This code generates a gaussian (random phase) realization of a
/    large, discrete field, given the power spectrum of perturbations.
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#ifdef USE_JBMEM
#include "jbMem.h"
#endif
#define DEFINE_STORAGE
#include "global_data.h"
#include "CosmologyParameters.h"
#include "PowerSpectrumParameters.h"
#undef DEFINE_STORAGE
#include "Parameters.h"

/* function prototypes */

int InterpretCommandLine(int argc, char *argv[], char *myname,
			 char *ParameterFile[], char *SubGridParameterFile[]);
int SetParameterDefaults(parmstruct *Parameters);
int ReadParameterFile(FILE *fptr, parmstruct *Parameters);
int InitializePowerSpectrum();
int MakePowerSpectrumLookUpTable();
int GenerateRealization(parmstruct *Parameters, parmstruct *SubGridParameters);
int CosmologyReadParameters(FILE *fptr);
int ReadPowerSpectrumParameters(FILE *fptr);


int main(int argc, char *argv[])
{

  /* Initialize */

  debug                        = FALSE;
  char *myname                 = argv[0];
  parmstruct Parameters, *SubGridParameters = NULL;

  char *ParameterFile = NULL, *SubGridParameterFile = NULL;
  FILE *fptr;

  /* Interpret command-line arguments. */

  InterpretCommandLine(argc, argv, myname, &ParameterFile, 
		       &SubGridParameterFile);

  /* Set Parameter defaults. */

  SetParameterDefaults(&Parameters);

  /* Open parameter file. */

  if ((fptr = fopen(ParameterFile, "r")) == NULL) {
    fprintf(stderr, "%s: error opening ParameterFile %s\n", myname, 
	    ParameterFile);
    exit(EXIT_FAILURE);
  }

  /* Read parameters from ParameterFile. */

  if (ReadParameterFile(fptr, &Parameters) == FAIL) {
    fprintf(stderr, "Error while reading ParameterFile %s\n", ParameterFile);
    exit(EXIT_FAILURE);
  }

  /* Read cosmology parameters. */

  rewind(fptr);
  if (CosmologyReadParameters(fptr) == FAIL) {
    fprintf(stderr, "Error in ReadCosmologyParameters.\n");;
    return FAIL;
  }

  /* Read power spectrum parameters. */

  rewind(fptr);
  if (ReadPowerSpectrumParameters(fptr) == FAIL) {
    fprintf(stderr, "Error in ReadCosmologyParameters.\n");;
    return FAIL;
  }

  /* Close parameter file. */

  fclose(fptr);

  /* Read parameters from SubGridParameter file. */

  if (SubGridParameterFile != NULL) {
    if ((fptr = fopen(SubGridParameterFile, "r")) == NULL) {
      fprintf(stderr, "%s: error opening SubGridParameterFile %s\n", myname, 
	      SubGridParameterFile);
      exit(EXIT_FAILURE);
    }
    SubGridParameters = new parmstruct;
    SetParameterDefaults(SubGridParameters);
    ReadParameterFile(fptr, SubGridParameters);
    fclose(fptr);
  }

  /* Initialize the power spectrum (and set amplitude) at z=0. */

  Redshift = 0.0;
  InitializePowerSpectrum();

  /* Generate a look-up table at the initial redshift. */

  Redshift = InitialRedshift;
  MakePowerSpectrumLookUpTable();

  /* Generate the fields and particles. */

  GenerateRealization(&Parameters, SubGridParameters);

  /* done */

#ifdef USE_JBMEM
  long long bcg,bmg,bcl,bml;
  bcg = jb::mem.bytes_global();
  bmg = jb::mem.bytes_global_high();
  bcl = jb::mem.bytes_local_max();
  bml = jb::mem.bytes_local_max_high();
  printf ("mem-global-curr    %lld\n", bcg);
  printf ("mem-global-high    %lld\n", bmg);
  printf ("mem-local-max-curr %lld\n", bcl);
  printf ("mem-local-max-high %lld\n", bml);
#endif

  printf("successful completion.\n");
  exit(EXIT_SUCCESS);

}

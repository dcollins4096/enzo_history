/***********************************************************************
/
/  READS POWER SPECTRUM PARAMETERS FROM INPUT FILE
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
#include "PowerSpectrumParameters.h"

int ReadPowerSpectrumParameters(FILE *fptr)
{

  int i;
  char line[MAX_LINE_LENGTH], *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  /* Set defaults. */

  sigma8          = 0.6;
  PrimordialIndex = 1.0;
  Gamma           = 0.25;
  kcutoff         = 0.0;
  RandomSeed      = -123456789;
  kmin            = 1.0e-3;     /* in Mpc^-1 */
  kmax            = 1.0e+5;     /* in Mpc^-1 */
  NumberOfkPoints = 10000;
  PowerSpectrumType = 1;
  WDMPartMass = 1.0;
  WDMg_x = 1.5;


  for (i = 0; i < 2; i++)
    PSFileName[i] = NULL;

  for (i = 0; i < MAX_SPECIES; i++)
    PSLookUpTable[i] = NULL;

  Normalization   = 1.0;
  GrowthFactor    = 1.0;
  Redshift        = 0.0;
  TophatRadius    = 0.0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    int ret = 0;

    /* read parameters */
    
    ret += sscanf(line, "PowerSpectrumSigma8 = %f", &sigma8);
    ret += sscanf(line, "PowerSpectrumPrimordialIndex = %f", &PrimordialIndex);
    ret += sscanf(line, "PowerSpectrumGamma = %f", &Gamma);
    ret += sscanf(line, "PowerSpectrumRandomSeed = %d", &RandomSeed);
    ret += sscanf(line, "PowerSpectrumkcutoff = %f", &kcutoff);
    ret += sscanf(line, "PowerSpectrumkmin = %f", &kmin);
    ret += sscanf(line, "PowerSpectrumkmax = %f", &kmax);
    ret += sscanf(line, "PowerSpectrumNumberOfkPoints = %d", &NumberOfkPoints);
    ret += sscanf(line, "PowerSpectrumType = %d", &PowerSpectrumType);
    ret += sscanf(line, "PowerSpectrumWDMParticleMass = %f",&WDMPartMass);
    ret += sscanf(line, "PowerSpectrumWDMDegreesOfFreedom = %f",&WDMg_x);

    if (sscanf(line, "PowerSpectrumFileNameRedshiftZero = %s", dummy) == 1)
      PSFileName[0] = dummy;

    if (sscanf(line, "PowerSpectrumFileNameInitialRedshift = %s", dummy) == 1)
      PSFileName[1] = dummy;

    /* If the dummy char space was used, then make another. */

    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#' && 
	strstr(line, "PowerSpectrum"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }

  /* Error check. */

  if (RandomSeed > -1) {
    fprintf(stderr, "RandomSeed must be a negative integer.\n");
    return FAIL;
  }
  
  return SUCCESS;
}

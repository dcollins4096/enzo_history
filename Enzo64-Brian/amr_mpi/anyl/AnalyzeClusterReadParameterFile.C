/***********************************************************************
/
/  READS A PARAMETER FILE FOR ANALYZE CLUSTER
/
/  written by: Greg Bryan
/  date:       August, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "AnalyzeClusters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
 
 
int AnalyzeClusterReadParameterFile(char *filename, int &NumberOfCenters,
				    float *CenterList[],
				    AnalyzeClusterParameters *parm)
 
{
 
  int dim, ret, j;
  char line[MAX_LINE_LENGTH], *char_dummy = new char[MAX_LINE_LENGTH];
  float center[MAX_DIMENSION], float_dummy;
 
  /* Set default vaules. */
 
  float BoxSize = 1;
  if (ComovingCoordinates)
    BoxSize = ComovingBoxSize/HubbleConstantNow;
  char *CenterListName = NULL;
  center[0]   = FLOAT_UNDEFINED;
 
  parm->rinner      = 0.0001*BoxSize;
  parm->router      = 0.1*BoxSize;
  parm->npoints     = 16;
  parm->virial_dens = FLOAT_UNDEFINED;
  parm->MeanVelocityVirialFraction = 1.0;
  parm->ColdTemperatureCutoff = 15000;    // in K
  parm->ColdTemperatureCutoffVirialFraction = FLOAT_UNDEFINED;
  parm->VirialTemperatureNormalization = 1.0;  // M-T-z normalization
  parm->LowerDensityCutoff     = 1.0e14;  /* in solar masses/Mpc^3 */
  parm->UpperDensityCutoff     = 1.0e35;  /* a big number */
  parm->ComputeDiskInformation = FALSE;
  parm->DiskImageSize          = 100;
  parm->DiskRadius             = 0.2;  // as a fraction of virial radius
  parm->XrayLowerCutoffkeV     = 0.5;
  parm->XrayUpperCutoffkeV     = 2.5;
  parm->XrayTableFileName      = NULL;
  parm->ComputeClumpingFactor  = FALSE;
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    center[dim] = FLOAT_UNDEFINED;
 
  /* Open file. */
 
  FILE *fptr;
  if ((fptr = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "error opening file %s\n", filename);
    exit(EXIT_FAILURE);
  }
 
  /* Read file. */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    ret += sscanf(line, "Rinner = %"FSYM, &parm->rinner);
    ret += sscanf(line, "Router = %"FSYM, &parm->router);
    ret += sscanf(line, "CenterPosition = %"FSYM" %"FSYM" %"FSYM,
		  center, center+1, center+2);
    ret += sscanf(line, "NumberOfPoints = %"ISYM, &parm->npoints);
    ret += sscanf(line, "VirialDensity = %"FSYM, &parm->virial_dens);
    ret += sscanf(line, "MeanVelocityVirialFraction = %"FSYM,
		  &parm->MeanVelocityVirialFraction);
    ret += sscanf(line, "ColdTemperatureCutoff = %"FSYM,
		  &parm->ColdTemperatureCutoff);
    ret += sscanf(line, "ColdTemperatureCutoffVirialFraction = %"FSYM,
		  &parm->ColdTemperatureCutoffVirialFraction);
    ret += sscanf(line, "VirialTemperatureNormalization = %"FSYM,
		  &parm->VirialTemperatureNormalization);
    ret += sscanf(line, "LowerDensityCutoff = %"FSYM,
		  &parm->LowerDensityCutoff);
    ret += sscanf(line, "UpperDensityCutoff = %"FSYM,
		  &parm->UpperDensityCutoff);
    ret += sscanf(line, "ComputeDiskInformation = %"ISYM,
		  &parm->ComputeDiskInformation);
    ret += sscanf(line, "DiskImageSize = %"ISYM, &parm->DiskImageSize);
    ret += sscanf(line, "DiskRadius = %"FSYM, &parm->DiskRadius);
    ret += sscanf(line, "XrayLowerCutoffkeV = %"FSYM,
		  &parm->XrayLowerCutoffkeV);
    ret += sscanf(line, "XrayUpperCutoffkeV = %"FSYM,
		  &parm->XrayUpperCutoffkeV);
    ret += sscanf(line, "ComputeClumpingFactor = %"ISYM,
		  &parm->ComputeClumpingFactor);
 
    if (sscanf(line, "CenterListName = %s", char_dummy) == 1)
      CenterListName = char_dummy;
 
    if (sscanf(line, "XrayTableFileName = %s", char_dummy) == 1)
      parm->XrayTableFileName = char_dummy;
 
    /* If the dummy char space was used, then make another. */
 
    if (*char_dummy != 0) {
      char_dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s", line);
 
  }
 
  /* Error check. */
 
  if (parm->rinner/BoxSize > 1 || parm->router/BoxSize > 1) {
    fprintf(stderr, "Rinner or Router > BoxSize.\n");
    exit(EXIT_FAILURE);
  }
 
  /* Close file. */
 
  fclose(fptr);
 
  /* If a CenterListName was specified, read this file, otherwise copy
     Center into CenterList. */
 
  if (CenterListName == NULL) {
 
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CenterList[dim] = new float[1];
      CenterList[dim][0] = center[dim];
    }
    NumberOfCenters = 1;
 
  } else {
 
    /* Open CenterListName. */
 
    if ((fptr = fopen(CenterListName, "r")) == NULL) {
      fprintf(stderr, "error opening CenterListName %s\n", CenterListName);
      exit(EXIT_FAILURE);
    }
 
    /* Count lines and allocate space. */
	
    NumberOfCenters = 0;
    while (fscanf(fptr, "%"FSYM, &float_dummy) > 0) {
      NumberOfCenters++;
      fgets(char_dummy, MAX_LINE_LENGTH, fptr);
    }
    if (debug) printf("NumberOfCenters = %"ISYM"\n", NumberOfCenters);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      CenterList[dim] = new float[NumberOfCenters];
    rewind(fptr);
 
    /* Read data. */
 
    j = 0;
    while (fscanf(fptr, "%"FSYM" %"FSYM" %"FSYM, CenterList[0]+j,
		  CenterList[1]+j, CenterList[2]+j) == 3) {
      fgets(char_dummy, MAX_LINE_LENGTH, fptr);  // get rid of rest of line
      j++;
    }
    if (j != NumberOfCenters) {
      fprintf(stderr, "Counting error (%"ISYM"/%"ISYM") in %s\n", j,
	      NumberOfCenters, CenterListName);
      exit(EXIT_FAILURE);
    }
 
  } // end: if (CenterListName == NULL)
 
  return SUCCESS;
}

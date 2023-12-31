/***********************************************************************
/
/  READS RADIATIVE TRANSFER PARAMETERS FROM INPUT FILE
/
/  written by: Tom Abel
/  date:       April, 2004
/  modified1:
/
/  PURPOSE: 
/
/  NOTE: modeled after CosmologyReadParameters
/
************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int RadiativeTransferReadParameters(FILE *fptr)
{

  int i;

  char line[MAX_LINE_LENGTH], *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  /* Set defaults. */

  RadiationPressure           = FALSE;             // off
  PhotonTime                  = 0; 
  dtPhoton                    = 0.1;
  for (i = 0; i < 4; i++) {
    EscapedPhotonCount[i] = 0.0;
    TotalEscapedPhotonCount[i] = 0.0;
  }
  PhotonEscapeFilename   = NULL;
  GlobalRadiationSources = new RadiationSourceEntry;
  GlobalRadiationSources->NextSource = NULL;
  GlobalRadiationSources->PreviousSource = NULL;
  SourceClusteringTree = NULL;
  OldSourceClusteringTree = NULL;

  RadiativeTransferSourceRadius               = 0;
  RadiativeTransferPropagationSpeedFraction   = 1.0;
  RadiativeTransferPropagationDistance        = 0.1;
  RadiativeTransferCoupledRateSolver          = TRUE;
  RadiativeTransferOpticallyThinH2            = TRUE;
  RadiativeTransferSplitPhotonRadius          = FLOAT_UNDEFINED; // kpc
  RadiativeTransferRaysPerCell                = 5.1;
  RadiativeTransferInitialHEALPixLevel        = 3;
  RadiativeTransferPhotonEscapeRadius         = 0.0;   // kpc
  RadiativeTransferInterpolateField           = FALSE;
  RadiativeTransferSourceClustering           = FALSE;
  RadiativeTransferPhotonMergeRadius          = 10.0;
  RadiativeTransferTimestepVelocityLimit      = 100.0; // km/s
  RadiativeTransferPeriodicBoundary           = FALSE;
  RadiativeTransferHIIRestrictedTimestep      = FALSE;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    int ret = 0;

    /* read parameters */
    
    ret += sscanf(line, "RadiativeTransferRadiationPressure = %"ISYM, \
		  &RadiationPressure);
    ret += sscanf(line, "RadiativeTransferSourceRadius = %"FSYM, 
		  &RadiativeTransferSourceRadius);
    ret += sscanf(line, "RadiativeTransferPropagationSpeedFraction = %"FSYM, 
		  &RadiativeTransferPropagationSpeedFraction);
    ret += sscanf(line, "RadiativeTransferPropagationDistance = %"FSYM, 
		  &RadiativeTransferPropagationDistance);
    ret += sscanf(line, "RadiativeTransferCoupledRateSolver = %"ISYM, 
		  &RadiativeTransferCoupledRateSolver);
    ret += sscanf(line, "RadiativeTransferOpticallyThinH2 = %"ISYM, 
		  &RadiativeTransferOpticallyThinH2);
    ret += sscanf(line, "RadiativeTransferPeriodicBoundary = %"ISYM, 
		  &RadiativeTransferPeriodicBoundary);
    ret += sscanf(line, "RadiativeTransferSplitPhotonRadius = %"FSYM, 
		  &RadiativeTransferSplitPhotonRadius);
    ret += sscanf(line, "RadiativeTransferRaysPerCell = %"FSYM, 
		  &RadiativeTransferRaysPerCell);
    ret += sscanf(line, "RadiativeTransferTimestepVelocityLimit = %"FSYM, 
		  &RadiativeTransferTimestepVelocityLimit);
    ret += sscanf(line, "RadiativeTransferInitialHEALPixLevel = %"ISYM, 
		  &RadiativeTransferInitialHEALPixLevel);
    ret += sscanf(line, "RadiativeTransferPhotonEscapeRadius = %"FSYM, 
		  &RadiativeTransferPhotonEscapeRadius);
    ret += sscanf(line, "RadiativeTransferInterpolateField = %"ISYM, 
		  &RadiativeTransferInterpolateField);
    ret += sscanf(line, "RadiativeTransferSourceClustering = %"ISYM, 
		  &RadiativeTransferSourceClustering);
    ret += sscanf(line, "RadiativeTransferPhotonMergeRadius = %"FSYM, 
		  &RadiativeTransferPhotonMergeRadius);
    ret += sscanf(line, "RadiativeTransferHIIRestrictedTimestep = %"ISYM, 
		  &RadiativeTransferHIIRestrictedTimestep);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#' && 
	MyProcessorNumber == ROOT_PROCESSOR &&
	strstr(line, "RadiativeTransfer") && !strstr(line, "RadiativeTransfer "))
      fprintf(stderr, "warning: the following parameter line was not "
	      "interpreted:\n%s\n", line);
  }

  /* Error check */

  /* Check if H2 cooling is turned on for Lyman-Werner radiation. */

  if (RadiativeTransferOpticallyThinH2 && MultiSpecies < 2 &&
      MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(stderr, "Warning: optically thin Lyman-Werner radiation turned on "
	    "without H2 cooling.  Setting LW radiation OFF.\n");
    RadiativeTransferOpticallyThinH2 = FALSE;
  }

  delete [] dummy;
  
  return SUCCESS;
}

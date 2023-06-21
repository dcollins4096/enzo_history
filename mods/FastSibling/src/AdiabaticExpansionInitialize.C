/***********************************************************************
/
/  INITIALIZE AN ADIABATIC EXPANSION TEST
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"

/* function prototypes */

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int AdiabaticExpansionInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid)
{
  char *DensName = "Density";
  char *TEName   = "Total Energy";
  char *GEName   = "Gas Energy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";

  /* declarations */

  char line[MAX_LINE_LENGTH];
  int ret;

  /* Error check. */

  if (!ComovingCoordinates) {
    fprintf(stderr, "ComovingCoordinates must be TRUE!\n");
    return FAIL;
  }

  /* set default parameters */

  float AdiabaticExpansionOmegaBaryonNow     = 1.0;  // standard
  float AdiabaticExpansionOmegaCDMNow        = 0.0;  // no dark matter
  float AdiabaticExpansionInitialTemperature = 200;  // degrees K
  float AdiabaticExpansionInitialVelocity    = 100;  // km/s

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "AdiabaticExpansionOmegaBaryonNow = %f", 
		  &AdiabaticExpansionOmegaBaryonNow);
    ret += sscanf(line, "AdiabaticExpansionOmegaCDMNow = %f", 
		  &AdiabaticExpansionOmegaCDMNow);
    ret += sscanf(line, "AdiabaticExpansionInitialTemperature = %f", 
		  &AdiabaticExpansionInitialTemperature);
    ret += sscanf(line, "AdiabaticExpansionInitialVelocity = %f", 
		  &AdiabaticExpansionInitialVelocity);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "AdiabaticExpansion"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }

  /* Get the cosmology units so we can convert temperature later. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
        VelocityUnits;
  if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		&TimeUnits, &VelocityUnits, InitialTimeInCodeUnits) == FAIL) {
    fprintf(stderr, "Error in CosmologyGetUnits.\n");
    return FAIL;
  }

  /* Put inputs in a form that will be understood by InitializeUniformGrid. */

  float InitialVels[MAX_DIMENSION], InitialTotalEnergy, InitialGasEnergy;
  InitialGasEnergy = AdiabaticExpansionInitialTemperature/TemperatureUnits /
    (Gamma - 1.0);
  InitialTotalEnergy = InitialGasEnergy;
  InitialVels[0] = AdiabaticExpansionInitialVelocity/VelocityUnits*1.0e5;
  InitialTotalEnergy += 0.5*pow(InitialVels[0],2);
  for (int dim = 1; dim < MAX_DIMENSION; dim++)
    InitialVels[dim] = 0.0;

  /* set up grid */

  if (TopGrid.GridData->InitializeUniformGrid(
					      AdiabaticExpansionOmegaBaryonNow,
					      InitialTotalEnergy, 
					      InitialGasEnergy, InitialVels
					      ) == FAIL) {
    fprintf(stderr, "Error in InitializeUniformGrid.\n");
    return FAIL;
  }

  /* set up field names and units */

  int i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;

  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "AdiabaticExpansionOmegaBaryonNow     = %f\n", 
	    AdiabaticExpansionOmegaBaryonNow);
    fprintf(Outfptr, "AdiabaticExpansionOmegaCDMNow        = %f\n", 
	    AdiabaticExpansionOmegaCDMNow);
    fprintf(Outfptr, "AdiabaticExpansionInitialTemperature = %f\n", 
	    AdiabaticExpansionInitialTemperature);
    fprintf(Outfptr, "AdiabaticExpansionInitialVelocity    = %f\n\n", 
	    AdiabaticExpansionInitialVelocity);
  }

  return SUCCESS;
}

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
/  INITIALIZE AN MHD ADIABATIC EXPANSION TEST
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

int MHDAdiabaticExpansionInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
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
  float AdiabaticExpansionBx                 = 0.0;
  float AdiabaticExpansionBy                 = 0.0;
  float AdiabaticExpansionBz                 = 0.0;


  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "AdiabaticExpansionOmegaBaryonNow = %"FSYM, 
		  &AdiabaticExpansionOmegaBaryonNow);
    ret += sscanf(line, "AdiabaticExpansionOmegaCDMNow = %"FSYM, 
		  &AdiabaticExpansionOmegaCDMNow);
    ret += sscanf(line, "AdiabaticExpansionInitialTemperature = %"FSYM, 
		  &AdiabaticExpansionInitialTemperature);
    ret += sscanf(line, "AdiabaticExpansionInitialVelocity = %"FSYM, 
		  &AdiabaticExpansionInitialVelocity);
    ret += sscanf(line,"AdiabaticExpansionBx = %"FSYM,
                  &AdiabaticExpansionBx);
    ret += sscanf(line,"AdiabaticExpansionBy = %"FSYM,
                  &AdiabaticExpansionBy);
    ret += sscanf(line,"AdiabaticExpansionBz = %"FSYM,
                  &AdiabaticExpansionBz);

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
  float InitialMagneticField[MAX_DIMENSION];
  InitialGasEnergy = AdiabaticExpansionInitialTemperature/TemperatureUnits /
    (Gamma - 1.0);
  InitialTotalEnergy = InitialGasEnergy;
  InitialVels[0] = AdiabaticExpansionInitialVelocity/VelocityUnits*1.0e5;
  InitialTotalEnergy += 0.5*POW(InitialVels[0],2);  
  for (int dim = 1; dim < MAX_DIMENSION; dim++)
    InitialVels[dim] = 0.0;
  InitialMagneticField[0] = AdiabaticExpansionBx;
  InitialMagneticField[1] = AdiabaticExpansionBy;
  InitialMagneticField[2] = AdiabaticExpansionBz;
   InitialTotalEnergy = InitialTotalEnergy
    +0.5*(InitialMagneticField[0]*InitialMagneticField[0]
         + InitialMagneticField[1]* InitialMagneticField[1]
         + InitialMagneticField[2]*InitialMagneticField[2]);

  /* set up grid */

  if (TopGrid.GridData->InitializeUniformGrid(
					      AdiabaticExpansionOmegaBaryonNow,
					      InitialTotalEnergy, 
					      InitialGasEnergy,InitialVels,
					      InitialMagneticField 
                                                  ) == FAIL) {
    fprintf(stderr, "Error in MHDInitializeUniformGrid.\n");
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

  MHDcLabel[0] = "MagneticField_C_1";
  MHDcLabel[1] = "MagneticField_C_2";
  MHDcLabel[2] = "MagneticField_C_3";
                                                                                
  MHDLabel[0] = "MagneticField_F_1";
  MHDLabel[1] = "MagneticField_F_2";
  MHDLabel[2] = "MagneticField_F_3";
                                                                                
  MHDeLabel[0] = "ElectricField_1";
  MHDeLabel[1] = "ElectricField_2";
  MHDeLabel[2] = "ElectricField_3";
                                                                                
  CurrentLabel[0] = "Current_1";
  CurrentLabel[1] = "Current_2";
  CurrentLabel[2] = "Current_3";

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "AdiabaticExpansionOmegaBaryonNow     = %f\n", 
	    AdiabaticExpansionOmegaBaryonNow);
    fprintf(Outfptr, "AdiabaticExpansionOmegaCDMNow        = %f\n", 
	    AdiabaticExpansionOmegaCDMNow);
    fprintf(Outfptr, "AdiabaticExpansionInitialTemperature = %f\n", 
	    AdiabaticExpansionInitialTemperature);
    fprintf(Outfptr, "AdiabaticExpansionInitialVelocity    = %f\n", 
	    AdiabaticExpansionInitialVelocity);
    fprintf(Outfptr, "AdiabaticExpansionBx = %f\n",
            AdiabaticExpansionBx);
    fprintf(Outfptr, "AdiabaticExpansionBy = %f\n",
            AdiabaticExpansionBy);
    fprintf(Outfptr, "AdiabaticExpansionBz = %f\n\n",
            AdiabaticExpansionBz);
  }

  return SUCCESS;
}

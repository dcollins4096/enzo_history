//  INITIALIZE A MHD Caustic

// This routine intializes a new simulation based on the parameter file.
//

#include <string.h>
#include <stdio.h>
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

int MHDCausticInitialize(FILE *fptr, FILE *Outfptr, 
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

 

  /* set default parameters */

  float CausticVelocityAmplitude = 1/(2*3.1415926);
  float CausticInitialPressure   = 1;  
  float CausticInitialDensity    = 1;
  
  float CausticMagneticfieldx     = 0.0;
  float CausticMagneticfieldy     = 0.0;
  float CausticMagneticfieldz     = 0.0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "CausticVelocityAmplitude = %"FSYM,  
		  &CausticVelocityAmplitude);
    ret += sscanf(line, "CausticInitialPressure  = %"FSYM, 
		  &CausticInitialPressure );
    ret += sscanf(line, "CausticInitialDensity = %"FSYM, 
		  &CausticInitialDensity);
    ret += sscanf(line, "CausticMagneticfieldx = %"FSYM, 
		  &CausticMagneticfieldx);
    ret += sscanf(line, "CausticMagneticfieldy = %"FSYM, 
		  &CausticMagneticfieldy);
    ret += sscanf(line, "CausticMagneticfieldz  = %"FSYM, 
		  &CausticMagneticfieldz );



    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "Caustic"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }

  /* set up grid */

  if (TopGrid.GridData->MHDCausticInitializeGrid(
					  CausticVelocityAmplitude,
				      CausticInitialPressure,
					  CausticInitialDensity,
					  CausticMagneticfieldx,
					  CausticMagneticfieldy,
					  CausticMagneticfieldz
						       ) == FAIL) {
    fprintf(stderr, "Error in CausticInitializeGrid.\n");
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
    fprintf(Outfptr, "CausticVelocityAmplitude         = %f\n",
	    CausticVelocityAmplitude);
    fprintf(Outfptr, "CausticInitialPressure     = %f\n", 
	    CausticInitialPressure);
    fprintf(Outfptr, "CausticInitialDensity    = %f\n", 
	    CausticInitialDensity);
    fprintf(Outfptr, "CausticMagneticfieldx        = %f\n", 
	    CausticMagneticfieldx);
    fprintf(Outfptr, "CausticMagneticfieldy   = %f\n", 
	    CausticMagneticfieldy);
    fprintf(Outfptr, "CausticMagneticfieldz = %f\n\n", 
	   CausticMagneticfieldz);
  }

  return SUCCESS;
}

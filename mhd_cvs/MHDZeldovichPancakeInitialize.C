//  INITIALIZE A MHD ZELDOVICH PANCAKE

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

int MHDZeldovichPancakeInitialize(FILE *fptr, FILE *Outfptr, 
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
  //  return FAIL;
  }  

  if (!SelfGravity)
    fprintf(stderr, "MHDZeldovichPancake: gravity is off!?!\n");
  if (CellFlaggingMethod[0] < 2)
    fprintf(stderr, "MHDZeldovichPancake: check CellFlaggingMethod.\n");

  /* set default parameters */

  int   ZeldovichPancakeDirection          = 0;    // along the x-axis
  float ZeldovichPancakeCentralOffset      = 0.0;  // no offset
  float ZeldovichPancakeOmegaBaryonNow     = 1.0;  // standard
  float ZeldovichPancakeOmegaCDMNow        = 0.0;  // no dark matter
  float ZeldovichPancakeCollapseRedshift   = 1.0;  // free parameter
  float ZeldovichPancakeInitialTemperature = 100;  // whatever
  
  float ZeldovichPancakeMagneticfieldx     = 0.0;
  float ZeldovichPancakeMagneticfieldy     = 0.0;
  float ZeldovichPancakeMagneticfieldz     = 0.0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "ZeldovichPancakeDirection = %d", 
		  &ZeldovichPancakeDirection);
    ret += sscanf(line, "ZeldovichPancakeCentralOffset = %"FSYM, 
		  &ZeldovichPancakeCentralOffset);
    ret += sscanf(line, "ZeldovichPancakeOmegaBaryonNow = %"FSYM, 
		  &ZeldovichPancakeOmegaBaryonNow);
    ret += sscanf(line, "ZeldovichPancakeOmegaCDMNow = %"FSYM, 
		  &ZeldovichPancakeOmegaCDMNow);
    ret += sscanf(line, "ZeldovichPancakeCollapseRedshift = %"FSYM, 
		  &ZeldovichPancakeCollapseRedshift);
    ret += sscanf(line, "ZeldovichPancakeInitialTemperature = %"FSYM, 
		  &ZeldovichPancakeInitialTemperature);

	ret += sscanf(line, "ZeldovichPancakeMagneticfieldx = %"FSYM, 
		  &ZeldovichPancakeMagneticfieldx);
    ret += sscanf(line, "ZeldovichPancakeMagneticfieldy = %"FSYM, 
		  &ZeldovichPancakeMagneticfieldy);
    ret += sscanf(line, "ZeldovichPancakeMagneticfieldz = %"FSYM, 
		  &ZeldovichPancakeMagneticfieldz);



    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "MHDZeldovichPancake"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }

  /* set up grid */

  if (TopGrid.GridData->MHDZeldovichPancakeInitializeGrid(
					  ZeldovichPancakeDirection,
					  ZeldovichPancakeCentralOffset,
					  ZeldovichPancakeOmegaBaryonNow,
					  ZeldovichPancakeOmegaCDMNow,
					  ZeldovichPancakeCollapseRedshift,
					  ZeldovichPancakeInitialTemperature,
					  ZeldovichPancakeMagneticfieldx,
					  ZeldovichPancakeMagneticfieldy,
					  ZeldovichPancakeMagneticfieldz
						       ) == FAIL) {
    fprintf(stderr, "Error in ZeldovichPancakeInitializeGrid.\n");
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
    fprintf(Outfptr, "ZeldovichPancakeDirection          = %d\n",
	    ZeldovichPancakeDirection);
    fprintf(Outfptr, "ZeldovichPancakeCentralOffset      = %f\n", 
	    ZeldovichPancakeCentralOffset);
    fprintf(Outfptr, "ZeldovichPancakeOmegaBaryonNow     = %f\n", 
	    ZeldovichPancakeOmegaBaryonNow);
    fprintf(Outfptr, "ZeldovichPancakeOmegaCDMNow        = %f\n", 
	    ZeldovichPancakeOmegaCDMNow);
    fprintf(Outfptr, "ZeldovichPancakeCollapseRedshift   = %f\n", 
	    ZeldovichPancakeCollapseRedshift);
    fprintf(Outfptr, "ZeldovichPancakeInitialTemperature = %f\n\n", 
	    ZeldovichPancakeInitialTemperature);
  }

  return SUCCESS;
}

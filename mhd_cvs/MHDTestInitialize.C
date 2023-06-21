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

int MHDTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData, ExternalBoundary &Exterior)
{

  char *DensName = "Density";
  
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  
  DataLabel[0] = DensName;
  DataLabel[1] = TEName;
  DataLabel[2] = Vel1Name;
  DataLabel[3] = Vel2Name;
  DataLabel[4] = Vel3Name;
  
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  
  float Density0 = .01, 
    Density1 = 1.,
    Energy0 = 178.571,  
    Energy1 = 0, 
    Velocity10 = 0, 
    Velocity11 = 0,
    Velocity20 = 0,
    Velocity21 = 0,
    Velocity30 = 0,
    Velocity31 = 0,
    B10 = 140,
    B20 = 0,
    B30 = 0,
    Radius = 3;

  MHDLabel[0] = "MagneticField_1";
  MHDLabel[1] = "MagneticField_2";
  MHDLabel[2] = "MagneticField_3";
  MHDeLabel[0] = "ElectricField_1";
  MHDeLabel[1] = "ElectricField_2";
  MHDeLabel[2] = "ElectricField_3";

  MHDUnits[0] = "Gauss";
  MHDUnits[1] = "Gauss";
  MHDUnits[2] = "Gauss";

  MHDeUnits[0] = "Gauss";
  MHDeUnits[1] = "Gauss";
  MHDeUnits[2] = "Gauss";
  

  if( TopGrid.GridData->MHDTestInitializeGrid(Density0, Density1,
					       Energy0,  Energy1,
					       Velocity10, Velocity11,
					       Velocity20, Velocity21,
					       Velocity30, Velocity31,
					       B10,B20,B30, Radius) == FAIL )
    {
      fprintf(stderr, "Shit.  Error in MHDBlastInitializeGrid.\n");
      return FAIL;
    }
  
  return SUCCESS;
}









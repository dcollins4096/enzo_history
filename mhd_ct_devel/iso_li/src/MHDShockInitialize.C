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

int MHDShockInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData, ExternalBoundary &Exterior)
{

  char line[MAX_LINE_LENGTH];
  int ret=0;
  int ShockTubeDirection = -1;

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

#ifdef HAOXU
if(DualEnergyFormalism == 1){
 char *GEName = "GasEnergy";
 DataLabel[5] = GEName;
 DataUnits[5] = NULL;
}

#endif //HAOXU


  MHDcLabel[0] = "MagneticField_C_1";
  MHDcLabel[1] = "MagneticField_C_2";
  MHDcLabel[2] = "MagneticField_C_3";

  MHDLabel[0] = "MagneticField_F_1";
  MHDLabel[1] = "MagneticField_F_2";
  MHDLabel[2] = "MagneticField_F_3";

  MHDeLabel[0] = "ElectricField_1";
  MHDeLabel[1] = "ElectricField_2";
  MHDeLabel[2] = "ElectricField_3";

  MHDcUnits[0] = "Four Pi Gauss";
  MHDcUnits[1] = "FourPiGauss";
  MHDcUnits[2] = "FourPiGauss";

  MHDUnits[0] = "Four Pi Gauss";
  MHDUnits[1] = "FourPiGauss";
  MHDUnits[2] = "FourPiGauss";

  MHDeUnits[0] = "FourPiGauss";
  MHDeUnits[1] = "FourPiGauss";
  MHDeUnits[2] = "FourPiGauss";

  CurrentLabel[0] = "Current_1";
  CurrentLabel[1] = "Current_2";
  CurrentLabel[2] = "Current_3";


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret += sscanf(line, "ShockTubeDirection = %d", &ShockTubeDirection);

  }

  if( ret == 0 ) {
    fprintf(stderr, "Shit!  need shock tube direction \n");
    return FAIL;

  }
  if( TopGrid.GridData->MHDShockInitializeGrid(ShockTubeDirection) == FAIL )
    {
      fprintf(stderr, "Shit.  Error in MHDShockInitializeGrid.\n");
      return FAIL;
    }
  
  return SUCCESS;
}









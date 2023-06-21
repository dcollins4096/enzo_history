

#include <math.h>
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

int MHDOrszagTangInit(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		      TopGridData &MetaData, ExternalBoundary &Exterior){


  float Pi = 3.14159265;
  float Density = 25.0/(36*Pi); 
  float Pressure = 5.0/(12*Pi);
  float V0= 1.0;
  float B0=1/sqrt(4*Pi);

  //The Isothermal initial conditions come from Mignone, JCP, 2007
  if( EquationOfState == 1 ){
    Density = 1;
    V0= IsothermalSoundSpeed ;
    B0 = IsothermalSoundSpeed * sqrt(3./5.);
  }

  int count=0;
  char *DensName = "Density";
  
  char *TEName = "Total_Energy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *GEName   = "GasEnergy";
  char *BXName = "MagneticField_C_1";
  char *BYName = "MagneticField_C_2";
  char *BZName = "MagneticField_C_3";

  DataLabel[count++] = DensName;
  if( EquationOfState == 0 )
    DataLabel[count++] = TEName;
  if (DualEnergyFormalism){
    DataLabel[count++] = GEName;
    //<dbg>
    fprintf(stderr,"MHDORszagTangInit: Dual Energy Formailism not completely installed.\n");
    return FAIL;
  }  
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;

  for (int i = 0; i < count; i++)
    DataUnits[i] = NULL;

  if( MHD_Used == TRUE ){
    DataLabel[count++] = BXName;
    DataLabel[count++] = BYName;
    DataLabel[count++] = BZName;

  }

#ifdef MHD

  MHDLabel[0] = "MagneticField_F_1";
  MHDLabel[1] = "MagneticField_F_2";
  MHDLabel[2] = "MagneticField_F_3";

  MHDeLabel[0] = "ElectricField_1";
  MHDeLabel[1] = "ElectricField_2";
  MHDeLabel[2] = "ElectricField_3";

  MHDcUnits[0] = "FourPiGauss";
  MHDcUnits[1] = "FourPiGauss";
  MHDcUnits[2] = "FourPiGauss";

  MHDeUnits[0] = "FourPiGauss";
  MHDeUnits[1] = "FourPiGauss";
  MHDeUnits[2] = "FourPiGauss";

  CurrentLabel[0] = "Current_1";
  CurrentLabel[1] = "Current_2";
  CurrentLabel[2] = "Current_3";
  
#endif //MHD

  if( TopGrid.GridData->MHDOrszagTangInitGrid(Density,Pressure,V0,B0) == FAIL ){
    fprintf(stderr, " Error in Grid_MHDOrzagTangInitGrid\n");
    return FAIL;
  }
  return SUCCESS;


}




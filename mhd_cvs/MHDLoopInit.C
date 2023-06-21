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

int MHDLoopInit(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData, ExternalBoundary &Exterior)
{

  float Pi = 3.14159265;
  float Density = 25.0/(36*Pi); 
  float Pressure = 5.0/(12*Pi);
  float Vx=0.0, Vy = 0.0, Vz = 0.0;
  float B0=1e-3;
  float R0 = 0.3;
  float Center[3];
  int CurrentAxis = 2;
  for(int dim =0;dim<3;dim++)
    Center[dim] = 0.5*(DomainLeftEdge[dim]+DomainRightEdge[dim]);

  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  #ifdef HAOXU
  char *GEName   = "GasEnergy";
#endif
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  
#ifdef ATHENA
  if( EquationOfState == 0 ){
    DataLabel[0] = DensName;
    DataLabel[1] = TEName;
    DataLabel[2] = Vel1Name;
    DataLabel[3] = Vel2Name;
    DataLabel[4] = Vel3Name;
  }else if (EquationOfState == 1){
    DataLabel[0] = DensName;
    DataLabel[1] = Vel1Name;
    DataLabel[2] = Vel2Name;
    DataLabel[3] = Vel3Name;
  }
#else //Athena
  DataLabel[0] = DensName;
  DataLabel[1] = TEName;
#ifdef HAOXU
   if (DualEnergyFormalism)
    DataLabel[5] = GEName;
#endif

  DataLabel[2] = Vel1Name;
  DataLabel[3] = Vel2Name;
  DataLabel[4] = Vel3Name;
#endif //Athena

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

  MHDcUnits[0] = "FourPiGauss";
  MHDcUnits[1] = "FourPiGauss";
  MHDcUnits[2] = "FourPiGauss";

  MHDUnits[0] = "FourPiGauss";
  MHDUnits[1] = "FourPiGauss";
  MHDUnits[2] = "FourPiGauss";

  MHDeUnits[0] = "FourPiGauss";
  MHDeUnits[1] = "FourPiGauss";
  MHDeUnits[2] = "FourPiGauss";

  CurrentLabel[0] = "Current_1";
  CurrentLabel[1] = "Current_2";
  CurrentLabel[2] = "Current_3";
  
  char line[MAX_LINE_LENGTH];
  int ret=0;

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    ret += sscanf(line, "MHDLoopDensity = %"PSYM, &Density);
    ret += sscanf(line, "MHDLoopPressure = %"PSYM, &Pressure);
    ret += sscanf(line, "MHDLoopVx = %"PSYM, &Vx);
    ret += sscanf(line, "MHDLoopVy = %"PSYM, &Vy);
    ret += sscanf(line, "MHDLoopVz = %"PSYM, &Vz);
    ret += sscanf(line, "MHDLoopB0 = %"PSYM, &B0);
    ret += sscanf(line, "MHDLoopR0 = %"PSYM, &R0);
    ret += sscanf(line, "MHDLoopCurrentAxis = %d", &CurrentAxis);
    ret += sscanf(line, "MHDLoopCenter = %"PSYM" %"PSYM" %"PSYM, 
		  Center, Center+1, Center+2);
    
  }  
  if( TopGrid.GridData->MHDLoopInitGrid(Density,Pressure,Vx,Vy, Vz, B0, R0,Center, CurrentAxis) == FAIL ){
    fprintf(stderr, " Shit, Man, the Loops all funky.\n");
    return FAIL;
  }

  return SUCCESS;


}





/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- EVOLVE A CONSTANT FIELD
/
/  written by: Daniel Reynolds
/  date:       November, 2006
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
#ifdef RAD_HYDRO

// This routine intializes a new simulation based on the parameter file.

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


// function prototypes
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits, 
		      float *TemperatureUnits, float *TimeUnits, 
		      float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int InitializeRateData(FLOAT Time);




int RadHydroStreamTestInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData)
{
  if (debug)
    fprintf(stdout,"Entering RadHydroStreamTestInitialize routine\n");

  char *DensName  = "Density";
  char *TEName    = "Total_Energy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName   = "Grey_Radiation_Energy";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *DeName    = "Electron_Density";

  // local declarations
  int dim;

  // make sure it is 3D
  if (MetaData.TopGridRank != 3) {
    printf("Cannot do Rad-Hydro Tests in %"ISYM" dimension(s)\n", 
	   MetaData.TopGridRank);
    return FAIL;
  }    


  // Setup and parameters:
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient radiation energy
  float RadHydroDensity         = 1.0;
  float RadHydroRadiationEnergy = 1.0e-10;


  // set the fluid boundaries to be periodic
  for (dim=0; dim<MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }


  // set up the grid
  if (TopGrid.GridData->RadHydroStreamTestInitializeGrid(
       RadHydroDensity, RadHydroRadiationEnergy) == FAIL) {
      fprintf(stderr, "Error in RadHydroStreamTestInitializeGrid.\n");
      return FAIL;
    }

  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = IEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  DataLabel[BaryonField++] = RadName;
  DataLabel[BaryonField++] = DeName;
  DataLabel[BaryonField++] = HIName;
  DataLabel[BaryonField++] = HIIName;
  DataLabel[BaryonField++] = HeIName;
  DataLabel[BaryonField++] = HeIIName;
  DataLabel[BaryonField++] = HeIIIName;

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;


  return SUCCESS;

}

#endif

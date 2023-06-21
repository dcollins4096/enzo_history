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




int RadHydroConstTestInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData)
{
  if (debug)
    fprintf(stdout,"Entering RadHydroConstTestInitialize routine\n");

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
  char line[MAX_LINE_LENGTH];
  int  dim, ret;

  // make sure it is 3D
  if (MetaData.TopGridRank != 3) {
    printf("Cannot do Rad-Hydro Tests in %"ISYM" dimension(s)\n", 
	   MetaData.TopGridRank);
    return FAIL;
  }    


  // Setup and parameters:
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient gas velocity - free parameter
  //  3. ambient gas temperature
  //  4. ambient radiation energy
  //  5. Hydrogen mass fraction 
  //  6. initial fraction HII
  //  7. initial fraction HeII
  //  8. initial fraction HeIII
  //  9. OmegaBaryonNow
  // 10. Number of chemical species
  // 11. mesh spacing
  float RadHydroX0Velocity           = 0.0;
  float RadHydroX1Velocity           = 0.0;
  float RadHydroX2Velocity           = 0.0;
  float RadHydroDensity              = 10.0;
  float RadHydroTemperature          = 1e3;
  float RadHydroRadiationEnergy      = 10.0;
  float RadHydroHydrogenMassFraction = 1.0;
  float RadHydroInitialFractionHII   = 0.0;
  float RadHydroInitialFractionHeII  = 0.0;
  float RadHydroInitialFractionHeIII = 0.0;
  float RadHydroOmegaBaryonNow       = 1.0;
  int   RadHydroChemistry            = 0;

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
   
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	
	ret = 0;

	// read relevant problem parameters
	ret += sscanf(line, "RadHydroVelocity = %"FSYM" %"FSYM" %"FSYM,
		      &RadHydroX0Velocity, &RadHydroX1Velocity, 
		      &RadHydroX2Velocity);
	ret += sscanf(line, "RadHydroChemistry = %"ISYM, 
		      &RadHydroChemistry);
	ret += sscanf(line, "RadHydroDensity = %"FSYM, 
		      &RadHydroDensity);
	ret += sscanf(line, "RadHydroTemperature = %"FSYM, 
		      &RadHydroTemperature);
	ret += sscanf(line, "RadHydroRadiationEnergy = %"FSYM, 
		      &RadHydroRadiationEnergy);
	if (RadHydroChemistry > 0)
	  ret += sscanf(line, "RadHydroInitialFractionHII = %"FSYM, 
			&RadHydroInitialFractionHII);
	if (RadHydroChemistry > 1) {
	  ret += sscanf(line, "RadHydroHydrogenMassFraction = %"FSYM, 
			&RadHydroHydrogenMassFraction);
	  ret += sscanf(line, "RadHydroInitialFractionHeII = %"FSYM, 
			&RadHydroInitialFractionHeII);
	  ret += sscanf(line, "RadHydroInitialFractionHeIII = %"FSYM, 
			&RadHydroInitialFractionHeIII);
	}
	ret += sscanf(line, "RadHydroOmegaBaryonNow = %"FSYM, 
		      &RadHydroOmegaBaryonNow);
	
      } // end input from parameter file
      
      fclose(RHfptr);
    }
  }


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
  if (TopGrid.GridData->RadHydroConstTestInitializeGrid(
       RadHydroChemistry, RadHydroDensity, RadHydroX0Velocity, 
       RadHydroX1Velocity, RadHydroX2Velocity, RadHydroTemperature, 
       RadHydroRadiationEnergy, RadHydroHydrogenMassFraction, 
       RadHydroInitialFractionHII, RadHydroInitialFractionHeII, 
       RadHydroInitialFractionHeIII, RadHydroOmegaBaryonNow) == FAIL) {
      fprintf(stderr, "Error in RadHydroConstTestInitializeGrid.\n");
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
  if (RadHydroChemistry > 0) {
    DataLabel[BaryonField++] = DeName;
    DataLabel[BaryonField++] = HIName;
    DataLabel[BaryonField++] = HIIName;
  }
  if (RadHydroChemistry > 1) {
    DataLabel[BaryonField++] = HeIName;
    DataLabel[BaryonField++] = HeIIName;
    DataLabel[BaryonField++] = HeIIIName;
  }

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;


  // Write parameters to output file
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "RadHydroVelocity = %"FSYM" %"FSYM" %"FSYM"\n", 
	    RadHydroX0Velocity, RadHydroX1Velocity, RadHydroX2Velocity);
    fprintf(Outfptr, "RadHydroChemistry = %"ISYM"\n", RadHydroChemistry);
    fprintf(Outfptr, "RadHydroDensity = %"FSYM"\n", RadHydroDensity);
    fprintf(Outfptr, "RadHydroTemperature = %"FSYM"\n", RadHydroTemperature);
    fprintf(Outfptr, "RadHydroRadiationEnergy = %"FSYM"\n", 
	    RadHydroRadiationEnergy);
    fprintf(Outfptr, "RadHydroHydrogenMassFraction = %"FSYM"\n", 
	    RadHydroHydrogenMassFraction);
    fprintf(Outfptr, "RadHydroInitialFractionHII = %"FSYM"\n", 
	    RadHydroInitialFractionHII);
    fprintf(Outfptr, "RadHydroInitialFractionHeII = %"FSYM"\n", 
	    RadHydroInitialFractionHeII);
    fprintf(Outfptr, "RadHydroInitialFractionHeIII = %"FSYM"\n", 
	    RadHydroInitialFractionHeIII);
    fprintf(Outfptr, "RadHydroOmegaBaryonNow = %"FSYM"\n", 
	    RadHydroOmegaBaryonNow);
  }

  return SUCCESS;

}

#endif

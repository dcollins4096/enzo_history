//  GRID CLASS (INITIALIZE THE GRID FOR A MHD Caustic Test)

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

#define TOLERANCE 1.0e-6

/* function prototypes */

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);


int grid::MHDCausticInitializeGrid(float CausticVelocityAmplitude,
				      float CausticInitialPressure,
					  float CausticInitialDensity,
					  float CausticMagneticfieldx,
					  float CausticMagneticfieldy,
					  float CausticMagneticfieldz)
{
  /* declarations */

  float nx;
  int   dim, field, i, index;
  const float Pi = 3.14159;

 
  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  int vel = NumberOfBaryonFields - 1;
    FieldType[NumberOfBaryonFields++] = Velocity2;
    FieldType[NumberOfBaryonFields++] = Velocity3;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

    /* Determine the size of the fields. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Allocate space for the fields. */

  for (field = 0; field < NumberOfBaryonFields; field++)
    BaryonField[field] = new float[size];

  
  /* Allocate space for the Magnetic fields */
    CenteredB[0] = new float[size];
    Current[0] = new float[size];
    MagneticField[0] = new float[MagneticSize[0]];
    ElectricField[0] = new float[ElectricSize[0]];
    CenteredB[1] = new float[size]; 
	Current[1] = new float[size];
	MagneticField[1] = new float[MagneticSize[1]];
    ElectricField[1] = new float[ElectricSize[1]]; 
	CenteredB[2] = new float[size]; 
	Current[2] = new float[size]; 
	MagneticField[2] = new float[MagneticSize[2]];
    ElectricField[2] = new float[ElectricSize[2]];
    DivB = new float[size];

     //Initialize Magnetic and Electric Fields

  for(i=0;i < size; i++) {
	  CenteredB[0][i] = CausticMagneticfieldx;
	  CenteredB[1][i] = CausticMagneticfieldy;
	  CenteredB[2][i] = CausticMagneticfieldz;
  }
  for(i=0;i<MagneticSize[0];i++) MagneticField[0][i] = CausticMagneticfieldx;
  for(i=0;i<MagneticSize[1];i++) MagneticField[1][i] = CausticMagneticfieldy;
  for(i=0;i<MagneticSize[2];i++) MagneticField[2][i] = CausticMagneticfieldz;
    
  for(i=0;i<ElectricSize[0];i++) ElectricField[0][i] = 0.0;
  for(i=0;i<ElectricSize[1];i++) ElectricField[1][i] = 0.0;
  for(i=0;i<ElectricSize[2];i++) ElectricField[2][i] = 0.0; 
  

   for (i = 0; i < size; i++) {
 
   
    /* Set density. */

    BaryonField[0][i] = CausticInitialDensity;

    /* Set total energy, gas energy. */

    BaryonField[1][i] = CausticInitialPressure / (Gamma-1.0);
    if (DualEnergyFormalism)
      BaryonField[2][i] = BaryonField[1][i]/BaryonField[0][i];
 
   BaryonField[1][i]
=BaryonField[1][i]+0.5*POW(CausticVelocityAmplitude*sin((((float)(i%GridDimension[0]))-DEFAULT_GHOST_ZONES)*(float)2*Pi/(GridDimension[0]-2*DEFAULT_GHOST_ZONES)),2)*BaryonField[0][i]
                + 0.5*POW(CenteredB[0][i],2)+0.5*POW(CenteredB[1][i],2)+0.5*POW(CenteredB[2][i],2);
    

/* wallshock
  BaryonField[1][i]  =BaryonField[1][i]+0.5*POW(CausticVelocityAmplitude,2)*BaryonField[0][i]
                + 0.5*POW(CenteredB[0][i],2)+0.5*POW(CenteredB[1][i],2)+0.5*POW(CenteredB[2][i],2);
*/
    /* Set velocities (some may be set to zero later -- see below). */


    BaryonField[vel][i] = CausticVelocityAmplitude*sin((((float)(i%GridDimension[0]))-DEFAULT_GHOST_ZONES)*(float)2*Pi/(GridDimension[0]-2*DEFAULT_GHOST_ZONES));
// for wallshock	BaryonField[vel][i] = CausticVelocityAmplitude;//*sin((((float)(i%GridDimension[0])-0.5)-3)*(float)2*Pi/(GridDimension[0]-6));
    BaryonField[vel+1][i]   = 0.0;
    BaryonField[vel+2][i]   = 0.0;
   
  }


  return SUCCESS;
}

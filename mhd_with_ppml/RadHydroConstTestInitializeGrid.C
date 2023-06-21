/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE PROTOSTELLAR CORE COLLAPSE TEST) 
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, June 2005.
/
/  PURPOSE: Sets the initial core region.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"


#define DEFAULT_MU 0.6 // for temperature field


// function prototypes
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits, 
		      float *TemperatureUnits, float *TimeUnits, 
		      float *VelocityUnits, FLOAT Time);



int grid::RadHydroConstTestInitializeGrid(int NumChemicals,
					  float DensityConstant, 
					  float VxConstant, 
					  float VyConstant, 
					  float VzConstant, 
					  float TempConstant, 
					  float EgConstant, 
					  float HydrogenMassFraction, 
					  float InitialFractionHII, 
					  float InitialFractionHeII, 
					  float InitialFractionHeIII, 
					  float OmegaBaryonNow)
{

  if (debug)
    fprintf(stdout,"Entering grid::RadHydroConstTestInitializeGrid routine\n");

  // ensure that we are performing 3D test problem (HYPRE requirement)
  if (GridRank != 3) {
    fprintf(stderr,"ERROR: Radiation Hydrodynamics tests MUST be 3D!\n");
    return FAIL;
  }
 
  // create necessary baryon fields
  int RhoNum, TENum, IENum, V0Num, V1Num, V2Num, EgNum, DeNum, 
    HINum, HIINum, HeINum, HeIINum, HeIIINum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++]  = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[IENum = NumberOfBaryonFields++] = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++] = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++] = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++] = Velocity3;
  FieldType[EgNum = NumberOfBaryonFields++] = RadiationFreq0;
  if (NumChemicals > 0) {
    FieldType[DeNum = NumberOfBaryonFields++]  = ElectronDensity;
    FieldType[HINum = NumberOfBaryonFields++]  = HIDensity;
    FieldType[HIINum = NumberOfBaryonFields++] = HIIDensity;
  }
  if (NumChemicals > 1) {
    FieldType[HeINum = NumberOfBaryonFields++]   = HeIDensity;
    FieldType[HeIINum = NumberOfBaryonFields++]  = HeIIDensity;    
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;    
  }

  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;  // no subgrids
 
  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Get cosmology units in order to convert temperature later
  float DensityUnits, LengthUnits, TempUnits, TimeUnits, VelUnits;
  DensityUnits = LengthUnits = TempUnits = TimeUnits = VelUnits = 1.0;
  if (ComovingCoordinates) 
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TempUnits, &TimeUnits, 
			  &VelUnits, InitialTimeInCodeUnits) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }
  if (debug) {
    fprintf(stdout,"  Internal Unit Conversion Factors:\n");
    fprintf(stdout,"        density = %g\n",DensityUnits);
    fprintf(stdout,"         length = %g\n",LengthUnits);
    fprintf(stdout,"    temperature = %g\n",TempUnits);
    fprintf(stdout,"           Time = %g\n",TimeUnits);
    fprintf(stdout,"       velocity = %g\n",VelUnits);
    fprintf(stdout,"   initial time = %g\n",InitialTimeInCodeUnits);
  }

  // compute size of fields
  int size = 1;
  for (int dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
 
  // allocate fields
  for (int field=0; field<NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];
 
  // set fluid density, total energy, [internal energy,] velocities, 
  // radiation energy, electron density, chemical species
  int i;
  float TEConstant, IEConstant, DeConstant, HIConstant, HIIConstant;
  float HeIConstant, HeIIConstant, HeIIIConstant;
  IEConstant = TempConstant/(TempUnits*DEFAULT_MU*(Gamma-1.0));
  TEConstant = IEConstant + 0.5/VelUnits/VelUnits*
    (VxConstant*VxConstant + VyConstant*VyConstant + VzConstant*VzConstant);
  // why does HII use OmegaMatterNow, OmegaBaryonNow and HubbleConstantNow, 
  // and HI, HeI, HeII, HeIII do not?
  if (ComovingCoordinates) 
    HIIConstant = InitialFractionHII * HydrogenMassFraction * 
      DensityConstant * sqrt(OmegaMatterNow) / 
      (OmegaBaryonNow*HubbleConstantNow);
  else
    HIIConstant = InitialFractionHII * HydrogenMassFraction * 
      DensityConstant;
  HIConstant = HydrogenMassFraction * DensityConstant - HIIConstant;
  HeIIConstant = InitialFractionHeII * DensityConstant * 
    4.0 * (1.0-HydrogenMassFraction);
  HeIIIConstant = InitialFractionHeIII * DensityConstant * 
    4.0 * (1.0-HydrogenMassFraction);
  HeIConstant = (1.0 - HydrogenMassFraction)*DensityConstant - 
    HeIIConstant - HeIIIConstant;
  DeConstant = HIIConstant + 0.25*HeIIConstant + 0.5*HeIIIConstant;
  for (i=0; i<size; i++) {
    BaryonField[RhoNum][i] = DensityConstant/DensityUnits;
    BaryonField[TENum][i]  = TEConstant;
    BaryonField[V0Num][i]  = VxConstant/VelUnits;
    BaryonField[V1Num][i]  = VyConstant/VelUnits;
    BaryonField[V2Num][i]  = VzConstant/VelUnits;
    BaryonField[EgNum][i]  = EgConstant;
    BaryonField[DeNum][i]  = DeConstant/DensityUnits;
  }
  if (DualEnergyFormalism)
    for (i=0; i<size; i++)
      BaryonField[IENum][i] = IEConstant;
  if (NumChemicals > 0)
    for (i=0; i<size; i++) {
      BaryonField[HINum][i]  = HIConstant/DensityUnits;
      BaryonField[HIINum][i] = HIIConstant/DensityUnits;
    }
  if (NumChemicals > 1)
    for (i=0; i<size; i++) {
      BaryonField[HeINum][i]   = HeIConstant/DensityUnits;
      BaryonField[HeIINum][i]  = HeIIConstant/DensityUnits;
      BaryonField[HeIIINum][i] = HeIIIConstant/DensityUnits;
    }

  if (debug) {
    fprintf(stdout,"\n  Initializing constant fields using CGS values:\n");
    fprintf(stdout,"        density = %g\n",DensityConstant);
    fprintf(stdout,"         energy = %g\n",TEConstant);
    fprintf(stdout,"     x-velocity = %g\n",VxConstant);
    fprintf(stdout,"     y-velocity = %g\n",VyConstant);
    fprintf(stdout,"     z-velocity = %g\n",VzConstant);
    fprintf(stdout,"      radiation = %g\n",EgConstant);
    fprintf(stdout,"      electrons = %g\n",DeConstant);
    fprintf(stdout,"            nHI = %g\n",HIConstant);
    fprintf(stdout,"           nHII = %g\n",HIIConstant);
    fprintf(stdout,"           nHeI = %g\n",HeIConstant);
    fprintf(stdout,"          nHeII = %g\n",HeIIConstant);
    fprintf(stdout,"         nHeIII = %g\n",HeIIIConstant);

    fprintf(stdout,"\n  Corresponding Enzo internal values:\n");
    fprintf(stdout,"        density = %g\n",BaryonField[RhoNum][1]);
    fprintf(stdout,"         energy = %g\n",BaryonField[TENum][1]);
    fprintf(stdout,"     x-velocity = %g\n",BaryonField[V0Num][1]);
    fprintf(stdout,"     y-velocity = %g\n",BaryonField[V1Num][1]);
    fprintf(stdout,"     z-velocity = %g\n",BaryonField[V2Num][1]);
    fprintf(stdout,"      radiation = %g\n",BaryonField[EgNum][1]);
    fprintf(stdout,"      electrons = %g\n",BaryonField[DeNum][1]);
    fprintf(stdout,"            nHI = %g\n",BaryonField[HINum][1]);
    fprintf(stdout,"           nHII = %g\n",BaryonField[HIINum][1]);
  }

  return SUCCESS;
}

/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE RAD-HYDRO IONIZATION TEST) 
/
/  written by: Daniel Reynolds
/  date:       July 2007
/
/  PURPOSE: 
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


// function prototypes
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);



int grid::RHIonizationTestInitializeGrid(int NumChemicals,
					 float DensityConstant, 
					 float VxConstant, 
					 float VyConstant, 
					 float VzConstant, 
					 float IEConstant, 
					 float EgConstant, 
					 float HydrogenMassFraction, 
					 float InitialFractionHII, 
					 float InitialFractionHeII, 
					 float InitialFractionHeIII, 
					 int   local)
{
#ifdef TRANSFER
//   if (debug)
//     fprintf(stdout,"Entering grid::RHIonizationTestInitializeGrid routine\n");

  // determine whether data should be allocated/initialized or not
  int NewData = TRUE;
  if ((ParallelRootGridIO == TRUE) && (local == 0))
    NewData = FALSE;

  // if grids allocated and already set up (i.e. restart), return
  if ((NumberOfBaryonFields > 5) && (BaryonField[5] != NULL))
    return SUCCESS;

  // create necessary baryon fields
  int RhoNum, TENum, IENum, V0Num, V1Num, V2Num, EgNum, DeNum, 
    HINum, HIINum, HeINum, HeIINum, HeIIINum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++]  = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[IENum = NumberOfBaryonFields++]  = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++]    = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++]    = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++]    = Velocity3;
  FieldType[EgNum = NumberOfBaryonFields++]    = RadiationFreq0;
  if (NumChemicals > 0) {
    FieldType[DeNum = NumberOfBaryonFields++]    = ElectronDensity;
    FieldType[HINum = NumberOfBaryonFields++]    = HIDensity;
    FieldType[HIINum = NumberOfBaryonFields++]   = HIIDensity;
  }
  if (NumChemicals == 3) {
    FieldType[HeINum = NumberOfBaryonFields++]   = HeIDensity;
    FieldType[HeIINum = NumberOfBaryonFields++]  = HeIIDensity;    
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
  }

  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;  // no subgrids

  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Get various units
  double MassUnits=1;
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr,"Error in GetUnits.\n");
    return FAIL;
  }
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(stdout,"  Internal Unit Conversion Factors:\n");
    fprintf(stdout,"         length = %g\n",LengthUnits);
    fprintf(stdout,"           mass = %lg\n",MassUnits);
    fprintf(stdout,"           time = %g\n",TimeUnits);
  }

  // compute size of fields
  int size = 1;
  for (int dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
 
  // allocate fields
  if (NewData == TRUE) {
    for (int field=0; field<NumberOfBaryonFields; field++)
      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[size];
    
    // set fluid density, total energy, [internal energy,] velocities, 
    // radiation energy, electron density, chemical species
    int i;
    float TEConstant = (IEConstant + 
			0.5*(VxConstant*VxConstant + 
			     VyConstant*VyConstant + 
			     VzConstant*VzConstant));
    float HIIConstant = InitialFractionHII*HydrogenMassFraction*DensityConstant;
    float HIConstant = HydrogenMassFraction*DensityConstant - HIIConstant;
    float HeIIConstant = InitialFractionHeII*DensityConstant*(1.0-HydrogenMassFraction);
    float HeIIIConstant = InitialFractionHeIII*DensityConstant*(1.0-HydrogenMassFraction);
    float HeIConstant = (1.0-HydrogenMassFraction)*DensityConstant 
      - HeIIConstant - HeIIIConstant;
    float DeConstant = HIIConstant + 0.25*HeIIConstant + 0.5*HeIIIConstant;
    float eUnits = VelocityUnits*VelocityUnits;
    float EUnits = DensityUnits*eUnits;

    for (i=0; i<size; i++) {
      BaryonField[RhoNum][i]   = DensityConstant/DensityUnits;
      BaryonField[TENum][i]    = TEConstant/eUnits;
      BaryonField[V0Num][i]    = VxConstant/VelocityUnits;
      BaryonField[V1Num][i]    = VyConstant/VelocityUnits;
      BaryonField[V2Num][i]    = VzConstant/VelocityUnits;
      BaryonField[EgNum][i]    = EgConstant/EUnits;
      if (NumChemicals > 0) {
	BaryonField[DeNum][i]    = DeConstant/DensityUnits;
	BaryonField[HINum][i]    = HIConstant/DensityUnits;
	BaryonField[HIINum][i]   = HIIConstant/DensityUnits;
      }
      if (NumChemicals == 3) {
	BaryonField[HeINum][i]   = HeIConstant/DensityUnits;
	BaryonField[HeIINum][i]  = HeIIConstant/DensityUnits;
	BaryonField[HeIIINum][i] = HeIIIConstant/DensityUnits;
      }
    }
    if (DualEnergyFormalism)
      for (i=0; i<size; i++)
	BaryonField[IENum][i] = IEConstant/eUnits;
    
    if (debug) {
      fprintf(stdout,"\n  Initializing constant fields using CGS values:\n");
      fprintf(stdout,"        density = %g\n",DensityConstant);
      fprintf(stdout,"   total energy = %g\n",TEConstant);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstant);
      fprintf(stdout,"     x-velocity = %g\n",VxConstant);
      fprintf(stdout,"     y-velocity = %g\n",VyConstant);
      fprintf(stdout,"     z-velocity = %g\n",VzConstant);
      fprintf(stdout,"      radiation = %g\n",EgConstant);
      if (NumChemicals > 0) {
	fprintf(stdout,"      electrons = %g\n",DeConstant);
	fprintf(stdout,"            nHI = %g\n",HIConstant);
	fprintf(stdout,"           nHII = %g\n",HIIConstant);
      }
      if (NumChemicals == 3) {
	fprintf(stdout,"           nHeI = %g\n",HeIConstant);
	fprintf(stdout,"          nHeII = %g\n",HeIIConstant);
	fprintf(stdout,"         nHeIII = %g\n",HeIIIConstant);
      }
      
      fprintf(stdout,"\n  Corresponding Enzo internal values:\n");
      fprintf(stdout,"        density = %g\n",BaryonField[RhoNum][1]);
      fprintf(stdout,"   total energy = %g\n",BaryonField[TENum][1]);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",BaryonField[IENum][1]);
      fprintf(stdout,"     x-velocity = %g\n",BaryonField[V0Num][1]);
      fprintf(stdout,"     y-velocity = %g\n",BaryonField[V1Num][1]);
      fprintf(stdout,"     z-velocity = %g\n",BaryonField[V2Num][1]);
      fprintf(stdout,"      radiation = %g\n",BaryonField[EgNum][1]);
      if (NumChemicals > 0) {
	fprintf(stdout,"      electrons = %g\n",BaryonField[DeNum][1]);
	fprintf(stdout,"            nHI = %g\n",BaryonField[HINum][1]);
	fprintf(stdout,"           nHII = %g\n",BaryonField[HIINum][1]);
      }
      if (NumChemicals == 3) {
	fprintf(stdout,"           nHeI = %g\n",BaryonField[HeINum][1]);
	fprintf(stdout,"          nHeII = %g\n",BaryonField[HeIINum][1]);
	fprintf(stdout,"         nHeIII = %g\n",BaryonField[HeIIINum][1]);
      }
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;

#endif

}

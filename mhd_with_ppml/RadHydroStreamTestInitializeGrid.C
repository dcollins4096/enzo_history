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



int grid::RadHydroStreamTestInitializeGrid(float DensityConstant, 
					   float EgConstant)
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
  FieldType[DeNum = NumberOfBaryonFields++]  = ElectronDensity;
  FieldType[HINum = NumberOfBaryonFields++]  = HIDensity;
  FieldType[HIINum = NumberOfBaryonFields++] = HIIDensity;
  FieldType[HeINum = NumberOfBaryonFields++]   = HeIDensity;
  FieldType[HeIINum = NumberOfBaryonFields++]  = HeIIDensity;    
  FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;    
  
  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;  // no subgrids
 
  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;


  // compute size of fields
  int size = 1;
  for (int dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
 
  // allocate fields
  for (int field=0; field<NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];
 
  // set fluid density, total energy, [internal energy,] velocities, 
  // radiation energy, electron density, chemical species
  int i, j, k;
  float TEConstant, IEConstant, DeConstant, HIConstant, HIIConstant;
  float HeIConstant, HeIIConstant, HeIIIConstant;
  IEConstant = DensityConstant*8.625e7/(Gamma-1.0)*sqrt(sqrt( 2.99792458e10*EgConstant/(4.0*3.1415926535897932*1.8044e-5)));
  TEConstant = IEConstant;
  HIIConstant = 0.0;
  HIConstant = 0.0;
  HeIIConstant = 0.0;
  HeIIIConstant = 0.0;
  HeIConstant = 0.0;
  DeConstant = 0.0;
  for (i=0; i<size; i++) {
    BaryonField[RhoNum][i] = DensityConstant;
    BaryonField[TENum][i]  = TEConstant;
    BaryonField[V0Num][i]  = 0.0;
    BaryonField[V1Num][i]  = 0.0;
    BaryonField[V2Num][i]  = 0.0;
    BaryonField[EgNum][i]  = EgConstant;
    BaryonField[DeNum][i]  = DeConstant;
  }
  if (DualEnergyFormalism)
    for (i=0; i<size; i++)
      BaryonField[IENum][i] = IEConstant;
  for (i=0; i<size; i++) {
      BaryonField[HINum][i]  = HIConstant;
      BaryonField[HIINum][i] = HIIConstant;
  }
  for (i=0; i<size; i++) {
      BaryonField[HeINum][i]   = HeIConstant;
      BaryonField[HeIINum][i]  = HeIIConstant;
      BaryonField[HeIIINum][i] = HeIIIConstant;
  }

  // adjust Radiation Energy Density BC at x0 left edge
  int idx;
  for (k=0; k<GridDimension[2]; k++)
    for (j=0; j<GridDimension[1]; j++)
      for (i=0; i<GridStartIndex[0]+1; i++) {
	idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
	BaryonField[EgNum][idx] = 1.0;
      }

  return SUCCESS;
}

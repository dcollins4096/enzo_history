/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A ZELDOVICH PANCAKE)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

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


int grid::ZeldovichPancakeInitializeGrid(int  ZeldovichPancakeDirection,
				      float ZeldovichPancakeCentralOffset,
				      float ZeldovichPancakeOmegaBaryonNow,
				      float ZeldovichPancakeOmegaCDMNow,
				      float ZeldovichPancakeCollapseRedshift,
				      float ZeldovichPancakeInitialTemperature)
{
  /* declarations */

  float Amplitude, AmplitudeVel, kx, xLagrange, xEulerian, xEulerianOld;
  int   dim, field, i, index;
  const float Pi = 3.14159;

  /* error check */

  if (ZeldovichPancakeDirection < 0 || ZeldovichPancakeDirection >= GridRank) {
    fprintf(stderr, "ZeldovichPancakeDirection is improperly set.\n");
    return FAIL;
  }
  if (ZeldovichPancakeOmegaCDMNow != 0) {
    fprintf(stderr, "Dark matter not yet supported.\n");
    return FAIL;
  }

  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  int vel = NumberOfBaryonFields - 1;
  if (GridRank > 1)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Get the cosmology units so we can convert temperature later. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
        VelocityUnits;
  if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		&TimeUnits, &VelocityUnits, InitialTimeInCodeUnits) == FAIL) {
    fprintf(stderr, "Error in CosmologyGetUnits.\n");
    return FAIL;
  }

  /* Determine the size of the fields. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Allocate space for the fields. */

  for (field = 0; field < NumberOfBaryonFields; field++)
    BaryonField[field] = new float[size];

  /* Find the stride between zones along the pancake direction. */

  int Divisor = 1;
  for (dim = 0; dim < ZeldovichPancakeDirection; dim++)
    Divisor *= GridDimension[dim];

  int NumberOfZones = GridEndIndex[ZeldovichPancakeDirection] -
                      GridStartIndex[ZeldovichPancakeDirection] + 1;

  /* Compute the amplitude of perturbations. */

  Amplitude    = (1+ZeldovichPancakeCollapseRedshift) / (1+InitialRedshift);
  AmplitudeVel = -sqrt(2.0/3.0)*(1.0+ZeldovichPancakeCollapseRedshift) /
                 ((1.0+InitialRedshift) * 2.0*Pi);
  kx           = 2*Pi/float(NumberOfZones);

  /* set density, total energy and velocity in problem dimension */

  for (i = 0; i < size; i++) {

    /* Determine the index along the pancake direction and convert this
       into a number from -0.5NumberOfZones to +0.5NumberOfZones. */

    index     = (i/Divisor % GridDimension[ZeldovichPancakeDirection]) -
                GridStartIndex[ZeldovichPancakeDirection];
    index     = (index+NumberOfZones) % NumberOfZones;
    xLagrange = float(index) + (0.5-ZeldovichPancakeCentralOffset) -
                0.5*float(NumberOfZones);

    /* Convert this Lagrangean position into an Eulerian one using
       a Newton-Raphson method. */

    xEulerian    = xLagrange;
    xEulerianOld = FLOAT_UNDEFINED;
    while (fabs((xEulerian-xEulerianOld)/xEulerian) > TOLERANCE) {
      xEulerianOld = xEulerian;
      xEulerian += (xLagrange - xEulerian + Amplitude*sin(kx*xEulerian)/kx)
	          /(1                     - Amplitude*cos(kx*xEulerian)   );
    }

    /* Set density. */

    BaryonField[0][i] = ZeldovichPancakeOmegaBaryonNow/
                        (1 - Amplitude*cos(kx*xEulerian));

    /* Set total energy, gas energy. */

    BaryonField[1][i] = ZeldovichPancakeInitialTemperature/TemperatureUnits *
          pow(BaryonField[0][i]/ZeldovichPancakeOmegaBaryonNow,Gamma-1) 
	                  / (Gamma-1);
    if (DualEnergyFormalism)
      BaryonField[2][i] = BaryonField[1][i];
    if (HydroMethod != Zeus_Hydro)
      BaryonField[1][i]  += 0.5*pow(AmplitudeVel*sin(kx*xEulerian),2);

    /* Set velocities (some may be set to zero later -- see below). */

    if (HydroMethod == Zeus_Hydro) xEulerian -= 0.5; // only approximate
    BaryonField[vel][i]     = AmplitudeVel * sin(kx*xEulerian);
    if (GridRank > 1)
      BaryonField[vel+1][i] = AmplitudeVel * sin(kx*xEulerian);
    if (GridRank > 2)
      BaryonField[vel+2][i] = AmplitudeVel * sin(kx*xEulerian);
  }

  /* set transverse velocities (i.e. erase any incorrectly set velocities). */

  for (field = vel; field < vel+GridRank; field++)
    if (field != ZeldovichPancakeDirection+vel)
      for (i = 0; i < size; i++)
	BaryonField[field][i] = 0.0;

  return SUCCESS;
}

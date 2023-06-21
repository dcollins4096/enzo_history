/***********************************************************************
/
/  GRID CLASS (FIND INERTIA TENSOR FOR GAS AND DARK MATTER)
/
/  written by: Greg Bryan
/  date:       April, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
#include <stdlib.h>
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
 
/* function prototypes */
 
 
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
 
 
int grid::FindIntertiaTensor(float SphereCenter[MAX_DIMENSION],
			  float SphereRadius,
			  float InertiaTensor[MAX_DIMENSION][MAX_DIMENSION][2],
			  float InertiaTensorWeight[MAX_DIMENSION])
{
 
  int i, j, k, index, dim;
  float radius2, delx, dely, delz;
  const double SolarMass = 1.989e33, Mpc = 3.086e24;
 
  /* Quick check to see if sphere overlaps this grid. */
 
  for (dim = 0; dim < GridRank; dim++)
    if (SphereCenter[dim] - SphereRadius > GridRightEdge[dim] ||
	SphereCenter[dim] + SphereRadius < GridLeftEdge[dim]   )
      return SUCCESS;
 
  /* Find the units if we are using comoving coordinates. */
 
  float DensityConversion = 1, VelocityConversion = 1;
 
  if (ComovingCoordinates) {
 
    float DensityUnits, LengthUnits, VelocityUnits, TimeUnits,TemperatureUnits;
 
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }
 
    /* Convert cgs units to more reasonable values.
       density:   M(solar)/Mpc^3
       velocity:  km/s */
 
    DensityConversion = float(double(DensityUnits) / SolarMass * POW(Mpc, 3));
    VelocityConversion = float(double(VelocityUnits) / 1.0e5);
 
  }
 
  /* Compute cell volume. */
 
  float BoxSize = 1, CellVolume = 1;
  if (ComovingCoordinates)
    BoxSize = ComovingBoxSize/HubbleConstantNow;
  for (dim = 0; dim < GridRank; dim++)
    CellVolume *= CellWidth[dim][0]*BoxSize;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				   Vel3Num, TENum);
 
  /* Loop over grid. */
 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    if (GridRank > 1)
      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - SphereCenter[2];
    else
      delz = 0;
 
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      if (GridRank > 0)
	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - SphereCenter[1];
      else
	dely = 0;
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
 
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	
	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - SphereCenter[0];
 
	radius2 = delx*delx + dely*dely + delz*delz;
 
	if (radius2 <= SphereRadius*SphereRadius)
	  for (dim = 0; dim < GridRank; dim++) {
	    MeanVelocity[dim][0] += BaryonField[DensNum][index]*CellVolume *
	                 BaryonField[Vel1Num+dim][index]*VelocityConversion;
	    MeanVelocityWeight[dim][0] +=
	                            BaryonField[DensNum][index]*CellVolume;
	  }
 
      }
    }
  }
 
  /* Loop over particles. */
 
  for (i = 0; i < NumberOfParticles; i++) {
 
    radius2 = 0;
    for (dim = 0; dim < GridRank; dim++)
      radius2 += (ParticlePosition[dim][i] - SphereCenter[dim])*
	         (ParticlePosition[dim][i] - SphereCenter[dim]);
 
    if (radius2 <= SphereRadius*SphereRadius)
      for (dim = 0; dim < GridRank; dim++) {
	MeanVelocity[dim][1] += ParticleMass[i]*CellVolume *
	                        ParticleVelocity[dim][i]*VelocityConversion;
	MeanVelocityWeight[dim][1] += ParticleMass[i]*CellVolume;
      }
 
  }
 
  return SUCCESS;
}

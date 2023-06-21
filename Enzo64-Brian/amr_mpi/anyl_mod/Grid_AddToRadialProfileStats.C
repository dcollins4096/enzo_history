/***********************************************************************
/
/  GRID CLASS (ADD THE CONTENTS OF THIS GRID TO THE GIVEN RADIAL PROFILE)
/
/  written by: Greg Bryan
/  date:       February, 1996
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
#include "StarParticleData.h"


/* function prototypes */

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, float Time);
int CosmologyComputeExpansionFactor(float time, float *a, float *dadt);
int FindField(int field, int farray[], int numfields);


int grid::AddToRadialProfileStats(float SphereCenter[MAX_DIMENSION], 
			     float SphereRadius,
			     float MeanVelocity[MAX_DIMENSION][3],
			     int NumberOfBins,
			     float ProfileRadius[], 
			     float ProfileValue[][MAX_PROFILES],
			     float ProfileWeight[][MAX_PROFILES],
			     char *ProfileName[MAX_PROFILES],
			     AnalyzeClusterParameters *parameters)
{
  
  int i, j, k, index, m, n, dim;
  float radius2, delx, dely, delz, InertiaWeight, radial_vel, AngWeight,
    velx, vely, velz, vel[MAX_DIMENSION], mu, number_density, 
    meanvelgas,gasdisp,metaldisp,radveldisp,densitydisp,pressuredisp;
  const double SolarMass = 1.989e33, Mpc = 3.086e24, kboltz = 1.38e-16,
               mh = 1.67e-24;

  /* Quick check to see if sphere overlaps this grid. */

  for (dim = 0; dim < GridRank; dim++)
    if (SphereCenter[dim] - SphereRadius > GridRightEdge[dim] ||
	SphereCenter[dim] + SphereRadius < GridLeftEdge[dim]   )
      return SUCCESS;

  /* Find the units if we are using comoving coordinates. */

  float DensityConversion = 1, VelocityConversion = 1, a = 1, dadt;
  float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
        TemperatureUnits = 1;

  if (ComovingCoordinates) {

    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }

    CosmologyComputeExpansionFactor(Time, &a, &dadt);
 
    /* Convert cgs units to more reasonable values.
       density:   M(solar)/Mpc^3 
       velocity:  km/s */

    DensityConversion = float(double(DensityUnits) / SolarMass * POW(Mpc, 3));
    VelocityConversion = float(double(VelocityUnits) / 1.0e5);

  }


  /* Compute cell volume and total grid size. */

  float BoxSize = 1, CellVolume = 1, DomainWidth[MAX_DIMENSION];
  if (ComovingCoordinates) 
    BoxSize = ComovingBoxSize/HubbleConstantNow*a/(1+InitialRedshift);
  int size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    CellVolume *= CellWidth[dim][0]*BoxSize;
    size *= GridDimension[dim];
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
  }
  float DConv = CellVolume*DensityConversion;


  /* Compute the temperature. */

  float *temperature = new float[size];
  this->ComputeTemperatureField(temperature);

  /*  Compute pressure */

  float *pressure = new float[size];
  if (DualEnergyFormalism)
    this->ComputePressureDualEnergyFormalism(Time, pressure);
  else
    this->ComputePressure(Time, pressure);



  /* Set free-free constant, 0.88 is n(e)*n(i) for ionized H/He gas,
     1.15 is a rough value for the mean Gaunt factor, in 10^44 erg/s. */

  float xray_const1 = 0.88 * 1.15 * 1.43e-27 * POW(DensityUnits / mh, 2) *
                     CellVolume * POW(Mpc, 3) * 1e-44;

  /* This is just CellVolume in cm^3 * 1e-44. The 10^-23 accounts for the
     units used in ComputeXrayEmissivity. */

  float xray_const2 = CellVolume * POW(Mpc, 3) * 1e-44 * 1e-23;

  /* Find fields: density, total energy, velocity1-3. */

  if (NumberOfBaryonFields > 0) {
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);
  
  /* Check for metallicty. */

  int MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);


  /* Loop over grid. */

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    if (GridRank > 1) {
      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - SphereCenter[2];
      delz = min(delz, DomainWidth[2]-delz);
    }
    else
      delz = 0;

    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      if (GridRank > 0) {
	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - SphereCenter[1];
	dely = min(dely, DomainWidth[1]-dely);
      }
      else
	dely = 0;
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];

      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	if (BaryonField[NumberOfBaryonFields][index] != 0.0)
	  continue;
	
	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - SphereCenter[0];
	delx = min(delx, DomainWidth[0]-delx);

	radius2 = delx*delx + dely*dely + delz*delz;

	if (radius2 <= SphereRadius*SphereRadius)
	  for (n = 0; n < NumberOfBins; n++)
	    if (radius2 <= ProfileRadius[n+1]*ProfileRadius[n+1]) {

	      float gas_mass = 
		BaryonField[DensNum][index]*CellVolume*DensityConversion;

	      if (gas_mass == 0)
		break;


	      /* get gas info */
	      for (dim = 0; dim < GridRank; dim++)
		vel[dim] = BaryonField[Vel1Num+dim][index];
	      if (HydroMethod == Zeus_Hydro) {
		vel[0] = 0.5*(vel[0] + BaryonField[Vel1Num][index+1]);
		vel[1] = 0.5*(vel[1] + 
			      BaryonField[Vel2Num][index+GridDimension[0]]);
		vel[2] = 0.5*(vel[2] + 
		BaryonField[Vel3Num][index+GridDimension[0]*GridDimension[1]]);
	      }
	      for (dim = 0; dim < GridRank; dim++)
		vel[dim] = VelocityConversion*vel[dim] - MeanVelocity[dim][0];

	      /* 131)  gas velocity dispersion */
	      gasdisp = 0.0;  // set value to zero
	      for (dim = 0; dim < GridRank; dim++)
		gasdisp += vel[dim] * vel[dim];    // this is velocity squared
	      gasdisp = sqrt(gasdisp);  // this is velocity

	      // v-mu
	      gasdisp = gasdisp - (ProfileValue[n][130]/ProfileWeight[n][130]);
	      gasdisp *= gasdisp;  // (v-mu)^2

	      ProfileValue[n][131] += gas_mass * gasdisp;
	      ProfileWeight[n][131] += gas_mass;   

	      if (ProfileName[131] == NULL) 
		ProfileName[131] = "v_sigma_gas (km/s)";

	      /* 132)  gas radial velocity dispersion */

	      radial_vel = (delx*vel[0] + dely*vel[1] + delz*vel[2]) / 
		           sqrt((radius2 == 0)? 1.0 : radius2);

	      // v-mu
	      radveldisp = radial_vel - (ProfileValue[n][6]/ProfileWeight[n][6]);

	      radveldisp *= radveldisp;  // (v-mu)^2

	      ProfileValue[n][132] += gas_mass * radveldisp;
	      ProfileWeight[n][132] += gas_mass;   	      

	      if (ProfileName[132] == NULL) 
		ProfileName[132] = "vr_sigma_gas (km/s)";


	      /* 133) RMS mach number */


	      if (temperature[index] > 20000.0)
		mu = 0.59;
	      else
		mu = 1.22;

	      // calculate sound speed squared in (km/s)^2
	      double sound_speed2 = Gamma*kboltz*temperature[index]/
		(mh*mu) / 1.0e10;
	      
	      double velocity2 = 0.0;
	      for (dim = 0; dim < GridRank; dim++)
		velocity2 += vel[dim] * vel[dim];

	      ProfileValue[n][133] += gas_mass * (sound_speed2/velocity2);
	      ProfileWeight[n][133] += gas_mass;   	      

	      if (ProfileName[133] == NULL) 
		ProfileName[133] = "rms_mach_num";

	      /* 134) Turbulent energy */

	      ProfileValue[n][134] += gas_mass * (velocity2 - radial_vel*radial_vel)*1.0e10;
	      ProfileWeight[n][134] += gas_mass;   	      

	      if (ProfileName[134] == NULL) 
		ProfileName[134] = "turbulent_energy (ergs/g)";


	      /* 135) metallicity fraction dispersion */

	      if (MetalNum != -1) {

		metaldisp = 0.0;

		metaldisp = BaryonField[MetalNum][index]/
		  BaryonField[DensNum][index] -
		  ProfileValue[n][83]/ProfileWeight[n][83];

		ProfileValue[n][135] += metaldisp*metaldisp*gas_mass;

		ProfileWeight[n][135] += gas_mass;

		if (!ProfileName[135]) ProfileName[135] = "Metal dispersion";
	      }


	      /* 137) Density dispersion */
	      densitydisp = BaryonField[DensNum][index]*DensityConversion
		- ProfileValue[n][0]/ProfileWeight[n][0];

	      ProfileValue[n][137] += densitydisp*densitydisp*gas_mass;
	      ProfileWeight[n][137] += gas_mass;
	      if (!ProfileName[137]) ProfileName[137] = "Density dispersion";


	      /* 138) Pressure dispersion */
	      double thispressure;

	      thispressure = double(pressure[index])*double(DensityUnits)
		*double(VelocityUnits*VelocityUnits);

	      pressuredisp = thispressure-ProfileValue[n][136]/ProfileWeight[n][136];

	      ProfileValue[n][138] += pressuredisp*pressuredisp*gas_mass;
	      ProfileWeight[n][138] += gas_mass;
	      if (!ProfileName[138]) ProfileName[138] = "Pressure dispersion";


	      /* Done. */

	      break;
	    }

      } // end i loop
    } // end j loop
  } // end k loop

  } // end: if (NumberOfBaryonFields > 0)


  return SUCCESS;
}
 

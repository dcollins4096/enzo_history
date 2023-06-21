/***********************************************************************
/
/  GRID CLASS (ADD THE CONTENTS OF THIS GRID TO THE GIVEN DISK PROFILE)
/
/  written by: Greg Bryan
/  date:       November, 1999
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
extern void CICDeposit(float *Image, float density, float x, float y, int Size);


int grid::AddToDiskProfileAngMomTrans(float SphereCenter[MAX_DIMENSION], 
			   float SphereRadius,
			   float MeanVelocity[MAX_DIMENSION][3],
			   int NumberOfBins,
			   float ProfileRadius[], 
			   float ProfileValue[][MAX_PROFILES],
			   float ProfileWeight[][MAX_PROFILES],
			   char *ProfileName[MAX_PROFILES],
			   AnalyzeClusterParameters *parameters,
			   float *DiskVector, float *DiskImage[],
			   int DiskImageSize, float DiskRadius)
{
  
  int i, j, k, index, n, dim;
  float radius2, delx, dely, delz, radial_vel,
        vel[MAX_DIMENSION], circ_vel, xpos, ypos, zpos,
        height, dradius2, length, gas_dens, delta_circ_vel, delta_radial_vel;
  const double SolarMass = 1.989e33, Mpc = 3.086e24, mh = 1.67e-24;

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

  fprintf(stderr,"in Grid_AddToDiskProfileAngMomTrans(1)\n");

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

  /* Set free-free constant, 0.88 is n(e)*n(i) for ionized H/He gas,
     1.15 is a rough value for the mean Gaunt factor, in 10^44 erg/s. */

  float xray_const = 0.88 * 1.15 * 1.43e-27 * POW(DensityUnits / mh, 2) *
                     CellVolume * POW(Mpc, 3) * 1e-44;

  /* Find fields: density, total energy, velocity1-3. */

  if (NumberOfBaryonFields > 0) {
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);
  
  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, 
      H2IINum, DINum, DIINum, HDINum;
  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
		     HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      return FAIL;
    }

  /* Compute the temperature. */

  float *temperature = new float[size];
  this->ComputeTemperatureField(temperature);

  /* Compute two vectors which are perpindicular to the disk vector. */

  float theta, phi, axis0[3], axis1[3], tang[3], plane[3], pi = 3.14159;
  theta = atan(sqrt(DiskVector[0]*DiskVector[0] + DiskVector[1]*DiskVector[1])
	       /DiskVector[2]);  // from -pi/2 to +pi/2
  phi   = atan2(DiskVector[1],DiskVector[0]);  // from -pi to +pi

  axis0[0] = sin(theta+0.5*pi)*cos(phi);
  axis0[1] = sin(theta+0.5*pi)*sin(phi);
  axis0[2] = cos(theta+0.5*pi);

  axis1[0] =  (DiskVector[2]*axis0[1] - DiskVector[1]*axis0[2]);
  axis1[1] = -(DiskVector[2]*axis0[0] - DiskVector[0]*axis0[2]);
  axis1[2] =  (DiskVector[1]*axis0[0] - DiskVector[0]*axis0[1]);

  /*
  printf("DV = %g %g %g\n", DiskVector[0], DiskVector[1], DiskVector[2]);
  printf("a0 = %g %g %g\n", axis0[0], axis0[1], axis0[2]);
  printf("a1 = %g %g %g\n", axis1[0], axis1[1], axis1[2]);
  printf("%g %g %g\n", 
	 DiskVector[0]*axis0[0]+DiskVector[1]*axis0[1]+DiskVector[2]*axis0[2],
	 DiskVector[0]*axis1[0]+DiskVector[1]*axis1[1]+DiskVector[2]*axis1[2],
	 axis1[0]*axis0[0]+axis1[1]*axis0[1]+axis1[2]*axis0[2]);
  */

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

	/* If the subgrid field is set, then ignore this cell. */

	if (BaryonField[NumberOfBaryonFields][index] != 0.0)
	  continue;

	float gas_mass = 
	  BaryonField[DensNum][index]*CellVolume*DensityConversion;

	/*
	if(i==j && j==k) fprintf(stderr,"mass:  %e %d %d %d\n",gas_mass,i,j,k);

	if(i==j && j==k) fprintf(stderr,"%e %e %e %e\n",
				 gas_mass/CellVolume, parameters->LowerDensityCutoff,
				 temperature[index],  parameters->ColdTemperatureCutoff);
	*/

	/* if the cell is not dense and cold, then ignore. */

	if (gas_mass/CellVolume < parameters->LowerDensityCutoff ||
	    temperature[index] > parameters->ColdTemperatureCutoff)
	  continue;
	
	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - SphereCenter[0];
	delx = min(delx, DomainWidth[0]-delx);

	radius2 = delx*delx + dely*dely + delz*delz;
	if (radius2 > SphereRadius*SphereRadius)
	  continue;

	/* Compute radius in plane. */

	height = DiskVector[0]*delx + DiskVector[1]*dely + DiskVector[2]*delz;
	dradius2 = radius2 - height*height;

	if (dradius2 <= SphereRadius*SphereRadius)
	  for (n = 0; n < NumberOfBins; n++)
	    if (dradius2 <= ProfileRadius[n+1]*ProfileRadius[n+1]) {

	      if (gas_mass == 0)
		break;

	      /* calculate radial velocity and delta_radial velocity of this gas packet */
	      radial_vel = (delx*vel[0] + dely*vel[1] + delz*vel[2]) / 
		           sqrt((radius2 == 0)? 1.0 : radius2);

	      /* this is delta_radial velocity -- this cell's v_r minus the mean v_r */
	      delta_radial_vel = radial_vel - (ProfileValue[n][104]/ProfileWeight[n][104]);

	      /* calculate circular velocity  First, compute radial vector
		 in disk, then compute tangential vector in disk, then
		 normalize it and dot with velocity vector. */

              plane[0] = delx - height*DiskVector[0];  // radial vector in disk
              plane[1] = dely - height*DiskVector[1];
              plane[2] = delz - height*DiskVector[2];
	      tang[0] =  (DiskVector[2]*plane[1] - DiskVector[1]*plane[2]);
	      tang[1] = -(DiskVector[2]*plane[0] - DiskVector[0]*plane[2]);
	      tang[2] =  (DiskVector[1]*plane[0] - DiskVector[0]*plane[1]);
	      length = sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]
			    +tiny_number);
	      circ_vel = (tang[0]*vel[0] + tang[1]*vel[1] + tang[2]*vel[2])/
		length;

	      /* delta_circular_velocity -- v_circ - <v_circ> */
	      
	      delta_circ_vel = circ_vel - (ProfileValue[n][105]/ProfileWeight[n][105]);

	      /* 106) gas radial velocity. times circular velocity * cell mass */
	      ProfileValue[n][106] += radial_vel * circ_vel * gas_mass;
	      if (ProfileName[106] == NULL) ProfileName[106] = "vr_vcirc_gas (km/s)";

	      /* 107) delta gas radial velocity times delta circular velocity * cell mass  */
	      ProfileValue[n][107] += delta_radial_vel * delta_circ_vel * gas_mass;
	      if (ProfileName[107] == NULL) ProfileName[107] = "dvr_dvcirc_gas (km/s)";

	      /* 108) delta gas radial velocity times circular velocity * cell mass */
	      ProfileValue[n][108] += delta_radial_vel * circ_vel * gas_mass;
	      if (ProfileName[108] == NULL) ProfileName[108] = "dvr_vcirc_gas (km/s)";

	      /* 109) gas radial velocity times delta circular velocity * cell mass */
	      ProfileValue[n][109] += radial_vel * delta_circ_vel * gas_mass;
	      if (ProfileName[109] == NULL) ProfileName[109] = "vr_dvcirc_gas (km/s)";

	      /* Done. */

	      break;
	    }

      } // end i loop
    } // end j loop
  } // end k loop

  /* Clean up. */

  delete temperature;

  } // end: if (NumberOfBaryonFields > 0)


  return SUCCESS;
}

/* 
void CICDeposit(float *Image, float density, float x, float y, int Size)
{
  int i, j, ip1, jp1;
  float dx, dy;

  i = int(x + 0.5) - 1;
  j = int(y + 0.5) - 1;

  if (i < 0 || j < 0 || i >= Size || j >= Size)
    return;

  dx = float(i) + 1.5 - x;
  dy = float(j) + 1.5 - y;

  ip1 = min(i+1, Size-1);
  jp1 = min(j+1, Size-1);

  //fprintf(stderr,"in CICDeposit: %e %e %d %d %d %d\n",dx,dy,i,j,ip1,jp1);

  Image[i  +j  *Size] += density*     dx *     dy;
  Image[ip1+j  *Size] += density*(1.0-dx)*     dy;
  Image[i  +jp1*Size] += density*     dx *(1.0-dy);
  Image[ip1+jp1*Size] += density*(1.0-dx)*(1.0-dy);

}
*/

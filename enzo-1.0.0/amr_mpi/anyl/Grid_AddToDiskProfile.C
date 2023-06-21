/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
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
		      float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
void CICDeposit(float *Image, float density, float x, float y, int Size);


int grid::AddToDiskProfile(float SphereCenter[MAX_DIMENSION], 
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
        height, dradius2, length, gas_dens;
  const double SolarMass = 1.989e33, Mpc = 3.086e24, mh = 1.67e-24;

  /* Quick check to see if sphere overlaps this grid. */

  for (dim = 0; dim < GridRank; dim++)
    if (SphereCenter[dim] - SphereRadius > GridRightEdge[dim] ||
	SphereCenter[dim] + SphereRadius < GridLeftEdge[dim]   )
      return SUCCESS;

  /* Find the units if we are using comoving coordinates. */

  FLOAT a = 1, dadt;
  float DensityConversion = 1, VelocityConversion = 1;
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

#if 0
  printf("DV = %g %g %g\n", DiskVector[0], DiskVector[1], DiskVector[2]);
  printf("a0 = %g %g %g\n", axis0[0], axis0[1], axis0[2]);
  printf("a1 = %g %g %g\n", axis1[0], axis1[1], axis1[2]);
  printf("%g %g %g\n", 
	 DiskVector[0]*axis0[0]+DiskVector[1]*axis0[1]+DiskVector[2]*axis0[2],
	 DiskVector[0]*axis1[0]+DiskVector[1]*axis1[1]+DiskVector[2]*axis1[2],
	 axis1[0]*axis0[0]+axis1[1]*axis0[1]+axis1[2]*axis0[2]);
#endif

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

	      /* 100) gas density (in Msolar/Mpc^3). */

	      ProfileValue[n][100] += gas_mass;
	      ProfileWeight[n][100] += CellVolume;
	      if (ProfileName[100] == NULL) 
		ProfileName[100] = "d_gas (Ms/Mpc^3)";

	      /* 101) gas rms velocity (in km/s, mass weighted). 
                    (first set vel[] to velocity, correctly adjusted for face-
	             centering if applicable). */

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

	      for (dim = 0; dim < GridRank; dim++)
		ProfileValue[n][101] += gas_mass * vel[dim] * vel[dim];
	      ProfileWeight[n][101] += gas_mass;
	      if (ProfileName[101] == NULL) 
		ProfileName[101] = "v_rms_gas_3d (km/s)";

	      /* 102) gas temperature (in K).  Mass-weighted in cube is
		    less than 40 Mpc (i.e. not a cluster sim), otherwise
	            weight by x-ray luminosity. */

	      ProfileValue[n][102] += temperature[index] * gas_mass;
	      ProfileWeight[n][102] += gas_mass;
	      if (ProfileName[102] == NULL) 
		ProfileName[102] = "temp_gas_mass (K)";
	      
	      /* 103) number of samples. */

	      ProfileValue[n][103] += 1.0;
	      ProfileWeight[n][103] = 0;
	      if (ProfileName[103] == NULL) ProfileName[103] = "N_gas";

	      /* 104) gas radial velocity. */

	      radial_vel = (delx*vel[0] + dely*vel[1] + delz*vel[2]) / 
		           sqrt((radius2 == 0)? 1.0 : radius2);
	      ProfileValue[n][104] += radial_vel*gas_mass;
	      ProfileWeight[n][104] += gas_mass;
	      if (ProfileName[104] == NULL) ProfileName[104] = "vr_gas (km/s)";

	      /* 105) gas circular velocity. First, compute radial vector
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
	      ProfileValue[n][105] += circ_vel*gas_mass;
	      ProfileWeight[n][105] += gas_mass;
	      if (ProfileName[105] == NULL) ProfileName[105] = "vcirc (km/s)";

	      /* Add to image. */

	      /* project against the two in-disk axis vectors to get
		 index in image. */

	      xpos = delx*axis0[0] + dely*axis0[1] + delz*axis0[2];
	      ypos = delx*axis1[0] + dely*axis1[1] + delz*axis1[2];

	      xpos = (xpos/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;
	      ypos = (ypos/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;
	      zpos = (height/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;
	      
	      gas_dens = gas_mass/
		POW(BoxSize*DiskRadius*SphereRadius/DiskImageSize, 2);

	      if (min(xpos, min(ypos, zpos)) > 0 && 
		  max(xpos, max(ypos, zpos)) < DiskImageSize) {
		CICDeposit(DiskImage[0], gas_dens, xpos, ypos, DiskImageSize);
		CICDeposit(DiskImage[1], gas_dens, xpos, zpos, DiskImageSize);
		CICDeposit(DiskImage[2], gas_dens, ypos, zpos, DiskImageSize);
	      }

	      /* Done. */

	      break;
	    }

      } // end i loop
    } // end j loop
  } // end k loop

  /* Clean up. */

  delete temperature;

  } // end: if (NumberOfBaryonFields > 0)

#ifdef UNUSED

  float velx, vely, velz, InertiaWeight, AngWeight;

  /* Loop over particles. */

  dely = delz = 0;
 
  for (i = 0; i < NumberOfParticles; i++) {

                      delx = ParticlePosition[0][i] - SphereCenter[0];
    if (GridRank > 1) dely = ParticlePosition[1][i] - SphereCenter[1];
    if (GridRank > 2) delz = ParticlePosition[2][i] - SphereCenter[2];

    radius2 = delx*delx + dely*dely + delz*delz;

    if (radius2 <= SphereRadius*SphereRadius)
      for (n = 0; n < NumberOfBins; n++)
	if (radius2 <= ProfileRadius[n+1]*ProfileRadius[n+1]) {

	  float part_mass = ParticleMass[i]*DensityConversion*CellVolume;

	  if (part_mass == 0)
	    break;

	  /* a) Dark matter particles */

	  if (ParticleNumber[i] < STAR_PARTICLE_NUMBER_START) {

	    /* 30) dm density (in Msolar/Mpc^3). */

	    ProfileValue[n][30] += part_mass;
	    ProfileWeight[n][30] -= CellVolume;
	    if (ProfileName[30] == NULL) ProfileName[30] = "d_dm (Ms/Mpc^3)";

	    /* 31) dm rms velocity (in km/s, mass weighted). */

	    for (dim = 0; dim < GridRank; dim++)
	      ProfileValue[n][31] += 
		(ParticleVelocity[dim][i]*VelocityConversion - 
		 MeanVelocity[dim][1])*
		(ParticleVelocity[dim][i]*VelocityConversion - 
		 MeanVelocity[dim][1])*
		part_mass;
	    ProfileWeight[n][31] += part_mass;
	    if (ProfileName[31] == NULL) 
	      ProfileName[31] = "v_rms_dm_3d (km/s)";

	    /* 32) number of samples. */

	    ProfileValue[n][32] += 1.0;
	    ProfileWeight[n][32] = 0;
	    if (ProfileName[32] == NULL) ProfileName[32] = "N_dm";

	    /* 33) dm radial velocity. */

	   velx = ParticleVelocity[0][i]*VelocityConversion-MeanVelocity[0][1];
	   vely = ParticleVelocity[1][i]*VelocityConversion-MeanVelocity[1][1];
	   velz = ParticleVelocity[2][i]*VelocityConversion-MeanVelocity[2][1];
	    radial_vel = (delx*velx + dely*vely + delz*velz) / 
	               sqrt(radius2 + 1.0e-20);
	    ProfileValue[n][33] += radial_vel*part_mass;
	    ProfileWeight[n][33] += part_mass;
	    if (ProfileName[33] == NULL) ProfileName[33] = "vr_dm (km/s)";

	    /* 34) dm radial velocity dispersion. */

	    ProfileValue[n][34] += radial_vel*radial_vel*part_mass;
	    ProfileWeight[n][34] += part_mass;
	    if (ProfileName[34] == NULL) ProfileName[34] = "vr_rms_dm (km/s)";

	    /* 35-37) angular momentum vector. */

	    AngWeight = part_mass * BoxSize;
	    ProfileValue[n][35] += AngWeight * ( vely*delz - velz*dely);
	    ProfileValue[n][36] += AngWeight * (-velx*delz + velz*delx);
	    ProfileValue[n][37] += AngWeight * ( velx*dely - vely*delx);
	    for (m = 35; m < 38; m++)
	      ProfileWeight[n][m] += part_mass;
	    ProfileName[35] = "Lx_dm (km/s*Mpc)";
	    ProfileName[36] = "Ly_dm (km/s*Mpc)";
	    ProfileName[37] = "Lz_dm (km/s*Mpc)";

	    /* 50-55) inertia tensor. */

	    InertiaWeight = part_mass * BoxSize * BoxSize;
	    ProfileValue[n][50] += delx * delx * InertiaWeight;
	    ProfileValue[n][51] += delx * dely * InertiaWeight;
	    ProfileValue[n][52] += delx * delz * InertiaWeight;
	    ProfileValue[n][53] += dely * dely * InertiaWeight;
	    ProfileValue[n][54] += dely * delz * InertiaWeight;
	    ProfileValue[n][55] += delz * delz * InertiaWeight;
	    for (m = 50; m < 56; m++)
	      ProfileWeight[n][m] += part_mass;
	    if (ProfileName[50] == NULL) {
	      ProfileName[50] = "I_xx_dm";
	      ProfileName[51] = "I_xy_dm";
	      ProfileName[52] = "I_xz_dm";
	      ProfileName[53] = "I_yy_dm";
	      ProfileName[54] = "I_yz_dm";
	      ProfileName[55] = "I_zz_dm";
	    }

	  } else {

	  /* b) star particles */

	    /* 60) sp density (in Msolar/Mpc^3). */

	    ProfileValue[n][60] += part_mass;
	    ProfileWeight[n][60] -= CellVolume;
	    if (ProfileName[60] == NULL) ProfileName[60] = "d_sp (Ms/Mpc^3)";

	    /* 61) sp rms velocity (in km/s, mass weighted). */

	    for (dim = 0; dim < GridRank; dim++)
	      ProfileValue[n][61] += 
		(ParticleVelocity[dim][i]*VelocityConversion - 
		 MeanVelocity[dim][1])*
		(ParticleVelocity[dim][i]*VelocityConversion - 
		 MeanVelocity[dim][1])*
		part_mass;
	    ProfileWeight[n][61] += part_mass;
	    if (ProfileName[61] == NULL) 
	      ProfileName[31] = "v_rms_sp_3d (km/s)";

	    /* 62) number of samples. */

	    ProfileValue[n][62] += 1.0;
	    ProfileWeight[n][62] = 0;
	    if (ProfileName[62] == NULL) ProfileName[62] = "N_sp";

	    /* 63) sp radial velocity. */

	   velx = ParticleVelocity[0][i]*VelocityConversion-MeanVelocity[0][1];
	   vely = ParticleVelocity[1][i]*VelocityConversion-MeanVelocity[1][1];
	   velz = ParticleVelocity[2][i]*VelocityConversion-MeanVelocity[2][1];
	    radial_vel = (delx*velx + dely*vely + delz*velz) / 
	               sqrt(radius2 + 1.0e-20);
	    ProfileValue[n][63] += radial_vel*part_mass;
	    ProfileWeight[n][63] += part_mass;
	    if (ProfileName[63] == NULL) ProfileName[63] = "vr_sp (km/s)";

	    /* 64) sp radial velocity dispersion. */

	    ProfileValue[n][64] += radial_vel*radial_vel*part_mass;
	    ProfileWeight[n][64] += part_mass;
	    if (ProfileName[64] == NULL) ProfileName[64] = "vr_rms_sp (km/s)";

	    /* 65-67) angular momentum vector. */

	    AngWeight = part_mass * BoxSize;
	    ProfileValue[n][65] += AngWeight * ( vely*delz - velz*dely);
	    ProfileValue[n][66] += AngWeight * (-velx*delz + velz*delx);
	    ProfileValue[n][67] += AngWeight * ( velx*dely - vely*delx);
	    for (m = 65; m < 68; m++)
	      ProfileWeight[n][m] += part_mass;
	    ProfileName[65] = "Lx_sp (km/s*Mpc)";
	    ProfileName[66] = "Ly_sp (km/s*Mpc)";
	    ProfileName[67] = "Lz_sp (km/s*Mpc)";

	  } // end if (star particle)

	  /* Done. */

	  break; 
	}  // end of radius loop
  } // end: loop over particles

#endif /* UNUSED */

  return SUCCESS;
}
 
void CICDeposit(float *Image, float density, float x, float y, int Size)
{
  int i, j, ip1, jp1;
  float dx, dy;

  i = int(x + 0.5) - 1;
  j = int(y + 0.5) - 1;

  if (i < 0 || j < 0 || i >= Size || j >= Size)
    return;

  dx = i + 1.5 - x;
  dy = j + 1.5 - y;

  ip1 = min(i+1, Size-1);
  jp1 = min(j+1, Size-1);

  Image[i  +j  *Size] += density*     dx *     dy;
  Image[ip1+j  *Size] += density*(1.0-dx)*     dy;
  Image[i  +jp1*Size] += density*     dx *(1.0-dy);
  Image[ip1+jp1*Size] += density*(1.0-dx)*(1.0-dy);

}

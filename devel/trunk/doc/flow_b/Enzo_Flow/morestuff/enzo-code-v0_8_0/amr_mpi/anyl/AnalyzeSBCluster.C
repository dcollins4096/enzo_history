/***********************************************************************
/
/  ANALYZE THE SANTA BARBARA CLUSTER
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "AnalyzeClusters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#undef DEFINE_STORAGE

/* function prototypes */

int InitializeNew(  char *filename, HierarchyEntry &TopGrid, TopGridData &tgd, 
		    ExternalBoundary &Exterior);
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
int EvolveHierarchy(                HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, LevelHierarchyEntry *Array[]);
void ExtractSection(HierarchyEntry &TopGrid, TopGridData &tgd, 
		    LevelHierarchyEntry *Array[], ExternalBoundary *Exterior,
		    int ExtractStart[], int ExtractEnd[], int ExtractLevel);
int ProjectToPlane(TopGridData &MetaData, LevelHierarchyEntry *LevelArray[],
		   int ProjectStart[], int ProjectEnd[], 
		   float ProjectStartCoordinates[], 
		   float ProjectEndCoordinates[], int ProjectLevel,
		   int ProjectionDimension, char *ProjectionFileName);
int InterpretCommandLine(int argc, char *argv[], char *myname,
			 int &restart, int &debug, int &extract,
			 int &project, int &ProjectionDimension,
			 char *RestartFile[], char *ParameterFile[],
			 int RegionStart[], int RegionEnd[], 
			 float RegionStartCoordinates[],
			 float RegionEndCoordinates[], int &Level);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);

main(int argc, char *argv[])
{

  /* Main declarations */

  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

  int level;

  /* Initialize */

  debug                = TRUE;
//  GFAlreadyInitialized = FALSE;
  char *myname         = argv[0];
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;

  /* Error check */

  if (argc != 2) {
    fprintf(stderr, "usage: %s amr_saved_filename\n", myname);
    exit(EXIT_FAILURE);
  }

  /* Read the saved file. */

  if (ReadAllData(argv[1], &TopGrid, MetaData, &Exterior) == FAIL) {
    fprintf(stderr, "Error in ParameterFile %s.\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

  /* ------------------------------------------------------------ */
  /* general declarations. */

  int i, j, profile;
  float pi = 3.14159;
  LevelHierarchyEntry *Temp2, *Temp;

  float BoxSize = 1;
  if (ComovingCoordinates) 
    BoxSize = ComovingBoxSize/HubbleConstantNow;

#define MAX_BINS 16
  float InnerEdge = 0.01/BoxSize, OuterEdge = 10.0/BoxSize;   // in Mpc;
  float MeanVelocityOuterEdge = 2.77/BoxSize; // in Mpc (fix by generalizing!!)

  /* ------------------------------------------------------------ */
  /* Find the highest density spot. */

  float MaxPosition[MAX_DIMENSION], MaxDensity = -huge_number;
  if (debug) printf("->finding maximum\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->FindMaximumBaryonDensity(MaxPosition, &MaxDensity);
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* ------------------------------------------------------------ */
  /* Zero the solution (on this grid) which is underneath any subgrid
     (so we get only the high resolution solution from the subgrid). */

  if (debug) printf("->zeroing redundant solution\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp2 = LevelArray[level+1];
      while (Temp2 != NULL) {
	Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						 ZERO_ALL_FIELDS);
	Temp2 = Temp2->NextGridThisLevel;
      }
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* ------------------------------------------------------------ */
  /* Set up radial grid. */

  float ProfileValue[MAX_BINS][MAX_PROFILES], ProfileRadius[MAX_BINS+1],
        ProfileWeight[MAX_BINS][MAX_PROFILES];
  char  *ProfileName[MAX_PROFILES];

  for (i = 0; i < MAX_PROFILES; i++) {
    for (j = 0; j < MAX_BINS; j++) {
      ProfileValue[j][i] = 0;
      ProfileWeight[j][i] = 0;
    }
    ProfileName[i] = NULL;
  }

  /* Compute radius (ProfileRadius is inner edge of bin). */

  ProfileRadius[0] = 0;
  float dlogRadius = (log10(OuterEdge) - log10(InnerEdge)) /
    float(MAX_BINS-1);
  for (i = 1; i < MAX_BINS+1; i++)
    ProfileRadius[i] = POW(10, log10(InnerEdge) + float(i-1)*dlogRadius);

  /* ------------------------------------------------------------ */
  /* Calculate the mean velocity. */

  float MeanVelocity[MAX_DIMENSION][3], MeanVelocityWeight[MAX_DIMENSION][3];
  for (i = 0; i < MAX_DIMENSION; i++) {
    MeanVelocity[i][0] = 0;
    MeanVelocity[i][1] = 0;
    MeanVelocityWeight[i][0] = 0;
    MeanVelocityWeight[i][1] = 0;
  }

  if (debug) printf("->finding mean velocity\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->FindMeanVelocity(MaxPosition, MeanVelocityOuterEdge, 
				       MeanVelocity, MeanVelocityWeight);
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* Compute mean velocities within MeanVelocityOuterEdge for gas [0], dm [1],
     and total [2]. */

  for (i = 0; i < MAX_DIMENSION; i++) {
    if ((MeanVelocityWeight[i][0] + MeanVelocityWeight[i][1]) > 0)
      MeanVelocity[i][2] = (MeanVelocity[i][0] + MeanVelocity[i][1])/
	(MeanVelocityWeight[i][0] + MeanVelocityWeight[i][1]);
    if (MeanVelocityWeight[i][0] > 0)
      MeanVelocity[i][0] /= MeanVelocityWeight[i][0];
    if (MeanVelocityWeight[i][1] > 0)
      MeanVelocity[i][1] /= MeanVelocityWeight[i][1];
  }

  /* ------------------------------------------------------------ */
  /* Calculate profiles. */

  if (debug) printf("->finding profiles\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->AddToRadialProfile(MaxPosition, OuterEdge, MeanVelocity,
					 MAX_PROFILES, MAX_BINS,
					 ProfileRadius,
					 ProfileValue, ProfileWeight,
					 ProfileName);
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* If the ProfileWeight is < 0, it is to be replaced by the effective
     volume in the annulus. */

  for (i = 0; i < MAX_PROFILES; i++)
    for (j = 0; j < MAX_BINS; j++)
      if (ProfileWeight[j][i] < 0)
	ProfileWeight[j][i] = POW(BoxSize, 3) * 4.0*pi/3.0 *
	  (POW(ProfileRadius[j+1], 3) - POW(ProfileRadius[j], 3));

  /* Normalize profile values by dividing each by it's cumulative weight. */

  for (i = 0; i < MAX_PROFILES; i++)
    for (j = 0; j < MAX_BINS; j++) {

      if (ProfileWeight[j][i] > 0)
	ProfileValue[j][i] /= ProfileWeight[j][i];

      /* If this is a radial velocity dispersion, subtract the mean velocity 
	 which MUST be the profile immediately before. */

      if (ProfileName[i] != NULL)
	if (strstr(ProfileName[i], "vr_rms") != NULL)
	  ProfileValue[j][i] -= POW(ProfileValue[j][i-1],2);

      /* If this is a velocity rms profile, convert from 3D rms^2 to 3D rms. */

      if (ProfileName[i] != NULL)
	if (strstr(ProfileName[i], "_rms") != NULL)
	  ProfileValue[j][i] = sqrt(ProfileValue[j][i]/1.0);
    }
  
  /* ------------------------------------------------------------ */
  /* Compute total density and find r200 (d200 is overdensity of 200 in
     M(solar)/Mpc^3). */

  if (debug) printf("->computing r200\n");
  float r200 = 0, 
        average_dens = 2.78e11*OmegaMatterNow*POW(HubbleConstantNow, 2);

  for (j = 0; j < MAX_BINS; j++) {

    /* Compute the mass, in M(solar), within this annalus. */

//    ProfileValue[j][18] = (ProfileValue[j][0] + ProfileValue[j][10]) *
//      (POW(ProfileRadius[j+1], 3) - POW(ProfileRadius[j], 3) ) * 
//      POW(BoxSize, 3) * 4.0*pi/3.0;
    ProfileValue[j][18] = ProfileValue[j][ 0]*ProfileWeight[j][ 0] +
                          ProfileValue[j][10]*ProfileWeight[j][10];
					 

    /* Keep a running sum of the total mass within radius j+1. */

    if (j > 0)
      ProfileValue[j][18] += ProfileValue[j-1][18];

    /* Compute the average overdensity within the sphere with radius [j+1]. */

    ProfileValue[j][19] = ProfileValue[j][18] / 
      (POW(ProfileRadius[j+1]*BoxSize, 3) * 4.0*pi/3.0) / average_dens;
  }

  ProfileName[18] = "d_total (Ms)";
  ProfileName[19] = "od_total";

  /* Find the radius at which the cumulative overdensity goes through 200. */

  for (j = 0; j < MAX_BINS-1; j++)
    if (ProfileValue[j+1][19] <= 200.0) {
      r200 = ProfileRadius[j+1] + (ProfileRadius[j+2] - ProfileRadius[j+1]) *
	(200.0                 - ProfileValue[j][19]) /
	(ProfileValue[j+1][19] - ProfileValue[j][19] + tiny_number);
      break;
    }
  if (r200 > OuterEdge)
    r200 = OuterEdge;

  /* ------------------------------------------------------------ */
  /* Compute properties averaged over r200. */

  float R200Radius[2], R200Value[1][MAX_PROFILES], 
        R200Weight[1][MAX_PROFILES];

  /* Allocate and clear r200 data. */

  R200Radius[0] = 0;
  R200Radius[1] = r200;
  for (i = 0; i < MAX_PROFILES; i++) {
    R200Value[0][i] = 0;
    R200Weight[0][i] = 0;
  }

  if (r200 > 0) {

    /* Find values within r200. */

    if (debug) printf("->finding r200 average\n");
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	Temp->GridData->AddToRadialProfile(MaxPosition, r200, MeanVelocity,
					   MAX_PROFILES, 1,
					   R200Radius, R200Value, R200Weight,
					   ProfileName);
	Temp = Temp->NextGridThisLevel;
      }
    }

    /* If the R200Weight is < 0, it is to be replaced by the annulus volume. */

    for (i = 0; i < MAX_PROFILES; i++)
      if (R200Weight[0][i] < 0) 
	R200Weight[0][i] = POW(BoxSize, 3) * 4.0*pi/3.0 *
	  (POW(R200Radius[1], 3) - POW(R200Radius[0], 3));

    /* Normalize profile values by dividing each by it's cumulative weight. */

    for (i = 0; i < MAX_PROFILES; i++) {

      if (R200Weight[0][i] > 0) 
	R200Value[0][i] /= R200Weight[0][i];

      /* If this is a radial velocity dispersion, subtract the mean velocity 
	 which MUST be the profile immediately before. */

      if (ProfileName[i] != NULL)
	if (strstr(ProfileName[i], "vr_rms") != NULL)
	  R200Value[j][i] -= POW(R200Value[j][i-1], 2);

      /* If this is a velocity rms profile, convert from 3D rms^2 to 3D rms. */

      if (ProfileName[i] != NULL)
	if (strstr(ProfileName[i], "_rms") != NULL)
	  R200Value[0][i] = sqrt(R200Value[0][i]/1.0);

    }

  } // end: if (r200 > 0)
  
  /* ------------------------------------------------------------ */
  /* Open the output file. */

  FILE *fptr = fopen("AnalyzeSBCluster.out", "w");

  /* Output global values. */

  fprintf(fptr, "# MaxBaryonPosition  = %g %g %g\n", MaxPosition[0], 
	  MaxPosition[1], MaxPosition[2]);
  fprintf(fptr, "# MaxBaryonValue     = %g\n", MaxDensity);
  fprintf(fptr, "# MeanVelocity (gas) = %g %g %g (km/s)\n", 
	  MeanVelocity[0][0], MeanVelocity[1][0], MeanVelocity[2][0]);
  fprintf(fptr, "# MeanVelocity (dm ) = %g %g %g (km/s)\n", 
	  MeanVelocity[0][1], MeanVelocity[1][1], MeanVelocity[2][1]);
  fprintf(fptr, "# MeanVelocity (tot) = %g %g %g (km/s)\n", 
	  MeanVelocity[0][2], MeanVelocity[1][2], MeanVelocity[2][2]);
  fprintf(fptr, "#  (Within r = %g Mpc\n#\n", MeanVelocityOuterEdge*BoxSize);
  fprintf(fptr, "# InnerEdge = %g Mpc  OuterEdge = %g Mpc\n", 
	  InnerEdge*BoxSize, OuterEdge*BoxSize);

  /* Output r200 information. */
  
  fprintf(fptr, "# r200              = %g (Mpc)\n", r200*BoxSize);
  for (profile = 0; profile < MAX_PROFILES; profile++)
    if (ProfileName[profile] != NULL)
      fprintf(fptr, "%12.5g ", R200Value[0][profile]);
  fprintf(fptr, "\n");  
  fprintf(fptr, "#\n");

  /* Print column info. */

  fprintf(fptr, "#\n# COLUMN   DESCRIPTION\n");
  int ColumnNumber = 1;
  for (profile = 0; profile < MAX_PROFILES; profile++)
    if (ProfileName[profile] != NULL)
      fprintf(fptr, "#   %3d    %s\n", ColumnNumber++, ProfileName[profile]);

  /* Print header. */
  
  fprintf(fptr, "#\n#%12.12s ", "r_cent (Mpc)");
  for (profile = 0; profile < MAX_PROFILES; profile++)
    if (ProfileName[profile] != NULL)
      fprintf(fptr, "%12.12s ", ProfileName[profile]);
  fprintf(fptr, "%12.12s %12.12s %12.12s ", "vol_gas (real)", "vol (calc)",
	  "vol_dm (real)");
  fprintf(fptr, "\n");

  /* Loop over all radial bins, printing each. */

  for (j = 0; j < MAX_BINS; j++) {

    fprintf(fptr, "%12.5g ",0.5*(ProfileRadius[j]+ProfileRadius[j+1])*BoxSize);

    for (profile = 0; profile < MAX_PROFILES; profile++)
      if (ProfileName[profile] != NULL)
	if (ProfileValue[j][profile] == 0)
	  fprintf(fptr, "%12.5g ", tiny_number);
	else
	  fprintf(fptr, "%12.5g ", ProfileValue[j][profile]);

    fprintf(fptr, "%12.5g %12.5g %12.5g ", ProfileWeight[j][0], 
	    (POW(ProfileRadius[j+1], 3) - POW(ProfileRadius[j], 3) ) * 
	    POW(BoxSize, 3) * 4.0*pi/3.0,
	    ProfileWeight[j][10]);
    fprintf(fptr, "\n");
  }

  fclose(fptr);

  exit(EXIT_SUCCESS);
}

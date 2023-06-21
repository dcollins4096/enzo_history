/***********************************************************************
/
/  COMPUTE A CLUSTER PROFILE
/
/  written by: Greg Bryan
/  date:       August, 1997
/  modified1:  Robert Harkness, July 2002
/  modified2:  James Bordner, June 2003,  merged HDF{4|5}_AnalyzeCluster.C
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
#include <assert.h>

#ifdef USE_HDF4
#include <df.h>
#endif
#ifdef USE_HDF5
#include "hdf4.h"
#include <hdf5.h>
#endif

#include "macros_and_parameters.h"
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

#define MAX_BINS 200

/* function prototypes */

int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int AnalyzeClusterReadParameterFile(char *filename, int &NumberOfCenters,
				    float *CenterList[],
				    AnalyzeClusterParameters *parameters);
int AnalyzeClusterComputeClumpingFactor(LevelHierarchyEntry *LevelArray[],
					TopGridData *MetaData,
					int NumberOfGridPoints,
					int NumberOfParticles,
					float SphereCenter[], 
					float SphereRadius, char *Name);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int CommunicationInitialize(int *argc, char **argv[]);
int CommunicationFinalize();
void my_exit(int status);

// HDF5 function prototypes

#include "extern_hdf5.h"

main(int argc, char *argv[])
{

#ifdef USE_HDF5
  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id, attr_type_id;

  hsize_t     OutDims[3];
  hsize_t     Dimm;
  hsize_t     Slab_Dims[4];
  hsize_t     Slab_Rank;

  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[4], file_count[4], file_block[4];

  hssize_t    mem_offset;
  hssize_t    file_offset[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int Part;
  int io_log = 1;
#endif

  CommunicationInitialize(&argc, &argv);

  /* Main declarations */

  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

  /* Initialize */

  debug                = TRUE;
  char *myname         = argv[0];
  AnalyzeClusterParameters parameters;

  /* general declarations. */

  int i, j, profile, dim, ret, level;
  float pi = 3.14159;
  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.67e-8;
  char Name[MAX_LINE_LENGTH];
  char DiskImageName[MAX_LINE_LENGTH];
#ifdef USE_HDF5
  char DiskDataName[MAX_LINE_LENGTH];
#endif
  LevelHierarchyEntry *Temp2, *Temp;

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;

  /* Error check */

  if (argc != 3) {
    fprintf(stderr, "usage: %s amr_file anyl_parameter_file\n", myname);
    my_exit(EXIT_FAILURE);
  }

  /* Read the saved file. */

  SetDefaultGlobalValues(MetaData);
  if (ReadAllData(argv[1], &TopGrid, MetaData, &Exterior) == FAIL) {
    fprintf(stderr, "Error in ParameterFile %s.\n", argv[1]);
    my_exit(EXIT_FAILURE);
  }
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

  /* Set top grid boundary conditions. */

  if (TopGrid.GridData->SetExternalBoundaryValues(&Exterior) == FAIL) {
    fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
    my_exit(EXIT_FAILURE);
  }

  /* Read the anyl parameter file. */

  int NumberOfCenters;
  float *CenterList[MAX_DIMENSION]; 

  AnalyzeClusterReadParameterFile(argv[2], NumberOfCenters, CenterList,
				  &parameters);
  int NumberOfPoints = parameters.npoints;

  /* From the time, compute the current redshift. */

  FLOAT a = 1, dadt;
  float CurrentRedshift = 0.0, OmegaCurvatureNow = 0.0, Esquared;

  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(LevelArray[0]->GridData->ReturnTime(),
					&a, &dadt) == FAIL) {
      fprintf(stderr, "Error in ComputeExpansionFactor.\n");
      my_exit(EXIT_FAILURE);
    }

    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;

    /* If unset, then set virial_dens according to spherical top-hot collapse
       (see Bryan & Norman 1998). */

    OmegaCurvatureNow = 1.0 - OmegaMatterNow - OmegaLambdaNow;
    Esquared = OmegaMatterNow    * POW(1.0+CurrentRedshift, 3) +
               OmegaCurvatureNow * POW(1.0+CurrentRedshift, 2) +
               OmegaLambdaNow;
    float x = OmegaMatterNow * POW(1.0+CurrentRedshift, 3) / Esquared - 1.0;
    if (parameters.virial_dens < 0) {
      if (OmegaCurvatureNow < 1.0e-4)
	parameters.virial_dens = 18.0*pi*pi + 82.0*x - 39.0*x*x;
      else if (OmegaLambdaNow < 1.0e-4)
	parameters.virial_dens = 18.0*pi*pi + 60.0*x - 32.0*x*x;
      else {
	printf("This cosmology not supported.\n");
	my_exit(EXIT_FAILURE);
      }
    }

  } // end: if (ComovingCoordinates)

  /* Set BoxSize. */

  float BoxSize = 1;
  if (ComovingCoordinates) 
    BoxSize = ComovingBoxSize/HubbleConstantNow;
  float InnerEdge = parameters.rinner/BoxSize, 
        OuterEdge = parameters.router/BoxSize;
  BoxSize *= a/(1+InitialRedshift);

  /* ------------------------------------------------------------ */
  /* Zero the solution (on this grid) which is underneath any subgrid
     (so we get only the high resolution solution from the subgrid). */

  if (debug) printf("->zeroing redundant solution & setting BCs\n");

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
      if (level > 0) {
	Temp->GridHierarchyEntry->ParentGrid->GridData->SetTime(
				Temp->GridData->ReturnTime());
	Temp->GridData->InterpolateBoundaryFromParent
	                   (Temp->GridHierarchyEntry->ParentGrid->GridData);
	CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, level);
      }
      Temp2 = LevelArray[level+1];
      while (Temp2 != NULL) {
	Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						 ZERO_UNDER_SUBGRID_FIELD);
	Temp2 = Temp2->NextGridThisLevel;
      }
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* ------------------------------------------------------------ */
  /* Loop over centers. */

  float Center[MAX_DIMENSION];

  for (int center = 0; center < NumberOfCenters; center++) {

    if (NumberOfCenters > 1) {
      debug = FALSE;
      printf("Computing center %d\n", center);
    }

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      Center[dim] = CenterList[dim][center];

    /* Set base name. */

    if (NumberOfCenters == 1)
      sprintf(Name, "%s", "AnalyzeCluster");
    else
      sprintf(Name, "%s%.3d", "AnalyzeCluster", center);

    /* ------------------------------------------------------------ */
    /* Find the highest density spot if Center not specified. */

    float MaxDensity = -huge_number;

    if (Center[0] == FLOAT_UNDEFINED) {
      if (debug) printf("->finding maximum\n");
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	Temp = LevelArray[level];
	while (Temp != NULL) {
	  Temp->GridData->FindMaximumBaryonDensity(Center, &MaxDensity);
	  Temp = Temp->NextGridThisLevel;
	}
      }
    }

    /* ------------------------------------------------------------ */
    /* Set up radial grid. */

    float ProfileValue[MAX_BINS][MAX_PROFILES], ProfileRadius[MAX_BINS+1],
          ProfileWeight[MAX_BINS][MAX_PROFILES];
    char  *ProfileName[MAX_PROFILES];

    for (i = 0; i < MAX_PROFILES; i++) {
      for (j = 0; j < MAX_BINS; j++) {
	ProfileValue[j][i] = 0.0;
	ProfileWeight[j][i] = 0.0;
      }
      ProfileName[i] = NULL;
    }

    /* Compute radii (ProfileRadius is inner edge of bin). */

    ProfileRadius[0] = 0;
    float dlogRadius = (log10(OuterEdge) - log10(InnerEdge)) /
      float(NumberOfPoints-1);
    for (i = 1; i < NumberOfPoints+1; i++)
      ProfileRadius[i] = POW(10, log10(InnerEdge) + float(i-1)*dlogRadius);

    /* ------------------------------------------------------------ */
    /* Calculate rvirial. */
    /* Compute total density and find rvir (dvir is overdensity of virial_dens
       M(solar)/Mpc^3). */

    /* First, compute profiles. */

    if (debug) printf("->computing r_virial (now with critical_dens)\n");
    float MeanVelocity[MAX_DIMENSION][3], MeanVelocityWeight[MAX_DIMENSION][3];

    for (i = 0; i < MAX_DIMENSION; i++) {
      MeanVelocity[i][0] = 0;
      MeanVelocity[i][1] = 0;
      MeanVelocityWeight[i][0] = 0;
      MeanVelocityWeight[i][1] = 0;
    }

    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->AddToRadialProfile(Center, OuterEdge, MeanVelocity,
					       NumberOfPoints, ProfileRadius,
					       ProfileValue, ProfileWeight,
					       ProfileName, &parameters) 
	    == FAIL) {
	  fprintf(stderr, "Error in grid->AddToRadialProfile.\n");
	  exit(EXIT_FAILURE);
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

    /* Compute rvir based on critical density (as of Feb 21/2000). */

    float rvir = 0, mvir = 0, mvir_gas = 0, spin_gas = 0, spin_dm = 0,
          mvir_star = 0,
       average_dens = 2.78e11*OmegaMatterNow*POW(HubbleConstantNow, 2) *
	  POW((1+InitialRedshift)/a, 3),
       critical_dens = 2.78e11*POW(HubbleConstantNow, 2) * Esquared;
    
    for (j = 0; j < NumberOfPoints; j++) {

      /* Compute the mass (gas+dm+star), in M(solar), within this annulus. */

      ProfileValue[j][38] = ProfileValue[j][0] + ProfileValue[j][30] +
	                    ProfileValue[j][60];

      /* Keep a running sum of the total mass within radius j+1. */

      if (j > 0)
	ProfileValue[j][38] += ProfileValue[j-1][38];

      /* Compute the average overdensity within the sphere with radius [j+1].*/

      ProfileValue[j][39] = ProfileValue[j][38] / 
	(POW(ProfileRadius[j+1]*BoxSize, 3) * 4.0*pi/3.0) / critical_dens;
    }

    /* Find the radius at which the cumulative overdensity goes through 
       virial_dens (typically 200). */

    int j1 = 1;

    while (ProfileValue[j1][39] < 1e-19 && j1 < NumberOfPoints)
      j1++;

    for (j = j1; j < NumberOfPoints; j++)
      if (ProfileValue[j][39] <= parameters.virial_dens) {
	rvir = ProfileRadius[j] + (ProfileRadius[j+1] - ProfileRadius[j]) *
	  (parameters.virial_dens - ProfileValue[j-1][39]) /
	  (ProfileValue[j][39] - ProfileValue[j-1][39] + tiny_number);
	break;
      }
    if (rvir > OuterEdge || rvir == 0) {
      printf("warning: rvir (%g) > OuterEdge (%g)\n", rvir, OuterEdge);
      rvir = OuterEdge;
    }
    j1 = j;
    mvir = POW(rvir*BoxSize, 3)*4.0/3.0*pi*critical_dens*
           parameters.virial_dens;

    /* ------------------------------------------------------------ */
    /* Calculate the mean velocity. */

    for (i = 0; i < MAX_DIMENSION; i++) {
      MeanVelocity[i][0] = 0;
      MeanVelocity[i][1] = 0;
      MeanVelocityWeight[i][0] = 0;
      MeanVelocityWeight[i][1] = 0;
    }
    float MeanVelocityOuterEdge = parameters.MeanVelocityVirialFraction*rvir;

    if (debug) printf("->finding mean velocity\n");
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	Temp->GridData->FindMeanVelocity(Center, MeanVelocityOuterEdge,
					 MeanVelocity, MeanVelocityWeight);
	Temp = Temp->NextGridThisLevel;
      }
    }

    /* Compute mean velocities within MeanVelocityOuterEdge for gas [0], 
       dm [1], and total [2]. */

    for (i = 0; i < MAX_DIMENSION; i++) {
      if ((MeanVelocityWeight[i][0] + MeanVelocityWeight[i][1]) > 0)
	MeanVelocity[i][2] = (MeanVelocity[i][0] + MeanVelocity[i][1])/
	  (MeanVelocityWeight[i][0] + MeanVelocityWeight[i][1]);
      if (MeanVelocityWeight[i][0] > 0)
	MeanVelocity[i][0] /= MeanVelocityWeight[i][0];
      if (MeanVelocityWeight[i][1] > 0)
	MeanVelocity[i][1] /= MeanVelocityWeight[i][1];
    }

    /* Set Tvirial (in K) and (if required) the cold temperature cutoff
       (from Bryan & Norman, 1998). */

    float mu = 1.22;    // Set the mean mass per particle either to fully
    if (mvir > 1.0e13)  //   ionized or fully neutral depending on (arbitrary)
      mu = 0.59;        //    cutoff.
    float Tvir = parameters.VirialTemperatureNormalization * 2.73e7 * mu *
                 POW(mvir/1.0e15, 2.0/3.0) * 
                 POW(HubbleConstantNow*HubbleConstantNow*
		     parameters.virial_dens * Esquared, 1.0/3.0);
    if (parameters.ColdTemperatureCutoffVirialFraction > 0)
      parameters.ColdTemperatureCutoff = 
	parameters.ColdTemperatureCutoffVirialFraction * Tvir;

    /* ------------------------------------------------------------ */
    /* Calculate profiles. */

    for (i = 0; i < MAX_PROFILES; i++) {
      for (j = 0; j < MAX_BINS; j++) {
	ProfileValue[j][i] = 0;
	ProfileWeight[j][i] = 0;
      }
      ProfileName[i] = NULL;
    }

    if (debug) printf("->finding profiles\n");
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	Temp->GridData->AddToRadialProfile(Center, OuterEdge, MeanVelocity,
					   NumberOfPoints, ProfileRadius,
					   ProfileValue, ProfileWeight,
					   ProfileName, &parameters);
	Temp = Temp->NextGridThisLevel;
      }
    }

    /* Compute the cummulative mass and overdensity. */    

    for (j = 0; j < NumberOfPoints; j++) {

      /* Compute the mass (gas+dm+star), in M(solar), within this annalus. */

      ProfileValue[j][38] = ProfileValue[j][0] + ProfileValue[j][30] +
	                    ProfileValue[j][60];

      /* Keep a running sum of the total mass within radius j+1. */

      ProfileValue[j][29] = ProfileValue[max(j-1,0)][29] + ProfileValue[j][0];
      ProfileValue[j][57] = ProfileValue[max(j-1,0)][57] + ProfileValue[j][30];
      ProfileValue[j][68] = ProfileValue[max(j-1,0)][68] + ProfileValue[j][60];
      if (j > 0)
	ProfileValue[j][38] += ProfileValue[j-1][38];

      /* Compute the average overdensity within the sphere with radius [j+1].*/
      
      ProfileValue[j][39] = ProfileValue[j][38] / 
	(POW(ProfileRadius[j+1]*BoxSize, 3) * 4.0*pi/3.0) / critical_dens;

      /* Compute the dynamical time within the sphere with radius j+1. */

      ProfileValue[j][22] = sqrt(3*pi/(16*GravConst*
				       max(ProfileValue[j][39],1e-20)*
				       critical_dens*SolarMass/(Mpc*Mpc*Mpc)));
    }

    ProfileName[22] = "T_dyn (s)";
    ProfileName[39] = "od_total";
    ProfileName[29] = "m_gas (Ms)";
    ProfileName[38] = "m_total (Ms)";
    ProfileName[57] = "m_dm (Ms)";
    if (ProfileName[60] != NULL)
      ProfileName[68] = "m_star (Ms)";
    
    /* If the ProfileWeight is < 0, it is to be replaced by the effective
       volume in the annulus. */

    for (i = 0; i < MAX_PROFILES; i++)
      for (j = 0; j < NumberOfPoints; j++)
	if (ProfileWeight[j][i] < 0)
	  ProfileWeight[j][i] = POW(BoxSize, 3) * 4.0*pi/3.0 *
	    (POW(ProfileRadius[j+1], 3) - POW(ProfileRadius[j], 3));

    /* Normalize profile values by dividing each by it's cumulative weight. */

    for (i = 0; i < MAX_PROFILES; i++)
      for (j = 0; j < NumberOfPoints; j++) {

	if (ProfileWeight[j][i] > 0.0){
	  ProfileValue[j][i] /= ProfileWeight[j][i];
	  /*        fprintf(stderr,"profile is greater than zero! %i\n ",i);*/

	}

	/* If this is a radial velocity dispersion, subtract the mean velocity 
	   which MUST be the profile immediately before. */

	if (ProfileName[i] != NULL)
	  if (strstr(ProfileName[i], "vr_rms") != NULL)
	    ProfileValue[j][i] -= POW(ProfileValue[j][i-1],2);

	/* If this is a velocity rms profile, convert from 3D ms to 3D rms. */

	if (ProfileName[i] != NULL)
	  if (strstr(ProfileName[i], "_rms") != NULL)
	    ProfileValue[j][i] = sqrt(max(fabs(ProfileValue[j][i]),0.0)/1.0);
      }

    /* Compute the dark matter relaxation time due to the finite
       number of particles (from (8-71) in Binney & Tremaine). */
    
    for (j = 0; j < NumberOfPoints; j++)
      if (ProfileValue[j][32] > 1) {
	ProfileValue[j][56] = 0.34 * POW(ProfileValue[j][31]*1e5, 3) /
	  (POW(GravConst, 2) * 
	   (ProfileValue[j][30]*ProfileWeight[j][30] / 
	    ProfileValue[j][32]) * SolarMass *
	   ProfileValue[j][30] * SolarMass / (Mpc*Mpc*Mpc) *
	   log(0.4*ProfileValue[j][32]) );
	ProfileName[56] = "T_relax (s)";
      }

    /* Compute the cumulative luminosity and luminosity-weighted
       temperatures. */

    float CumLum = 0, CumTemp = 0, CumClump = 0;
    for (j = 0; j < NumberOfPoints; j++) {
      CumLum += ProfileValue[j][3];
      CumTemp += ProfileValue[j][84]*ProfileValue[j][3];
      ProfileValue[j][85] = CumLum;
      ProfileValue[j][86] = CumTemp/CumLum;
    }
    ProfileName[85] = "xray_cumulative_luminosity (erg/s)";
    ProfileName[86] = "temp_gas_cumulative_xray_weighted (K)";

    /* Compute the clumping factor and the cumulative clumping factor. */

    for (j = 0; j < NumberOfPoints; j++) {
      if (ProfileValue[j][0] > 0)
	ProfileValue[j][87] /= ProfileValue[j][0]*ProfileValue[j][0];
      CumClump += ProfileValue[j][87]*(ProfileValue[j][29] - 
	       ( (j==0) ? 0 : ProfileValue[j-1][29]));
      ProfileValue[j][88] = CumClump/ProfileValue[j][29];
    }
    ProfileName[88] = "cumulative clumping factor";

    /* ------------------------------------------------------------ */
    /* Compute properties averaged over rvir. */

    float RvirRadius[2], RvirValue[1][MAX_PROFILES], 
          RvirWeight[1][MAX_PROFILES];

    /* Allocate and clear rvir data. */

    RvirRadius[0] = 0;
    RvirRadius[1] = rvir;
    for (i = 0; i < MAX_PROFILES; i++) {
      RvirValue[0][i] = 0;
      RvirWeight[0][i] = 0;
    }

    if (rvir > 0) {

      /* Find values within rvir. */

      if (debug) printf("->finding rvir average\n");
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	LevelHierarchyEntry *Temp = LevelArray[level];
	while (Temp != NULL) {
	  Temp->GridData->AddToRadialProfile(Center, rvir, MeanVelocity, 1,
					     RvirRadius, RvirValue, RvirWeight,
					     ProfileName, &parameters);
	  Temp = Temp->NextGridThisLevel;
	}
      }

      /* Set gas, dm and star masses within rvir. */

      mvir_gas = RvirValue[0][29] = RvirValue[0][0];
      RvirValue[0][57] = RvirValue[0][30];
      mvir_star = RvirValue[0][60] = RvirValue[0][60];

      /* If the RvirWeight is < 0, it is replaced by the annulus volume. */

      for (i = 0; i < MAX_PROFILES; i++)
	if (RvirWeight[0][i] < 0) 
	  RvirWeight[0][i] = POW(BoxSize, 3) * 4.0*pi/3.0 *
	    (POW(RvirRadius[1], 3) - POW(RvirRadius[0], 3));

      /* Normalize profile values by dividing each by cumulative weight. */

      for (i = 0; i < MAX_PROFILES; i++) {

	if (RvirWeight[0][i] > 0) 
	  RvirValue[0][i] /= RvirWeight[0][i];

	/* If this is a radial velocity dispersion, subtract the mean velocity 
	   which MUST be the profile immediately before. */

	if (ProfileName[i] != NULL)
	  if (strstr(ProfileName[i], "vr_rms") != NULL) {
	    RvirValue[0][i] -= POW(RvirValue[0][i-1], 2);
	  }

	/* If this is an rms profile, convert from 3D rms^2 to 3D rms. */

	if (ProfileName[i] != NULL)
	  if (strstr(ProfileName[i], "_rms") != NULL)
	    RvirValue[0][i] = sqrt(max(RvirValue[0][i],0.0)/1.0);

      }

      /* Compute the spin parameter.  Note: this is usually defined as 
	 L E^1/2 / GM^5/2, but here is l e^1/2 / GM  where l (e) is the
	 specific angular momentum (specific energy).*/

      float ang_mom, SpinUnits = Mpc * 1.0e5 * 1.0e5 /
	                         (GravConst * SolarMass);
      ang_mom = sqrt(RvirValue[0][17]*RvirValue[0][17] +
		     RvirValue[0][18]*RvirValue[0][18] +
		     RvirValue[0][19]*RvirValue[0][19] );
      spin_gas = SpinUnits * ang_mom * RvirValue[0][1] / mvir;
      ang_mom = sqrt(RvirValue[0][35]*RvirValue[0][35] +
		     RvirValue[0][36]*RvirValue[0][36] +
		     RvirValue[0][37]*RvirValue[0][37] );
      spin_dm = SpinUnits * ang_mom * RvirValue[0][31] / mvir;

      /* Compute clumping factor. */

      RvirValue[0][87] /= RvirValue[0][0]*RvirValue[0][0];

    } // end: if (rvir > 0)

    /* ------------------------------------------------------------ */
    /* If requested, compute special clumping factor. */

    if (parameters.ComputeClumpingFactor)
      if (AnalyzeClusterComputeClumpingFactor(LevelArray, &MetaData,
		 nint(RvirValue[0][4]), nint(RvirValue[0][32]), 
                 Center, rvir, Name) == FAIL) {
	fprintf(stderr, "Error in AnalyzeClusterComputeClumpingFactor\n");
	my_exit(EXIT_FAILURE);
      }

    /* ------------------------------------------------------------ */
    /* If requested, compute disk properties. */

    if (rvir > 0 && parameters.ComputeDiskInformation) {

      /* Allocate space. */

      int const nimages = 3;
      int image_size = POW(parameters.DiskImageSize, 2);
      float *DiskImage[nimages];

      for (j = 0; j < nimages; j++) {
	DiskImage[j] = new float[image_size];
	for (i = 0; i < image_size; i++)
	  DiskImage[j][i] = 0;
      }

      /* Compute unit disk vector. */

      float DiskVector[3], length;
      length = sqrt(RvirValue[0][80]*RvirValue[0][80] +
		    RvirValue[0][81]*RvirValue[0][81] +
		    RvirValue[0][82]*RvirValue[0][82]);

      for (dim = 0; dim < 3; dim++)
	DiskVector[dim] = RvirValue[0][80+dim]/length;

      /* Get info from grids. */

      if (debug) printf("->finding disk info\n");

      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	LevelHierarchyEntry *Temp = LevelArray[level];
	while (Temp != NULL) {
	  Temp->GridData->AddToDiskProfile(Center, rvir, MeanVelocity,
				       NumberOfPoints,
				       ProfileRadius,
				       ProfileValue, ProfileWeight,
				       ProfileName, &parameters,
				       DiskVector, DiskImage,
				       parameters.DiskImageSize, 
				       parameters.DiskRadius);
	  Temp = Temp->NextGridThisLevel;
	}
      }

      /* If the RvirWeight is < 0, it is replaced by the annulus volume. */

      for (i = 100; i < 110; i++)
	for (j = 0; j < NumberOfPoints; j++)
	  if (ProfileWeight[j][i] < 0)
	    ProfileWeight[j][i] = POW(BoxSize, 2) * pi *
	      (POW(ProfileRadius[j+1], 3) - POW(ProfileRadius[j], 2));

      /* Normalize profile values by dividing each by cumulative weight. */

      for (i = 100; i < 110; i++)

	for (j = 0; j < NumberOfPoints; j++)
	  if (ProfileWeight[j][i] > 0) 
	    ProfileValue[j][i] /= ProfileWeight[j][i];

      /* Save disk image. */

#ifdef USE_HDF4
      int32 OutDims[2];
      OutDims[0] = OutDims[1] = parameters.DiskImageSize;
      if (DFSDsetdims(2, OutDims) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDsetdims.\n");
	return FAIL;
      }
#endif

#ifdef USE_HDF5
      OutDims[0] = nimages;
      OutDims[1] = parameters.DiskImageSize;
      OutDims[2] = parameters.DiskImageSize;

      Slab_Rank = 3;
      Slab_Dims[0] = OutDims[0];
      Slab_Dims[1] = OutDims[1];
      Slab_Dims[2] = OutDims[2];

      Dimm = image_size;  // OutDims[1] * OutDims[2]
#endif

      strcpy(DiskImageName, Name);
      strcat(DiskImageName, ".DiskImage");
#ifdef USE_HDF5
      strcpy(DiskDataName, "DiskImage");

      int ii = sizeof(float);

      switch(ii)
      {

        case 4:
          mem_type_id = HDF5_R4;
          file_type_id = HDF5_FILE_R4;
          break;

        case 8:
          mem_type_id = HDF5_R8;
          file_type_id = HDF5_FILE_R8;
          break;

        default:
          mem_type_id = HDF5_R4;
          file_type_id = HDF5_FILE_R4;

      }

// Data in memory is considered 1D, stride 1, with zero offset

      mem_stride = 1;      // contiguous elements
      mem_count = Dimm;    // number of elements in field
      mem_offset = 0;      // zero offset in buffer
      mem_block = 1;       // single element blocks

// 1D memory model

      mem_dsp_id = H5Screate_simple(1, &Dimm, NULL);

      if (io_log) printf("H5Screate mem_dsp_id %d\n",mem_dsp_id);
      assert( mem_dsp_id != h5_error );

      h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);

      if (io_log) printf("H5Sselect mem slab: %d\n",h5_status);
      assert( h5_status != h5_error );

// Each image indexed by first index of Slab_Dims
// Number of images = nimages
// Size of each image is Dimm

      file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);

      if (io_log) printf("H5Screate file_dsp_id %d\n",file_dsp_id);
      assert( file_dsp_id != h5_error );



      if (io_log) printf("Calling H5Fcreate with Name = %s\n",DiskImageName);

      file_id = H5Fcreate(DiskImageName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      if (io_log) printf("H5Fcreate File id: %d\n",file_id);
      assert( file_id != h5_error );

      if (io_log) printf("Calling H5Dcreate with Name = %s\n",Name);

      dset_id =  H5Dcreate(file_id, DiskDataName, file_type_id, file_dsp_id, H5P_DEFAULT);

      if (io_log) printf("H5Dcreate Dataset id: %d\n",dset_id);
      assert( dset_id != h5_error );
#endif /* USE_HDF5 */

      float32 *float_temp = new float32[image_size];

      for (j = 0; j < nimages; j++) {

#ifdef USE_HDF5
      Part = j;
#endif

	for (i = 0; i < image_size; i++)
	  float_temp[i] = float32(DiskImage[j][i]);
#ifdef USE_HDF4
	if (j == 0)
	  ret = DFSDputdata(DiskImageName, 2, OutDims, (VOIDP) float_temp);
	if (j != 0)
	  ret = DFSDadddata(DiskImageName, 2, OutDims, (VOIDP) float_temp);

	if (ret == HDF_FAIL) {
	  fprintf(stderr, "Error in DFSDput/adddata.\n");
	  return FAIL;
	}
#endif

#ifdef USE_HDF5

// Offset[0] is the component Part of Npart components.  Data for each
// Part are contiguous in the file, so stride = 1.

        file_stride[0] = 1;      // contiguous elements
        file_count[0] = 1;       // one component per call
        file_offset[0] = Part;   // component Part of Npart
        file_block[0] = 1;       // single element blocks

        for ( dim = 1; dim < Slab_Rank; dim++ )
        {
          file_stride[dim] = 1;                   // contiguous elements
          file_count[dim] = OutDims[dim];         // field dimensions
          file_offset[dim] = 0;                   // complete field, no offset
          file_block[dim] = 1;                    // single element blocks
        }

        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);

        if (io_log) printf("H5Sselect file slab: %d\n",h5_status);
        assert( h5_status != h5_error );

        h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) float_temp);

        if (io_log) printf("H5Dwrite %d\n",h5_status);
        assert( h5_status != h5_error );
#endif /* USE_HDF5 */

      }

#ifdef USE_HDF5
      h5_status = H5Sclose(mem_dsp_id);

      if (io_log) printf("H5Sclose %d\n",h5_status);
      assert( h5_status != h5_error );

      h5_status = H5Dclose(dset_id);

      if (io_log) printf("H5Dclose %d\n",h5_status);
      assert( h5_status != h5_error );

      h5_status = H5Fclose(file_id);

      if (io_log) printf("H5Fclose %d\n",h5_status);
      assert( h5_status != h5_error );
#endif /* USE_HDF5 */

      /* Clean up. */

      delete float_temp;
      for (j = 0; j < nimages; j++)
	delete DiskImage[j];
      
    }
  
    /* ------------------------------------------------------------ */
    /* Output the data: naming convention:
       file 0: Name
       file 1: Name.Inertial
       file 2: Name.Species
       file 3: Name.DarkMatter
       file 4: Name.StarParticles
       file 5: Name.Disk            */

#define NUMBER_OF_FILES 6

    FILE *fptrs[NUMBER_OF_FILES];
    char *FileName[NUMBER_OF_FILES];
    int   ColumnNumber[NUMBER_OF_FILES];

    /* Generate output file names. */

    for (i = 0; i < NUMBER_OF_FILES; i++) {
      ColumnNumber[i] = 3;
      FileName[i] = new char[MAX_LINE_LENGTH];
      strcpy(FileName[i], Name);
    }
    strcat(FileName[1], ".Inertial");
    strcat(FileName[2], ".Species");
    strcat(FileName[3], ".DarkMatter");
    strcat(FileName[4], ".StarParticles");
    strcat(FileName[5], ".Disk");

    /* Open output files. */

    for (i = 0; i < NUMBER_OF_FILES; i++) {
      fptrs[i] = NULL;
      if ((i != 4 || StarParticleCreation > 0) &&
	  (i != 5 || parameters.ComputeDiskInformation == TRUE))
	fptrs[i] = fopen(FileName[i], "w");
    }

    /* Generate an array of file pointers so each profile knows which
       file it should be sent to. */

    int ProfileFile[MAX_PROFILES];
    for (profile = 0; profile < 8; profile++)
      ProfileFile[profile] = 0;
    for (profile = 8; profile < 17; profile++)
      ProfileFile[profile] = 2;
    for (profile = 17; profile < 30; profile++)
      ProfileFile[profile] = 0;
    for (profile = 25; profile < 28; profile++)
      ProfileFile[profile] = 2;
    for (profile = 40; profile < 46; profile++)
      ProfileFile[profile] = 1;
    for (profile = 30; profile < 40; profile++)
      ProfileFile[profile] = 3;
    for (profile = 50; profile < 56; profile++)
      ProfileFile[profile] = 1;
    ProfileFile[56] = ProfileFile[57] = 3;
    if (StarParticleCreation > 0)
      for (profile = 60; profile < 69; profile++)
	ProfileFile[profile] = 4;
    for (profile = 80; profile < 89; profile++)
      ProfileFile[profile] = 0;
    ProfileFile[83] = 2;
    for (profile = 100; profile < 110; profile++)
      ProfileFile[profile] = 5;
    for (profile = 120; profile < 130; profile++)
      ProfileFile[profile] = 2;

    /* Output global values. */

    fprintf(fptrs[0], "# Center             = %.10g %.10g %.10g\n", Center[0], 
	  Center[1], Center[2]);
    fprintf(fptrs[0], "# MaxBaryonValue     = %g\n", MaxDensity);
    fprintf(fptrs[0], "# MeanVelocity (gas) = %g %g %g (km/s)\n", 
	  MeanVelocity[0][0], MeanVelocity[1][0], MeanVelocity[2][0]);
    fprintf(fptrs[0], "# MeanVelocity (dm ) = %g %g %g (km/s)\n", 
	  MeanVelocity[0][1], MeanVelocity[1][1], MeanVelocity[2][1]);
    fprintf(fptrs[0], "# MeanVelocity (tot) = %g %g %g (km/s)\n", 
	  MeanVelocity[0][2], MeanVelocity[1][2], MeanVelocity[2][2]);
    fprintf(fptrs[0], "#  (Within r = %g Mpc)\n#\n",
	    MeanVelocityOuterEdge*BoxSize);
    fprintf(fptrs[0], "# L (gas)            = %g %g %g (Mpc km/s)\n",
	    RvirValue[0][17], RvirValue[0][18], RvirValue[0][19]);
    fprintf(fptrs[0], "# L (dm)             = %g %g %g (Mpc km/s)\n",
	    RvirValue[0][35], RvirValue[0][36], RvirValue[0][37]);
    fprintf(fptrs[0], "# InnerEdge = %g Mpc  OuterEdge = %g Mpc\n\n", 
	  InnerEdge*BoxSize, OuterEdge*BoxSize);

    /* Output general information about simulation. */

    fprintf(fptrs[0], "# OmegaMatterNow     = %g\n", OmegaMatterNow);
    fprintf(fptrs[0], "# OmegaLambdaNow     = %g\n", OmegaLambdaNow);
    fprintf(fptrs[0], "# HubbleConstantNow  = %g\n", HubbleConstantNow);
    fprintf(fptrs[0], "# CurrentRedshift    = %g\n", CurrentRedshift);
    fprintf(fptrs[0], "# AMR file name      = %s\n\n", argv[1]);

    /* Output rvir information. */

    fprintf(fptrs[0], "# rvir               = %g (Mpc)\n", rvir*BoxSize);
    fprintf(fptrs[0], "# mvir               = %g (M_solar)\n", mvir);
    fprintf(fptrs[0], "# VirialDensity      = %g (to critical density)\n", 
	    parameters.virial_dens);
    fprintf(fptrs[0], "# mvir (gas,dm,star) = %g %g %g (M_solar)\n", 
	    mvir_gas, mvir - mvir_gas - mvir_star, mvir_star);
    fprintf(fptrs[0], "# spin (gas, dm)     = (%g, %g)\n", spin_gas, spin_dm);
    fprintf(fptrs[0], "# Tvir predicted (K) = %g\n", Tvir);
    if (parameters.XrayTableFileName != NULL) {
      fprintf(fptrs[0], "# XrayLowerCutoffkeV = %g\n", 
	      parameters.XrayLowerCutoffkeV);
      fprintf(fptrs[0], "# XrayUpperCutoffkeV = %g\n", 
	      parameters.XrayUpperCutoffkeV);
      fprintf(fptrs[0], "# XrayTableFileName  = %s\n", 
	      parameters.XrayTableFileName);
    }
    fprintf(fptrs[0], "# ColdTempCutoff (K) = %g\n\n", 
	    parameters.ColdTemperatureCutoff);

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fprintf(fptrs[i], "#");

    for (profile = 0; profile < MAX_PROFILES; profile++)
      if (ProfileName[profile] != NULL)
	fprintf(fptrs[ProfileFile[profile]], "%12.5g ", RvirValue[0][profile]);

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fprintf(fptrs[i], "\n#\n");

    /* Print column info. */

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL) {
	fprintf(fptrs[i], "#\n# COLUMN   DESCRIPTION\n");
	fprintf(fptrs[i], "#   %3d    %s\n", 1, "bin central radius (Mpc)");
	fprintf(fptrs[i], "#   %3d    %s\n", 2, "bin outer radius (Mpc)");
      }

    for (profile = 0; profile < MAX_PROFILES; profile++)
      if (ProfileName[profile] != NULL) {
	fprintf(fptrs[ProfileFile[profile]], "#   %3d    %s\n", 
		ColumnNumber[ProfileFile[profile]]++,
		ProfileName[profile]);
      }

    /* Print header. */
  
    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fprintf(fptrs[i], "#\n#%12.12s %12.12s ", "r_cent (Mpc)", 
		"r_edge (Mpc)");

    for (profile = 0; profile < MAX_PROFILES; profile++)
      if (ProfileName[profile] != NULL)
	fprintf(fptrs[ProfileFile[profile]], "%12.12s ", ProfileName[profile]);
    fprintf(fptrs[0], "%12.12s %12.12s %12.12s ", "vol_gas (real)", 
	    "vol (calc)", "vol_dm (real)");

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fprintf(fptrs[i], "\n");

    /* Loop over all radial bins, printing each. */

    for (j = 0; j < NumberOfPoints; j++) {

      float rmid = 0.5*(ProfileRadius[j]+ProfileRadius[j+1])*BoxSize;
      for (i = 0; i < NUMBER_OF_FILES; i++)
	if (fptrs[i] != NULL)
	  fprintf(fptrs[i], "%12.5g %12.5g ", 
		  rmid, ProfileRadius[j+1]*BoxSize);

      for (profile = 0; profile < MAX_PROFILES; profile++)
	if (ProfileName[profile] != NULL)
	  fprintf(fptrs[ProfileFile[profile]], "%12.5g ", 
		  (ProfileValue[j][profile] == 0) ? 
		  tiny_number : ProfileValue[j][profile]);

      fprintf(fptrs[0], "%12.5g %12.5g %12.5g ", ProfileWeight[j][0], 
	      (POW(ProfileRadius[j+1], 3) - POW(ProfileRadius[j], 3) ) * 
	      POW(BoxSize, 3) * 4.0*pi/3.0,
	      ProfileWeight[j][30]);

      for (i = 0; i < NUMBER_OF_FILES; i++)
	if (fptrs[i] != NULL)
	  fprintf(fptrs[i], "\n");

    }

    /* close files */

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fclose(fptrs[i]);

  } // end: loop over centers

  my_exit(EXIT_SUCCESS);
}

void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}

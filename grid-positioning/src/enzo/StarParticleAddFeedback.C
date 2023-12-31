/***********************************************************************
/
/  ADD FEEDBACK TO RADIAL PROFILE OVER MULTIPLE GRIDS
/
/  written by: John Wise
/  date:       September, 2005
/  modified1:
/
/ PURPOSE: To apply feedback effects, we must consider multiple grids
/          since sometimes the feedback radius often exceeds the grid
/          boundaries.  This routine makes sure that all of the grids
/          have the same code time to ensure consistency.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "StarParticleData.h"

#define MAX_TEMPERATURE 1e8

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RemoveParticles(LevelHierarchyEntry *LevelArray[], int level, int ID);
#ifdef USE_MPI
#endif /* USE_MPI */

int StarParticleAddFeedback(TopGridData *MetaData, 
			    LevelHierarchyEntry *LevelArray[], int level, 
			    Star *&AllStars, bool* &AddedFeedback)
{

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  Star *cstar;
  int i, l, dim, temp_int, SkipMassRemoval, SphereContained, count;
  float influenceRadius, RootCellWidth, SNe_dt;
  double EjectaThermalEnergy, EjectaDensity, EjectaMetalDensity;
  FLOAT Time;
  LevelHierarchyEntry *Temp;

  if (AllStars == NULL)
    return SUCCESS;

  /* Get time and SNe timestep */

  Temp = LevelArray[level];
  Time = Temp->GridData->ReturnTime();
  if (LastSupernovaTime < 0)
    SNe_dt = 0.0;
  else
    SNe_dt = Time - LastSupernovaTime;
  LastSupernovaTime = Time;
  RootCellWidth = 1.0 / MetaData->TopGridDims[0];

  /* Set the units. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  /* Initialize the AddedFeedback flag array */
  
  int TotalNumberOfStars = 0;
  for (cstar = AllStars; cstar; cstar = cstar->NextStar)
    TotalNumberOfStars++;
  if (TotalNumberOfStars > 0)
    AddedFeedback = new bool[TotalNumberOfStars];

  count = 0;
  for (cstar = AllStars; cstar; cstar = cstar->NextStar, count++) {

    AddedFeedback[count] = false;

    if(cstar->ReturnFeedbackFlag() != MBH_THERMAL) 
      if (!cstar->ApplyFeedbackTrue(SNe_dt))
	continue;

    /* Compute some parameters */

    cstar->CalculateFeedbackParameters(influenceRadius, RootCellWidth, 
           SNe_dt, EjectaDensity, EjectaThermalEnergy, EjectaMetalDensity, 
	   DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
	   VelocityUnits);

    /* Determine if a sphere with enough mass (or equivalently radius
       for SNe) is enclosed within grids on this level */

    if (cstar->FindFeedbackSphere(LevelArray, level, influenceRadius, 
	       EjectaDensity, SphereContained, SkipMassRemoval,	DensityUnits, 
	       LengthUnits, TemperatureUnits, TimeUnits, 
	       VelocityUnits) == FAIL) {
            ENZO_FAIL("Error in star::FindFeedbackSphere");
    }

    if (SphereContained == FALSE)
      continue;
    
    /* Now set cells within the radius to their values after feedback.
       While walking through the hierarchy, look for particle to
       change their properties to post-feedback values. */

    int CellsModified = 0;

    if (SkipMassRemoval == FALSE)
      for (l = level; l < MAX_DEPTH_OF_HIERARCHY; l++)
	for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->
	      AddFeedbackSphere(cstar, l, influenceRadius, VelocityUnits, 
				TemperatureUnits, EjectaDensity, 
				EjectaMetalDensity, EjectaThermalEnergy, 
				CellsModified) == FAIL) {
	    	    ENZO_FAIL("Error in AddFeedbackSphere.");
	  }

    /* Only kill a Pop III star after it has gone SN */

    if (cstar->ReturnFeedbackFlag() == SUPERNOVA)
      cstar->SetFeedbackFlag(DEATH);

    /* We only color the fields once */


    AddedFeedback[count] = true;

#ifdef UNUSED
    temp_int = CellsModified;
    MPI_Reduce(&temp_int, &CellsModified, 1, MPI_INT, MPI_SUM, ROOT_PROCESSOR,
	       MPI_COMM_WORLD);

    if (debug) {
      if (cstar->ReturnFeedbackFlag() != FORMATION)
	fprintf(stdout, "StarParticleAddFeedback[%"ISYM"][%"ISYM"]: "
		"Radius = %"GSYM" pc\n",
		cstar->ReturnID(), level, influenceRadius*LengthUnits/pc);
      if (cstar->ReturnFeedbackFlag() == SUPERNOVA || 
	  cstar->ReturnFeedbackFlag() == CONT_SUPERNOVA)
	fprintf(stdout, "StarParticleAddFeedback[%"ISYM"][%"ISYM"]: "
		"Energy = %"GSYM"  , skip = %"ISYM"\n",
		cstar->ReturnID(), level, EjectaThermalEnergy, SkipMassRemoval);
      fprintf(stdout, "StarParticleAddFeedback[%"ISYM"][%"ISYM"]: "
	      "changed %"ISYM" cells.\n", 
	      cstar->ReturnID(), level, CellsModified);
    }
#endif
    
  } // ENDFOR stars

  return SUCCESS;

}

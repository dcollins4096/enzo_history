/***********************************************************************
/
/  STAR PARTICLE INITIALIZATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE: Contains all routines to finalize the star particles.
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
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
#include "CommunicationUtilities.h"

#define NO_DEATH 0
#define KILL_STAR 1
#define KILL_ALL 2

int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids,
					 int TotalStarParticleCountPrevious[]);
int StarParticleAddFeedback(TopGridData *MetaData, 
			    LevelHierarchyEntry *LevelArray[], int level, 
			    Star* &AllStars, bool* &AddedFeedback);
int StarParticleAccretion(TopGridData *MetaData, 
			  LevelHierarchyEntry *LevelArray[], int level, 
			  Star *&AllStars);
int StarParticleSubtractAccretedMass(TopGridData *MetaData, 
				     LevelHierarchyEntry *LevelArray[], int level, 
				     Star *&AllStars);
int StarParticleDeath(LevelHierarchyEntry *LevelArray[], int level,
		      Star *&AllStars);
int CommunicationMergeStarParticle(HierarchyEntry *Grids[], int NumberOfGrids);
void DeleteStarList(Star * &Node);

int StarParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			 int level, Star *&AllStars,
			 int TotalStarParticleCountPrevious[])
{

  if (!StarParticleCreation && !StarParticleFeedback)
    return SUCCESS;

  int l, NumberOfStars;
  float TotalMass;
  Star *ThisStar, *MoveStar;
  LevelHierarchyEntry *Temp;
  FLOAT TimeNow;

  NumberOfStars = 0;
  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar)
    NumberOfStars++;
  bool *AddedFeedback = new bool[NumberOfStars];

  LCAPERF_START("StarParticleFinalize");

  /* Update the star particle counters. */

  CommunicationUpdateStarParticleCount(Grids, MetaData, NumberOfGrids,
				       TotalStarParticleCountPrevious);

  /* Update position and velocity of star particles from the actual
     particles */

  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar)
    ThisStar->UpdatePositionVelocity();


  /* Apply any stellar feedback onto the grids and add any gas to the
     accretion rates of the star particles */

  StarParticleAddFeedback(MetaData, LevelArray, level, 
			  AllStars, AddedFeedback);

  /* Update star particles for any accretion */

  StarParticleAccretion(MetaData, LevelArray, level, AllStars);

  /* Collect all sink particles and report the total mass to STDOUT */
  
  if (STARMAKE_METHOD(SINK_PARTICLE) && level == MaximumRefinementLevel) {
    TotalMass = 0.0;
    for (l = 0; l <= MaximumRefinementLevel; l++)
      for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel)
	TotalMass += Temp->GridData->ReturnTotalSinkMass();
#ifdef USE_MPI
    CommunicationReduceValues(&TotalMass, 1, MPI_SUM);
#endif
    if (debug)
      fprintf(stdout, "SinkParticle: Time = %"GOUTSYM", TotalMass = %"GSYM"\n", 
	      TimeNow, TotalMass);
  }

  /* Subtract gas from the grids that has accreted on to the star particles */

  StarParticleSubtractAccretedMass(MetaData, LevelArray, level, AllStars);  

  /* Check for any stellar deaths */

  StarParticleDeath(LevelArray, level, AllStars);

  /* 
     If the new particles are above a specified mass threshold,
     "activate" them.  Then check for any stellar deaths.

     Sync all star and normal particles that are stored in the grids
     to the global list (AllStars) so these changes are reflected
     there. 
  */

  int count = 0;
  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar, count++) {
    //TimeNow = LevelArray[ThisStar->ReturnLevel()]->GridData->ReturnTime();
    TimeNow = LevelArray[level]->GridData->ReturnTime();
//    if (debug) {
//      printf("AddedFeedback[%d] = %d\n", count, AddedFeedback[count]);
//      ThisStar->PrintInfo();
//    }
    if (AddedFeedback[count])
      ThisStar->ActivateNewStar(TimeNow);
    ThisStar->ResetAccretion();
    ThisStar->CopyToGrid();
    ThisStar->MirrorToParticle();

    // The pointers have been copied to the grid copy above, so we can
    // set the pointers in the global copy to NULL before deleting the stars.
    ThisStar->ResetAccretionPointers();

  } // ENDFOR stars

  /* Merge star particles */

  if (STARMAKE_METHOD(SINK_PARTICLE) && level == MaximumRefinementLevel) {  
    if (CommunicationMergeStarParticle(Grids, NumberOfGrids) == FAIL) {
      printf("CommunicationMergeStarParticle failed.\n");
      return FAIL;
    }
  }

  /* Delete the global star particle list, AllStars */

  DeleteStarList(AllStars);
  delete [] AddedFeedback;

  LCAPERF_STOP("StarParticleFinalize");
  return SUCCESS;

}

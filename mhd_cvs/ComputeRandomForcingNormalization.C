/***********************************************************************
/
/  COMPUTE RANDOM FORCING NORMALIZATION
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified1:
/
/  PURPOSE:
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
 
/* ======================================================================= */
/* Function prototypes. */
 
int   GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
			HierarchyEntry **Grids[]);
int   CommunicationAllSumValues(float *Values, int Number);
float CommunicationMinValue(float Value);
float CommunicationMaxValue(float Value);
 
/* This routine calculates the normalization for Random Forcing on
   level 0 (since this involves communication). */
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
				      int level, TopGridData *MetaData,
				      float * norm, float * bulkMomentum, 
				      float * pTopGridTimeStep)
{
 
  /* If level is above 0 then complain: forcing will only work on level 0
     grid(s). */
 
  if (level != 0) {
    fprintf(stderr, "Error in ComputeRandomForcingNormalization.\n");
    return FAIL;
  }
 
  /* Create an array (Grids) of all the grids on level 0. */
 
  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
 
  /* Loop over level 0 grids and compute sums over each of the grids. */
 
  int grid, GlobNum = 15;
  float * GlobVal = new float[GlobNum];
  for (int num = 0; num < GlobNum; num++)
    GlobVal[num] = 0.0;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->PrepareRandomForcingNormalization(GlobVal,
								 GlobNum)
	== FAIL) {
      fprintf(stderr, "Error in grid->PrepareRandomForcingNormalization.\n");
      return FAIL;
    }
 
  /* Communicate grid-specific sums and compute global sums;
     also communicate min/max Density values. */
 
  CommunicationAllSumValues(GlobVal, GlobNum-2);
  float minDens = CommunicationMinValue(*(GlobVal+GlobNum-2));
  float maxDens = CommunicationMaxValue(*(GlobVal+GlobNum-1));
 
  /* Compute normalization. */
 
  int numberOfGridZones = 1;
  for (int dim = 0; dim <  MetaData->TopGridRank; dim++)
    numberOfGridZones *= MetaData->TopGridDims[dim];
 
  float dt = LevelArray[level]->GridData->ReturnTimeStep();
  *pTopGridTimeStep = dt;
  if (RandomForcingEdot == 0.0)
    *norm = 0.0;
  else {
    float gv0 = GlobVal[0];
    float gv1 = GlobVal[1];
    if (gv0 < 1e-30 && gv0 > -1e-30){
      ERROR_MESSAGE
    }else{
      
      //From Alexei's fix, skipping the second order in norm^2 term
      //  *norm = 1.25*dt*RandomForcingEdot*numberOfGridZones/gv0;

      //From the original.
      *norm = ( sqrt(gv0*gv0 + gv1*dt*RandomForcingEdot*2.0*numberOfGridZones) - gv0 )/gv1;
    }
    for (int dim = 0; dim <  MetaData->TopGridRank; dim++)
      bulkMomentum[dim] = GlobVal[7+dim]/numberOfGridZones;
  }

  if (debug) printf("RandomForcingNormalization %.10g\n", *norm);
  if (debug) printf("RandomForcingGlobals: gv0 %g, gv1 %g, E_k %g M_m %g M_v %g px %g py %g pz %g norm %g \n",
		    GlobVal[0],
		    GlobVal[1],
		    0.5*GlobVal[4]/numberOfGridZones,
		    sqrt(GlobVal[2]/numberOfGridZones),
		    sqrt(GlobVal[3]/numberOfGridZones),
		    GlobVal[10]/numberOfGridZones,      // mean x-momentum
		    GlobVal[11]/numberOfGridZones,      // mean y-momentum
		    GlobVal[12]/numberOfGridZones,      // mean z-momentum
		    *norm);
 
  /* Output results for level 0 grids. */
 
  if (MetaData->CycleSkipGlobalDataDump != 0)
    if (MyProcessorNumber == ROOT_PROCESSOR &&
	MetaData->CycleNumber % MetaData->CycleSkipGlobalDataDump == 0.0) {
      FILE * Fptr;
      if ((Fptr = fopen("randomForcing.out", "a")) == NULL)
	ERROR_MESSAGE;
      fprintf( Fptr, "%d %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %10.5f\n",
	       MetaData->CycleNumber,
	       MetaData->Time,
	       0.50*GlobVal[4]/numberOfGridZones,   // kinetic energy
	       sqrt(GlobVal[2]/numberOfGridZones),  // mass weighted rms Mach
	       sqrt(GlobVal[3]/numberOfGridZones),  // volume weighed rms Mach
	       sqrt(GlobVal[5]/numberOfGridZones),  // rms Velocity
	       sqrt(GlobVal[6]/numberOfGridZones),  // Density variance
	       GlobVal[10]/numberOfGridZones,       // acquired x-momentum 
	       GlobVal[11]/numberOfGridZones,       // acquired y-momentum 
	       GlobVal[12]/numberOfGridZones,       // acquired z-momentum 
	       minDens, maxDens );                  // min/max Density
      fclose(Fptr);
    }
 
  /* clean up. */
 
  delete [] GlobVal;
 
  return SUCCESS;
}

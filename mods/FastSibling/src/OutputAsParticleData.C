/***********************************************************************
/
/  OUTPUTS GRID DATA AS PARTICLE DATA
/
/  written by: Greg Bryan
/  date:       August, 1996
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <df.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"

/* function prototypes */

int  DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid, 
				      TopGridData *MetaData, 
				      LevelHierarchyEntry *LevelArray[], 
				      int level);

int OutputAsParticleData(TopGridData &MetaData, 
			 LevelHierarchyEntry *LevelArray[],
			 int RegionStart[], int RegionEnd[], 
			 FLOAT RegionStartCoordinate[],
			 FLOAT RegionEndCoordinate[], int RegionLevel,
			 char *OutputFileName)
{

  /* Declarations */

  int i, j, k, dim, level, part, TotalRefineBy = 1;
  float TempCellWidth, BaseRadius;
  FLOAT RegionLeft[MAX_DIMENSION], RegionRight[MAX_DIMENSION];

  /* Set the Cell width of the root grid. */

  BaseRadius = (DomainRightEdge[0] - DomainLeftEdge[0])/
      float(MetaData.TopGridDims[0]);

  /* If undefined, set parameters. */

  if (RegionLevel == INT_UNDEFINED)
    RegionLevel = 0;

  for (level = 0; level < RegionLevel; level++)
    TotalRefineBy *= RefineBy;

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {

    /* If the start/end coordinate have been set, use them to set the
       indexes. */

    TempCellWidth = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
      float(MetaData.TopGridDims[dim]*TotalRefineBy);

    if (RegionStartCoordinate[dim] != FLOAT_UNDEFINED)
      RegionStart[dim] = nint((RegionStartCoordinate[dim] - 
				DomainLeftEdge[dim] ) / TempCellWidth );

    if (RegionEndCoordinate[dim] != FLOAT_UNDEFINED)
      RegionEnd[dim] = nint((RegionEndCoordinate[dim] - 
			      DomainLeftEdge[dim] ) / TempCellWidth ) - 1;

    /* If start/end indexes haven't been set, then set some default
       values. */

    if (RegionStart[dim] == INT_UNDEFINED)
      RegionStart[dim] = 0;
    if (RegionEnd[dim] == INT_UNDEFINED)
      RegionEnd[dim] = MetaData.TopGridDims[dim]*TotalRefineBy - 1;

    /* Find the position (this is the same as RegionStart/EndCoordinate
       if they are set). */

    RegionLeft[dim] = DomainLeftEdge[dim] + 
      (DomainRightEdge[dim] - DomainLeftEdge[dim])*
       FLOAT(RegionStart[dim])/FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);

    RegionRight[dim] = DomainLeftEdge[dim] + 
      (DomainRightEdge[dim] - DomainLeftEdge[dim])*
       FLOAT(RegionEnd[dim]+1)/FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);

  }

  if (debug)
    printf("OutputAsParticleData: Left = %g %g %g   Right = %g %g %g\n",
	   RegionLeft[0], RegionLeft[1], RegionLeft[2],
	   RegionRight[0], RegionRight[1], RegionRight[2]);

  /* Error check. */

  /* Initialize Particle List info */

  ListOfParticles *ListOfParticlesHead[NUM_PARTICLE_TYPES];
  int TotalNumberOfParticles[NUM_PARTICLE_TYPES];
  for (i = 0; i < NUM_PARTICLE_TYPES; i++) {
    ListOfParticlesHead[i] = NULL;
    TotalNumberOfParticles[i] = 0;
  }

  /* --------------------------------------------------------------- */
  /* Loop over all the levels */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    /* If SelfGravity, set all the particle mass fields. */

    LevelHierarchyEntry *Temp = LevelArray[level];
    if (SelfGravity)
      while (Temp != NULL) {
	DepositParticleMassField(Temp->GridHierarchyEntry);
	Temp = Temp->NextGridThisLevel;
      }

    /* Loop over all the grids. */

    Temp = LevelArray[level];
    while (Temp != NULL) {

      /* Allocate a new ListOfParticles for this grid. */

      for (i = 0; i < NUM_PARTICLE_TYPES; i++) {
	ListOfParticles *Temp = ListOfParticlesHead[i];
	ListOfParticlesHead[i] = new ListOfParticles;
	ListOfParticlesHead[i]->NextList = Temp;
      }

      /* Set particle density. */

      if (SelfGravity) {
	CopyOverlappingParticleMassFields(Temp->GridData, &MetaData, 
					  LevelArray, level);
	if (Temp->GridHierarchyEntry->ParentGrid != NULL)
	  Temp->GridHierarchyEntry->ParentGrid->GridData->DepositParticlePositions(Temp->GridData, Temp->GridHierarchyEntry->ParentGrid->GridData->ReturnTime(), GRAVITATING_MASS_FIELD_PARTICLES);
      }

      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */

      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);

      /* Zero the solution (on this grid) which is underneath any subgrid
	 (so we get only the high resolution solution from the subgrid). */

      LevelHierarchyEntry *Temp2 = LevelArray[level+1];
      if (level < RegionLevel)
	while (Temp2 != NULL) {
	  Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						   ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}

      /* Generate particle list for this grid. */

      Temp->GridData->OutputAsParticleData(RegionLeft, RegionRight,
					   ListOfParticlesHead, BaseRadius);
      for (i = 0; i < NUM_PARTICLE_TYPES; i++)
	TotalNumberOfParticles[i] += ListOfParticlesHead[i]->NumberOfParticles;

      /* Next grid on this level. */

      Temp = Temp->NextGridThisLevel;

    } // end loop over grids

  } // end loop over levels

  /* Write out the particles. */

  if (debug) {
    printf("TotalNumberOfParticles = %d %d %d\n", TotalNumberOfParticles[0],
	   TotalNumberOfParticles[1], TotalNumberOfParticles[2]);
    printf("Writing particle data to %s.\n", OutputFileName);
  }

  for (i = 0; i < NUM_PARTICLE_TYPES; i++) {

    char FileName[MAX_LINE_LENGTH];
    strcpy(FileName, OutputFileName);
    if (i == 0) strcat(FileName, ".gas");
    if (i == 1) strcat(FileName, ".dm");
    if (i == 2) strcat(FileName, ".star");

    if (TotalNumberOfParticles[i] > 0) {

      /* Allocate a field for these particles. */

      ListOfParticles FullList;
      FullList.NumberOfValues = ListOfParticlesHead[i]->NumberOfValues;
      if (debug) 
	printf("NumberOfValues[%d] = %d\n", i, FullList.NumberOfValues);
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	FullList.ParticlePosition[dim] = new float[TotalNumberOfParticles[i]];
	FullList.ParticleVelocity[dim] = new float[TotalNumberOfParticles[i]];
      }
      FullList.ParticleRadius = new float[TotalNumberOfParticles[i]];
      for (j = 0; j < FullList.NumberOfValues; j++)
	FullList.ParticleValue[j] = new float[TotalNumberOfParticles[i]];
      if (i == 1)
	FullList.ParticleIndex = new int[TotalNumberOfParticles[i]];

      /* Copy grid lists into the full list. */

      ListOfParticles *Temp = ListOfParticlesHead[i];
      int count = 0;
      while (Temp != NULL) {
//     printf("i=%d part=%d count=%d\n", i, Temp->NumberOfParticles, count);
	for (dim = 0; dim < MetaData.TopGridRank; dim++)
	  for (part = 0; part < Temp->NumberOfParticles; part++) {
	    FullList.ParticlePosition[dim][count+part] = 
	      Temp->ParticlePosition[dim][part];
	    FullList.ParticleVelocity[dim][count+part] = 
	      Temp->ParticleVelocity[dim][part];
	  }
	for (part = 0; part < Temp->NumberOfParticles; part++)
	  FullList.ParticleRadius[count+part] = Temp->ParticleRadius[part];
	for (j = 0; j < FullList.NumberOfValues; j++)
	  for (part = 0; part < Temp->NumberOfParticles; part++)
	    FullList.ParticleValue[j][count+part] = 
	      Temp->ParticleValue[j][part];
	if (i == 1)
	  for (part = 0; part < Temp->NumberOfParticles; part++)
	    FullList.ParticleIndex[count+part] = Temp->ParticleIndex[part];

	count += Temp->NumberOfParticles;
	Temp = Temp->NextList;
      }

      /* Set dimensions of HDF field */

      int32 TempInt = TotalNumberOfParticles[i];
      if (DFSDsetdims(1, &TempInt) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDsetdims.\n");
	return FAIL;
      }

      /* Create buffer. */

      float32 *buffer = new float32[TotalNumberOfParticles[i]];
      int ret;

      /* Write positions to HDF file. */

      for (j = 0; j < MetaData.TopGridRank; j++) {
	for (k = 0; k < TotalNumberOfParticles[i]; k++)
	  buffer[k] = FullList.ParticlePosition[j][k];
	if (j == 0) DFSDsetdatastrs("x_position", "", "e10.3", "");
	if (j == 1) DFSDsetdatastrs("y_position", "", "e10.3", "");
	if (j == 2) DFSDsetdatastrs("z_position", "", "e10.3", "");
	if (j == 0) 
	  ret = DFSDputdata(FileName, 1, &TempInt, (VOIDP) buffer);
	else
	  ret = DFSDadddata(FileName, 1, &TempInt, (VOIDP) buffer);
	if (ret == HDF_FAIL) {
	  fprintf(stderr, "Error writing %s pos (%d,%d).\n", FileName, i, j);
	  return FAIL;
	}
      }

      /* write velocities (note replication in baryons), sigh. */

      for (j = 0; j < MetaData.TopGridRank; j++) {
	for (k = 0; k < TotalNumberOfParticles[i]; k++)
	  buffer[k] = FullList.ParticleVelocity[j][k];
	if (j == 0) DFSDsetdatastrs("x_velocity", "km/s", "e10.3", "");
	if (j == 1) DFSDsetdatastrs("y_velocity", "km/s", "e10.3", "");
	if (j == 2) DFSDsetdatastrs("z_velocity", "km/s", "e10.3", "");
	if (DFSDadddata(FileName, 1, &TempInt, (VOIDP) buffer) == HDF_FAIL) {
	  fprintf(stderr, "Error writing %s vel (%d,%d).\n", FileName, i, j);
	  return FAIL;
	}
      }

      /* write radius. */

      DFSDsetdatastrs("radius", "box", "e10.3", "");
      for (k = 0; k < TotalNumberOfParticles[i]; k++)
	buffer[k] = FullList.ParticleRadius[k];
      if (DFSDadddata(FileName, 1, &TempInt, (VOIDP) buffer) == HDF_FAIL) {
	fprintf(stderr, "Error writing %s radius (%d).\n", FileName, i);
	return FAIL;
      }

      /* write other values. */

      for (j = 0; j < FullList.NumberOfValues; j++) {
	for (k = 0; k < TotalNumberOfParticles[i]; k++)
	  buffer[k] = FullList.ParticleValue[j][k];
	DFSDsetdatastrs("values", "", "e10.3", "");  // catch-all
	if (i == 0 && j == 0) DFSDsetdatastrs("baryon mass", "Msolar", "", "");
	if (i == 0 && j == 1) DFSDsetdatastrs("baryon density", "mean", "","");
	if (i == 0 && j == 2) DFSDsetdatastrs("temperature", "K", "", "");
	if (i >= 1 && j == 0) DFSDsetdatastrs("mass", "Msolar","","");
	if (i >= 1 && j == 1) DFSDsetdatastrs("density", "", "mean", "");
	if (DFSDadddata(FileName, 1, &TempInt, (VOIDP) buffer) == HDF_FAIL) {
	  fprintf(stderr, "Error writing %s value (%d,%d).\n", FileName, i, j);
	  return FAIL;
	}
      }

      /* write particle index. */

      if (i >= 1) {
	DFSDsetdatastrs("particle index", "", "", "");
	DFSDsetNT(DFNT_INT32);
	if (DFSDadddata(FileName, 1, &TempInt, (VOIDP) FullList.ParticleIndex) 
	    == HDF_FAIL) {
	  fprintf(stderr, "Error writing %s index (%d).\n", FileName, i);
	  return FAIL;
	}
	DFSDsetNT(DFNT_FLOAT32);
      }

      delete buffer;

    } // end: if (TotalNumberOfParticles[i] > 0)

  } // end: loop over dark matter/baryon particle lists

  return SUCCESS;
}

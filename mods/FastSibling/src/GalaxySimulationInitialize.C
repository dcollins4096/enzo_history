/***********************************************************************
/
/  INITIALIZE A GALAXY SIMULATION
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, March 2004
/
/  PURPOSE:
/
/    Set up a number of spherical objects
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <string.h>
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
#include "LevelHierarchy.h"
#include "TopGridData.h"

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int GalaxySimulationInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ColourName = "colour";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, disk, i;

  /* set default parameters */

  float GalaxySimulationGasMass[MAX_SPHERES],
        GalaxySimulationGalaxyMass[MAX_SPHERES],    
        GalaxySimulationDiskTemperature[MAX_SPHERES],
        GalaxySimulationAngularMomentum[MAX_SPHERES][MAX_DIMENSION],
        GalaxySimulationUniformVelocity[MAX_DIMENSION];
  float GalaxySimulationDiskRadius[MAX_SPHERES],
        GalaxySimulationDiskPosition[MAX_SPHERES][MAX_DIMENSION];
  float GalaxySimulationInitialTemperature,
        GalaxySimulationDiskScaleHeightz[MAX_SPHERES],
        GalaxySimulationDiskScaleHeightR[MAX_SPHERES],
        GalaxySimulationDarkMatterConcentrationParameter[MAX_SPHERES];
  int   GalaxySimulationNumberOfDisks,
        GalaxySimulationRefineAtStart,
        GalaxySimulationUseParticles,
        GalaxySimulationUseColour;
  
  //Default Values?

  GalaxySimulationNumberOfDisks   = 1;
  GalaxySimulationRefineAtStart   = TRUE;
  GalaxySimulationUseParticles    = FALSE;
  GalaxySimulationUseColour       = FALSE;
  GalaxySimulationInitialTemperature = 1000;
 
  for (disk = 0; disk < MAX_SPHERES; disk++) {
    GalaxySimulationDiskRadius[disk]     = 1.0;
    GalaxySimulationDiskTemperature[disk] = 1.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      GalaxySimulationDiskPosition[disk][dim] = 0.5*(DomainLeftEdge[dim] +
						     DomainRightEdge[dim]);
      GalaxySimulationAngularMomentum[disk][dim] = 0;
    }
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    GalaxySimulationUniformVelocity[dim] = 0;

  /* read input from file /exe/GalaxySimulation */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters from definition file Galaxy Simulation*/
    /* found at amr_mpi/exe/GalaxySimulation */

    ret += sscanf(line, "GalaxySimulationNumberOfDisks = %d",
		  &GalaxySimulationNumberOfDisks);
    ret += sscanf(line, "GalaxySimulationRefineAtStart = %d", 
		  &GalaxySimulationRefineAtStart);
    ret += sscanf(line, "GalaxySimulationUseParticles = %d", 
		  &GalaxySimulationUseParticles);
    ret += sscanf(line, "GalaxySimulationUseColour = %d", 
		  &GalaxySimulationUseColour);
    ret += sscanf(line, "GalaxySimulationInitialTemperature = %f", 
		  &GalaxySimulationInitialTemperature);
    ret += sscanf(line, "GalaxySimulationUniformVelocity = %f %f %f", 
                   &GalaxySimulationUniformVelocity[0], &GalaxySimulationUniformVelocity[1],
                   &GalaxySimulationUniformVelocity[2]);
    
    if (sscanf(line, "GalaxySimulationDiskRadius[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationDiskRadius[%d] = %"FSYM, &disk,
		    &GalaxySimulationDiskRadius[disk]);

    if (sscanf(line, "GalaxySimulationGalaxyMass[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationGalaxyMass[%d] = %"FSYM, &disk,
		    &GalaxySimulationGalaxyMass[disk]);

    if (sscanf(line, "GalaxySimulationGasMass[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationGasMass[%d] = %f", &disk,
		    &GalaxySimulationGasMass[disk]);

    if (sscanf(line, "GalaxySimulationDiskPosition[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationDiskPosition[%d] = %"FSYM" %"FSYM" %"FSYM, 
		    &disk, &GalaxySimulationDiskPosition[disk][0],
		    &GalaxySimulationDiskPosition[disk][1],
		    &GalaxySimulationDiskPosition[disk][2]);

    if (sscanf(line, "GalaxySimulationDiskScaleHeightz[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationDiskScaleHeightz[%d] = %"FSYM, &disk,
		    &GalaxySimulationDiskScaleHeightz[disk]);

    if (sscanf(line, "GalaxySimulationDiskScaleHeightR[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationDiskScaleHeightR[%d] = %"FSYM, &disk,
		    &GalaxySimulationDiskScaleHeightR[disk]);

    if (sscanf(line, "GalaxySimulationDarkMatterConcentrationParameter[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationDarkMatterConcentrationParameter[%d] = %"FSYM, &disk,
		    &GalaxySimulationDarkMatterConcentrationParameter[disk]);

    if (sscanf(line, "GalaxySimulationDiskTemperature[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationDiskTemperature[%d] = %f", &disk,
		    &GalaxySimulationDiskTemperature[disk]);

    if (sscanf(line, "GalaxySimulationAngularMomentum[%d]", &disk) > 0)
      ret += sscanf(line, "GalaxySimulationAngularMomentum[%d] = %f %f %f", 
		    &disk, &GalaxySimulationAngularMomentum[disk][0],
		    &GalaxySimulationAngularMomentum[disk][1],
		    &GalaxySimulationAngularMomentum[disk][2]);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "GalaxySimulation") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up grid */

  if (TopGrid.GridData->GalaxySimulationInitializeGrid(GalaxySimulationNumberOfDisks, 
						       GalaxySimulationDiskRadius,
						       GalaxySimulationGalaxyMass, 
						       GalaxySimulationGasMass,
						       GalaxySimulationDiskPosition, 
						       GalaxySimulationDiskScaleHeightz,
						       GalaxySimulationDiskScaleHeightR, 
						      GalaxySimulationDarkMatterConcentrationParameter,
						       GalaxySimulationDiskTemperature, 
						       GalaxySimulationInitialTemperature,
						       GalaxySimulationUseParticles,
						       GalaxySimulationUseColour, 
						       GalaxySimulationAngularMomentum,
						       GalaxySimulationUniformVelocity, 0) == FAIL){
    fprintf(stderr, "Error in GalaxySimulationInitializeGrid.\n");
    return FAIL;
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */

  if (GalaxySimulationRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->GalaxySimulationInitializeGrid(
	     GalaxySimulationNumberOfDisks, GalaxySimulationDiskRadius,
	     GalaxySimulationGalaxyMass, GalaxySimulationGasMass,
	     GalaxySimulationDiskPosition, GalaxySimulationDiskScaleHeightz,
	     GalaxySimulationDiskScaleHeightR, GalaxySimulationDarkMatterConcentrationParameter,
	     GalaxySimulationDiskTemperature, GalaxySimulationInitialTemperature,
	     GalaxySimulationUseParticles,
	     GalaxySimulationUseColour, GalaxySimulationAngularMomentum,
             GalaxySimulationUniformVelocity, level+1) == FAIL) {
	  fprintf(stderr, "Error in GalaxySimulationInitializeGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
				   *LevelArray[level-1]->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (GalaxySimulationRefineAtStart)

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if (GalaxySimulationUseColour)
    DataLabel[count++] = ColourName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "GalaxySimulationNumberOfDisks      = %d\n",
	    GalaxySimulationNumberOfDisks);
    fprintf(Outfptr, "GalaxySimulationRefineAtStart      = %d\n",
	    GalaxySimulationRefineAtStart);
    fprintf(Outfptr, "GalaxySimulationUseParticles       = %d\n",
	    GalaxySimulationUseParticles);
    fprintf(Outfptr, "GalaxySimulationUseColour          = %d\n",
	    GalaxySimulationUseColour);
    fprintf(Outfptr, "GalaxySimulationInitialTemperature = %f\n",
	    GalaxySimulationInitialTemperature);
    fprintf(Outfptr, "GalaxySimulationUniformVelocity    = %f %f %f\n",
	    GalaxySimulationUniformVelocity[0], GalaxySimulationUniformVelocity[1],
	    GalaxySimulationUniformVelocity[2]);

    for (disk = 0; disk < GalaxySimulationNumberOfDisks; disk++) {
      fprintf(Outfptr, "GalaxySimulationDiskRadius[%d] = %"GOUTSYM"\n", disk,
	      GalaxySimulationDiskRadius[disk]);

      fprintf(Outfptr, "GalaxySimulationGalaxyMass[%d] = %"GOUTSYM"\n", disk,
	      GalaxySimulationGalaxyMass[disk]);

      fprintf(Outfptr, "GalaxySimulationGasMass[%d] = %f\n", disk,
	      GalaxySimulationGasMass[disk]);

      fprintf(Outfptr, "GalaxySimulationDiskScaleHeightz[%d] = %f\n", disk,
	      GalaxySimulationDiskScaleHeightz[disk]);
      
      fprintf(Outfptr, "GalaxySimulationDiskScaleHeightR[%d] = %f\n", disk,
	      GalaxySimulationDiskScaleHeightR[disk]);

      fprintf(Outfptr, "GalaxySimulationDarkMatterConcentrationParameter[%d] = %f\n", disk,
	      GalaxySimulationDarkMatterConcentrationParameter[disk]);

      fprintf(Outfptr, "GalaxySimulationDiskTemperature[%d] = %f\n", disk,
	      GalaxySimulationDiskTemperature[disk]);
  
      fprintf(Outfptr, "GalaxySimulationDiskPosition[%d] = ", disk);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			GalaxySimulationDiskPosition[disk]);

      fprintf(Outfptr, "GalaxySimulationAngularMomentum[%d] = ", disk);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			GalaxySimulationAngularMomentum[disk]);
    }
  }

  return SUCCESS;

}

/***********************************************************************
/
/  INITIALIZE A MHD COLLAPSE TEST
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

int MHDCollapseTestInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  if( MHD_Used != TRUE ){
    fprintf(stderr, "Don't use MHDCollapseTest if MHD_Used is NOT on.\n");
    return FAIL;
  }
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ColourName = "colour";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int CollapseTestNumberOfSpheres = 1;
  int CollapseTestRefineAtStart   = TRUE;
  int CollapseTestUseParticles    = FALSE;
  int CollapseTestUseColour       = FALSE;
  float CollapseTestInitialTemperature = 1000;
  int   CollapseTestSphereType[MAX_SPHERES];
  float CollapseTestSphereBackgroundDensity = 1.0;
  float CollapseTestSphereDensity[MAX_SPHERES],
        CollapseTestSphereTemperature[MAX_SPHERES],
        CollapseTestSphereVelocity[MAX_SPHERES][MAX_DIMENSION],
        CollapseTestUniformVelocity[MAX_DIMENSION];
  FLOAT CollapseTestSphereRadius[MAX_SPHERES],
        CollapseTestSphereCoreRadius[MAX_SPHERES],
        CollapseTestSpherePosition[MAX_SPHERES][MAX_DIMENSION];
  float MHDCollapseTestMagneticFieldx = 0.0;
  float MHDCollapseTestMagneticFieldy = 0.0;
  float MHDCollapseTestMagneticFieldz = 0.0;

  for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    CollapseTestSphereRadius[sphere]     = 1.0;
    CollapseTestSphereCoreRadius[sphere] = 0.1;
    CollapseTestSphereDensity[sphere]    = 1.0;
    CollapseTestSphereTemperature[sphere] = 1.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CollapseTestSpherePosition[sphere][dim] = 0.5*(DomainLeftEdge[dim] +
						     DomainRightEdge[dim]);
      CollapseTestSphereVelocity[sphere][dim] = 0;
    }
    CollapseTestSphereType[sphere]       = 0;
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    CollapseTestUniformVelocity[dim] = 0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "CollapseTestNumberOfSpheres = %d",
		  &CollapseTestNumberOfSpheres);
    ret += sscanf(line, "CollapseTestRefineAtStart = %d", 
		  &CollapseTestRefineAtStart);
    ret += sscanf(line, "CollapseTestUseParticles = %d", 
		  &CollapseTestUseParticles);
    ret += sscanf(line, "CollapseTestUseColour = %d", 
		  &CollapseTestUseColour);
    ret += sscanf(line, "CollapseTestInitialTemperature = %"FSYM, 
		  &CollapseTestInitialTemperature);
    ret += sscanf(line, "CollapseTestUniformVelocity = %"FSYM" %"FSYM" %"FSYM,
		  CollapseTestUniformVelocity, CollapseTestUniformVelocity+1,
		  CollapseTestUniformVelocity+2);
    if (sscanf(line, "CollapseTestSphereType[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereType[%d] = %d", &sphere,
		    &CollapseTestSphereType[sphere]);
    if (sscanf(line, "CollapseTestSphereRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereRadius[%d] = %"PSYM, &sphere,
		    &CollapseTestSphereRadius[sphere]);
    if (sscanf(line, "CollapseTestSphereCoreRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereCoreRadius[%d] = %"PSYM, &sphere,
		    &CollapseTestSphereCoreRadius[sphere]);
					
    ret += sscanf(line, "CollapseTestSphereBackgroundDensity = %"FSYM, 
		  &CollapseTestSphereBackgroundDensity);
    if (sscanf(line, "CollapseTestSphereDensity[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereDensity[%d] = %"FSYM, &sphere,
		    &CollapseTestSphereDensity[sphere]);
    if (sscanf(line, "CollapseTestSphereTemperature[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereTemperature[%d] = %"FSYM, &sphere,
		    &CollapseTestSphereTemperature[sphere]);
    if (sscanf(line, "CollapseTestSpherePosition[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSpherePosition[%d] = %"PSYM" %"PSYM" %"PSYM, 
		    &sphere, &CollapseTestSpherePosition[sphere][0],
		    &CollapseTestSpherePosition[sphere][1],
		    &CollapseTestSpherePosition[sphere][2]);
    if (sscanf(line, "CollapseTestSphereVelocity[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereVelocity[%d] = %"FSYM" %"FSYM" %"FSYM,
		    &sphere, &CollapseTestSphereVelocity[sphere][0],
		    &CollapseTestSphereVelocity[sphere][1],
		    &CollapseTestSphereVelocity[sphere][2]);

     ret += sscanf(line, "MHDCollapseTestMagneticFieldx = %"FSYM,
                  &MHDCollapseTestMagneticFieldx);
     ret += sscanf(line, "MHDCollapseTestMagneticFieldy = %"FSYM,
                  &MHDCollapseTestMagneticFieldy);
     ret += sscanf(line, "MHDCollapseTestMagneticFieldz = %"FSYM,
                  &MHDCollapseTestMagneticFieldz);



    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "CollapseTest") 
	&& line[0] != '#')
      fprintf(stderr, "MHD Collapse warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up grid */

  if (TopGrid.GridData->MHDCollapseTestInitializeGrid(
	     CollapseTestNumberOfSpheres, CollapseTestSphereRadius,
	     CollapseTestSphereCoreRadius, 
	     CollapseTestSphereBackgroundDensity,
	     CollapseTestSphereDensity,
	     CollapseTestSphereTemperature,
	     CollapseTestSpherePosition, CollapseTestSphereVelocity,
	     CollapseTestSphereType, CollapseTestUseParticles,
             CollapseTestUniformVelocity, CollapseTestUseColour,
             CollapseTestInitialTemperature, 0,
             MHDCollapseTestMagneticFieldx,
             MHDCollapseTestMagneticFieldy,
             MHDCollapseTestMagneticFieldz) == FAIL) {
    fprintf(stderr, "Error in CollapseTestInitializeGrid.\n");
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

  if (CollapseTestRefineAtStart) {

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
	if (Temp->GridData->MHDCollapseTestInitializeGrid(
	     CollapseTestNumberOfSpheres, CollapseTestSphereRadius,
	     CollapseTestSphereCoreRadius, 
	     CollapseTestSphereBackgroundDensity,
	     CollapseTestSphereDensity,
	     CollapseTestSphereTemperature,
	     CollapseTestSpherePosition, CollapseTestSphereVelocity,
	     CollapseTestSphereType, CollapseTestUseParticles,
	     CollapseTestUniformVelocity, CollapseTestUseColour,
	     CollapseTestInitialTemperature, level+1,
             MHDCollapseTestMagneticFieldx,
             MHDCollapseTestMagneticFieldy,
             MHDCollapseTestMagneticFieldz) == FAIL) {
	  fprintf(stderr, "Error in MHDCollapseTestInitializeGrid.\n");
	  return FAIL;
	}
       	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency among levels. */
  
    //Projection should happen on the Magnetic field for initialization.
    int MHD_ProjectEtmp = MHD_ProjectE;
    int MHD_ProjectBtmp = MHD_ProjectB;
    MHD_ProjectE=FALSE;
    MHD_ProjectB=TRUE;


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

  MHD_ProjectE = MHD_ProjectEtmp;
  MHD_ProjectB = MHD_ProjectBtmp;

 
  } // end: if (CollapseTestRefineAtStart)

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
#ifdef ATHENA
  if( EquationOfState == 0 )
#endif //ATHENA
    DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if (CollapseTestUseColour)
    DataLabel[count++] = ColourName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  MHDcLabel[0] = "MagneticField_C_1";
  MHDcLabel[1] = "MagneticField_C_2";
  MHDcLabel[2] = "MagneticField_C_3";
                                                                                                                                                             
  MHDLabel[0] = "MagneticField_F_1";
  MHDLabel[1] = "MagneticField_F_2";
  MHDLabel[2] = "MagneticField_F_3";
                                                                                                                                                             
  MHDeLabel[0] = "ElectricField_1";
  MHDeLabel[1] = "ElectricField_2";
  MHDeLabel[2] = "ElectricField_3";
                                                                                                                                                             
  CurrentLabel[0] = "Current_1";
  CurrentLabel[1] = "Current_2";
  CurrentLabel[2] = "Current_3";


  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CollapseTestNumberOfSpheres    = %d\n",
	    CollapseTestNumberOfSpheres);
    fprintf(Outfptr, "CollapseTestRefineAtStart      = %d\n",
	    CollapseTestRefineAtStart);
    fprintf(Outfptr, "CollapseTestUseParticles       = %d\n",
	    CollapseTestUseParticles);
    fprintf(Outfptr, "CollapseTestUseColour          = %d\n",
	    CollapseTestUseColour);
    fprintf(Outfptr, "CollapseTestInitialTemperature = %f\n",
	    CollapseTestInitialTemperature);
    fprintf(Outfptr, "CollapseTestUniformVelocity    = %f %f %f\n",
	    CollapseTestUniformVelocity[0], CollapseTestUniformVelocity[1],
	    CollapseTestUniformVelocity[2]);
    for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++) {
      fprintf(Outfptr, "CollapseTestSphereType[%d] = %d\n", sphere,
	      CollapseTestSphereType[sphere]);
      fprintf(Outfptr, "CollapseTestSphereRadius[%d] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereRadius[sphere]);
      fprintf(Outfptr, "CollapseTestSphereCoreRadius[%d] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereCoreRadius[sphere]);
      fprintf(Outfptr, "CollapseTestSphereDensity[%d] = %f\n", sphere,
	      CollapseTestSphereDensity[sphere]);
      fprintf(Outfptr, "CollapseTestSphereTemperature[%d] = %f\n", sphere,
	      CollapseTestSphereTemperature[sphere]);
      fprintf(Outfptr, "CollapseTestSpherePosition[%d] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CollapseTestSpherePosition[sphere]);
      fprintf(Outfptr, "CollapseTestSphereVelocity[%d] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CollapseTestSphereVelocity[sphere]);
    }
  }

  return SUCCESS;

}

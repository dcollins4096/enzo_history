/***********************************************************************
/
/  INITIALIZE MHD 3D TEST
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

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

int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData, 
		 ExternalBoundary *Exterior, FLOAT WriteTime);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int MHD3DTestInitialize(FILE *fptr, FILE *Outfptr, 
			HierarchyEntry &TopGrid, TopGridData &MetaData) 
{
  char *DensName = "Density";
  char *PresName = "Pressure";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ColourName = "colour";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = FALSE;
  int MHD3DProblemType = 0;
  float  rhol = 1.0, rhou = 1.0, 
    vxl = 0, vxu = 0, 
    vyl = 0, vyu = 0, 
    pl = 1.0, pu = 1.0,
    Bxl = 0.0, Bxu = 0.0,
    Byl = 0.0, Byu = 0.0;
  
  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "RefineAtStart = %d", 
		  &RefineAtStart);
    ret += sscanf(line, "LeftVelocityX = %f",
		  &vxl);
    ret += sscanf(line, "LeftVelocityY = %f",
		  &vyl);
    ret += sscanf(line, "LeftPressure = %f", 
		  &pl);
    ret += sscanf(line, "LeftDensity = %f", 
		  &rhol);
    ret += sscanf(line, "LeftBx = %f",
		  &Bxl);
    ret += sscanf(line, "LeftBy = %f",
		  &Byl);
    ret += sscanf(line, "RightVelocityX = %f", 
		  &vxu);
    ret += sscanf(line, "RightVelocityY = %f", 
		  &vyu);
    ret += sscanf(line, "RightPressure = %f", 
		  &pu);
    ret += sscanf(line, "RightDensity = %f",
                  &rhou);
    ret += sscanf(line, "RightBx = %f",
		  &Bxu);
    ret += sscanf(line, "RightBy = %f",
		  &Byu);
    ret += sscanf(line, "MHD3DProblemType = %d",
		  &MHD3DProblemType);

  } // end input from parameter file

  float DensityUnits = 1, LengthUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  /* set up grid */

  if (TopGrid.GridData->MHD3DTestInitializeGrid(MHD3DProblemType,
						rhol, rhou,
						vxl,  vxu,
						vyl,  vyu,
						pl,   pu,
						Bxl,  Bxu,
						Byl,  Byu)  == FAIL) {
    fprintf(stderr, "Error in MHD3DTestInitializeGrid.\n");
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

  if (RefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      printf("In level %i\n", level);
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->MHD3DTestInitializeGrid(MHD3DProblemType,
						    rhol, rhou,
						    vxl,  vxu,
						    vyl,  vyu,
						    pl,   pu,
						    Bxl,  Bxu,
						    Byl,  Byu) == FAIL) {
	  fprintf(stderr, "Error in MHD3DTestInitializeGrid.\n");
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

    //WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
    //       &TopGrid, MetaData, Exterior, -1);

  } // end: if (RefineAtStart)


  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = GEName;
  }
  DataLabel[count++] = BxName;
  DataLabel[count++] = ByName;
  DataLabel[count++] = BzName;
  DataLabel[count++] = PhiName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  /*if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "RefineAtStart      = %d\n",
	    RefineAtStart);
    fprintf(Outfptr, "LeftDensity       = %f\n",
	    rhol);
    fprintf(Outfptr, "RightDensity          = %f\n",
	    rhor);
    fprintf(Outfptr, "LeftVelocityX = %f\n",
	    vxl);
    fprintf(Outfptr, "RightVelocityX = %f\n",
            vxr);
    fprintf(Outfptr, "LeftVelocityY = %f\n",
	    vyl);
    fprintf(Outfptr, "RightVelocityY = %f\n",
            vyr);
    fprintf(Outfptr, "LeftBx = %f\n",
	    Bxl);
    fprintf(Outfptr, "RightBx = %f\n",
	    Bxr);
    fprintf(Outfptr, "LeftBy = %f\n",
	    Byl);
    fprintf(Outfptr, "RightBy = %f\n",
	    Byr);
    fprintf(Outfptr, "LeftPressure = %f\n",
            pl);
    fprintf(Outfptr, "RightPressure = %f\n",
            pr);
	    }*/

  return SUCCESS;

}


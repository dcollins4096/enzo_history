/***********************************************************************
/
/  PROJECTS A SECTION OF THE GRID TO A PLANE
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// This function projects a specified region of the grid to a plane,
//   along one of the orthogonal directions.

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <df.h>
#ifdef USE_MPI
#include "mpi.h"
#endif
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

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
int  DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid, 
				      TopGridData *MetaData, 
				      LevelHierarchyEntry *LevelArray[], 
				      int level);
#ifdef USE_MPI
int CommunicationReduceValues(float *values, int number, MPI_Op ReduceOp);
#endif /* USE_MPI */

#define NUMBER_OF_PROJECTED_FIELDS 12


int ProjectToPlane(TopGridData &MetaData, LevelHierarchyEntry *LevelArray[],
		   int ProjectStartTemp[], int ProjectEndTemp[], 
		   FLOAT ProjectStartCoordinate[],
		   FLOAT ProjectEndCoordinate[], int ProjectLevel,
		   int ProjectionDimension, char *ProjectionFileName,
		   int ProjectionSmooth, ExternalBoundary *Exterior)
{

  /* Declarations */

  int i, j, dim, field, level, ret, size = 1;
  int ProjectDim[MAX_DIMENSION];
  float *ProjectedField[NUMBER_OF_PROJECTED_FIELDS], TempCellWidth;
  FLOAT ProjectLeft[MAX_DIMENSION], ProjectRight[MAX_DIMENSION];

  /* Set The GravityResolution to 1 to make the DM resolution the same at
     the gas. This is one of two reason's why this routine cannot be called
     during the evolution. */

  GravityResolution = 1;

  /* Copy ProjectStart and ProjectEnd into long int version so we can do
     deep projections. */

  long_int ProjectStart[MAX_DIMENSION], ProjectEnd[MAX_DIMENSION],
           TotalRefineBy = 1;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    ProjectStart[dim] = ProjectStartTemp[dim];
    ProjectEnd[dim]   = ProjectEndTemp[dim];
  }

  /* If undefined, set parameters. */

  if (ProjectLevel == INT_UNDEFINED)
    ProjectLevel = 0;

  for (level = 0; level < ProjectLevel; level++)
    TotalRefineBy *= RefineBy;

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {

    /* If the start/end coordinate have been set, use them to set the
       indexes. */

    TempCellWidth = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
      float(MetaData.TopGridDims[dim]*TotalRefineBy);

    if (ProjectStartCoordinate[dim] != FLOAT_UNDEFINED)
      ProjectStart[dim] = nlongint((ProjectStartCoordinate[dim] - 
				DomainLeftEdge[dim] ) / TempCellWidth );

    if (ProjectEndCoordinate[dim] != FLOAT_UNDEFINED)
      ProjectEnd[dim] = nlongint((ProjectEndCoordinate[dim] - 
			      DomainLeftEdge[dim] ) / TempCellWidth ) - 1;

    /* If start/end indexes haven't been set, then set some default
       values. */

    if (ProjectStart[dim] == INT_UNDEFINED)
      ProjectStart[dim] = 0;
    if (ProjectEnd[dim] == INT_UNDEFINED)
      ProjectEnd[dim] = MetaData.TopGridDims[dim]*TotalRefineBy - 1;

    /* Compute the dimension and the size. */

    ProjectDim[dim] = ProjectEnd[dim] - ProjectStart[dim] + 1;

    if (dim != ProjectionDimension)
      size *= ProjectDim[dim];

    /* Find the position (this is the same as ProjectStart/EndCoordinate
       if they are set). */

    ProjectLeft[dim] = DomainLeftEdge[dim] + 
      (DomainRightEdge[dim] - DomainLeftEdge[dim])*
       FLOAT(ProjectStart[dim])/FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);

    ProjectRight[dim] = DomainLeftEdge[dim] + 
      (DomainRightEdge[dim] - DomainLeftEdge[dim])*
       FLOAT(ProjectEnd[dim]+1)/FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);

  }

  if (debug)
    printf("ProjectToPlane: Left = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"   Right = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
	   ProjectLeft[0], ProjectLeft[1], ProjectLeft[2],
	   ProjectRight[0], ProjectRight[1], ProjectRight[2]);

  /* Error check. */

  if (ProjectionDimension < 0 || ProjectionDimension > MetaData.TopGridRank) {
    fprintf(stderr, "Invalid ProjectionDimension (%d).\n",ProjectionDimension);
    return FAIL;
  }

  /* Check to see if the file ProjectParameters exists.  If it does, read
     the parameters within it. */

  FILE *fptr;
  char line[MAX_LINE_LENGTH], 
    XrayTableFileName[MAX_LINE_LENGTH] = "lookup_metal0.3.data";
  float XrayLowerCutoffkeV = 0.5, XrayUpperCutoffkeV = 2.5;
  int XrayUseLookupTable = FALSE;

  if ((fptr = fopen("ProjectionParameters", "r")) != NULL) {
    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
      sscanf(line, "XrayLowerCutoffkeV = %f", &XrayLowerCutoffkeV);
      sscanf(line, "XrayUpperCutoffkeV = %f", &XrayUpperCutoffkeV);
      sscanf(line, "XrayTableFileName = %s", XrayTableFileName);
    }
    fclose(fptr);
    printf("XrayLowerCutoffkeV = %g\n", XrayLowerCutoffkeV);
    printf("XrayUpperCutoffkeV = %g\n", XrayUpperCutoffkeV);
    printf("XrayTableFileName = %s\n", XrayTableFileName);
    XrayUseLookupTable = TRUE;
  }

  /* Allocate plane. */

  for (field = 0; field < NUMBER_OF_PROJECTED_FIELDS; field++) {

    ProjectedField[field] = new float[size];

    for (i = 0; i < size; i++)
      ProjectedField[field][i] = 0;
  }

  /* --------------------------------------------------------------- */
  /* Loop over the levels down to the requested one. */

  for (level = 0; level <= ProjectLevel; level++) {

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

      /* Set particle density from parent (usually small effect). */

#ifdef UNUSED
      if (SelfGravity) {
	CopyOverlappingParticleMassFields(Temp->GridData, &MetaData, 
					  LevelArray, level);
	//      if (Temp->GridHierarchyEntry->ParentGrid != NULL)
	//	  Temp->GridHierarchyEntry->ParentGrid->GridData->DepositParticlePositions(Temp->GridData, Temp->GridHierarchyEntry->ParentGrid->GridData->ReturnTime(), GRAVITATING_MASS_FIELD_PARTICLES);
      }
#endif

      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */

      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);

      /* Set boundary (interpolate from parent is good enough) . */

      if (level > 0) {
        /* the following line make sure the times are identical so no
           time-interpolation is required. */
	Temp->GridHierarchyEntry->ParentGrid->GridData->SetTime(
				             Temp->GridData->ReturnTime());
	Temp->GridData->InterpolateBoundaryFromParent
	                   (Temp->GridHierarchyEntry->ParentGrid->GridData);
	CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, level);
      } else
	Temp->GridData->SetExternalBoundaryValues(Exterior);

      /* Zero the solution (on this grid) which is underneath any subgrid
	 (so we get only the high resolution solution from the subgrid). */

      LevelHierarchyEntry *Temp2 = LevelArray[level+1];
      if (level < ProjectLevel)
	while (Temp2 != NULL) {
	  Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						   ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}

      /* Project to plane. */

      Temp->GridData->ProjectToPlane(ProjectLeft, ProjectRight,
				     ProjectDim, ProjectedField,
				     ProjectionDimension,
				     ProjectionSmooth,
				     NUMBER_OF_PROJECTED_FIELDS, level,
				     XrayUseLookupTable, XrayLowerCutoffkeV,
				     XrayUpperCutoffkeV, XrayTableFileName);

      /* Next grid on this level. */

      Temp = Temp->NextGridThisLevel;

    } // end loop over grids

  } // end loop over levels

  /* Loop over the rest of levels if doing star particle (to complete
     star particle density fields -- the -1 means only star particle stuff). */

  for (level = ProjectLevel+1; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->ProjectToPlane(ProjectLeft, ProjectRight,
				     ProjectDim, ProjectedField,
				     ProjectionDimension,
				     ProjectionSmooth,
				     NUMBER_OF_PROJECTED_FIELDS, -1,
				     XrayUseLookupTable, XrayLowerCutoffkeV,
				     XrayUpperCutoffkeV, XrayTableFileName);
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* Sum the contribution from each processor (except for the level field,
     for which we take the max value). */
     
#ifdef USE_MPI
  for (field = 0; field < NUMBER_OF_PROJECTED_FIELDS; field++)
    if (field == 4)
      CommunicationReduceValues(ProjectedField[field], size, MPI_MAX);
    else
      CommunicationReduceValues(ProjectedField[field], size, MPI_SUM);
#endif /* USE_MPI */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

  /* Write out the projected planes. */

  if (debug) printf("Writing projected planes to %s.\n", ProjectionFileName);

  /* Set dimensions (reversed since this is c and we're using f77 order). */

  int32 OutDims[2];
  for (dim = 0, i = 0; dim < MetaData.TopGridRank; dim++)
    if (dim != ProjectionDimension)
      OutDims[1 - i++] = ProjectDim[dim];

  if (DFSDsetdims(2, OutDims) == HDF_FAIL) {
    fprintf(stderr, "Error in DFSDsetdims.\n");
    return FAIL;
  }

  /* Compute the scales. */

  float32 *float_temp;
  for (dim = 0, i = 0; dim < MetaData.TopGridRank; dim++)
    if (dim != ProjectionDimension) {

      float_temp = new float32[ProjectDim[dim]];

      for (j = 0; j < ProjectDim[dim]; j++)
	float_temp[j] = float32(ProjectLeft[dim] + 
	  (float(j)+0.5) * (ProjectRight[dim]-ProjectLeft[dim]) / 
	    float(ProjectDim[dim]));

      if (DFSDsetdimscale(2 - i++, ProjectDim[dim], (VOIDP) float_temp) == 
	  HDF_FAIL) {
	fprintf(stderr, "Error in DFSDsetdimscale.\n");
	return FAIL;
      }

      delete [] float_temp;
    }

  /* Write projected planes to HDF file. */

  float_temp = new float32[size];
  for (i = 0; i < NUMBER_OF_PROJECTED_FIELDS; i++) {

    /* Set name of field. */

    if (i == 0) 
      DFSDsetdatastrs("projected_gas_density", "Msolar/Mpc^2", "e10.3", "");
    if (i == 1) 
      DFSDsetdatastrs("projected_x-ray_luminosity_div1e23", "", "e10.3", "");
    if (i == 2) 
      DFSDsetdatastrs("projected_dm_density", "Msolar/Mpc^2", "e10.3", "");
    if (i == 3) 
      DFSDsetdatastrs("projected_x-ray_weighted_temperature", "K", "e10.3","");
    if (i == 4)
      DFSDsetdatastrs("projected_level", "", "f8.3","");
    if (i == 5)
      DFSDsetdatastrs("SZ y effect", "", "f8.3","");
    if (i == 6)
      DFSDsetdatastrs("DT/T Doppler effect", "", "f8.3","");
    if (i == 7)
      DFSDsetdatastrs("Metallicity", "", "f8.3","");
    if (i == 8)
      DFSDsetdatastrs("projected_star_density", "", "f8.3","");
    if (i == 9)
      DFSDsetdatastrs("density of low entropy gas", "Msolar/Mpc^2", "f8.3","");
    if (i == 10)
      DFSDsetdatastrs("mean entropy of gas", "keV cm^-2", "f8.3","");
    if (i == 11)
      DFSDsetdatastrs("inverse cooling time", "1/code units", "f8.3","");
    if (i == 12)
      DFSDsetdatastrs("OVII column density", "cm^-2", "f8.3","");

    /* for i == 3, divide by x-ray luminosity. */

    if (i == 3)
      for (j = 0; j < size; j++)
        ProjectedField[3][j] /= ProjectedField[1][j];

    /* for i == 7/10/11, divide by baryon density. */

    if (i == 7 || i == 10 || i == 11)
      for (j = 0; j < size; j++)
        ProjectedField[7][j] /= ProjectedField[0][j];

    /* Copy to float32 space. */

    for (j = 0; j < size; j++)
      float_temp[j] = float32(ProjectedField[i][j]);

    /* Write data. */

    if (i == 0)
      ret = DFSDputdata(ProjectionFileName, 2, OutDims, (VOIDP) float_temp);
    if (i != 0)
      ret = DFSDadddata(ProjectionFileName, 2, OutDims, (VOIDP) float_temp);

    if (ret == HDF_FAIL) {
      fprintf(stderr, "Error in DFSDput/adddata.\n");
      return FAIL;
    }

  }
  delete [] float_temp;

  } // end: if (MyProcessNumber == ROOT_PROCESSOR)

  return SUCCESS;
}

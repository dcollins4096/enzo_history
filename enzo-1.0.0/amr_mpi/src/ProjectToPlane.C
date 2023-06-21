/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  PROJECTS A SECTION OF THE GRID TO A PLANE
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:  Robert Harkness, July 2002
/  modified2:  James Bordner, June 2003   Merged Harkness HDF4_* HDF5_*
/
/  PURPOSE:
/
************************************************************************/

// This function projects a specified region of the grid to a plane,
//   along one of the orthogonal directions.

#ifdef USE_HDF4
#include <df.h>
#endif
#ifdef USE_HDF5
#include <hdf5.h>
#include "hdf4.h"
#include "extern_hdf5.h"
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "performance.h"
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

// function prototypes

#ifdef USE_HDF5
int HDF5_WriteStringAttr(hid_t dset_id, char *Alabel, char *String, 
			 FILE *log_fptr);
#endif


void WriteListOfFloats(FILE *fptr, int N, float floats[]);
int  DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid, 
				      TopGridData *MetaData, 
				      LevelHierarchyEntry *LevelArray[], 
				      int level);

#define NUMBER_OF_PROJECTED_FIELDS 9




int ProjectToPlane(TopGridData &MetaData, LevelHierarchyEntry *LevelArray[],
		   int ProjectStartTemp[], int ProjectEndTemp[], 
		   FLOAT ProjectStartCoordinate[],
		   FLOAT ProjectEndCoordinate[], int ProjectLevel,
		   int ProjectionDimension, char *ProjectionFileName,
		   int ProjectionSmooth, ExternalBoundary *Exterior)
{

  int i, j, dim, field, level, ret, size = 1;
  int ProjectDim[MAX_DIMENSION];
  float *ProjectedField[NUMBER_OF_PROJECTED_FIELDS], TempCellWidth;
  FLOAT ProjectLeft[MAX_DIMENSION], ProjectRight[MAX_DIMENSION];

#ifdef USE_HDF5
  FILE *log_fptr;

  char *dset_name;

  hid_t       file_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_type_id, FILE_type_id;
  hid_t       file_dsp_id;

  hsize_t     OutDims[2];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int         io_log = 1;

  int ii = sizeof(float32);

  switch(ii)
  {

    case 4:
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
      break;

    case 8:
      float_type_id = HDF5_R8;
      file_type_id = HDF5_FILE_R8;
      break;

    default:
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;

  }

  int jj = sizeof(FLOAT);

  switch(jj)
  {

    case 4:
      FLOAT_type_id = HDF5_R4;
      FILE_type_id = HDF5_FILE_R4;
      break;

    case 8:
      FLOAT_type_id = HDF5_R8;
      FILE_type_id = HDF5_FILE_R8;
      break;

    case 16:
      FLOAT_type_id = HDF5_R16;
      FILE_type_id = H5Tcopy(HDF5_FILE_B8);
                     H5Tset_size(FILE_type_id,16);
      break;

    default:
      printf("INCORRECT FLOAT DEFINITION\n");

  }
#endif /* USE_HDF5 */

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
      sscanf(line, "XrayLowerCutoffkeV = %"FSYM, &XrayLowerCutoffkeV);
      sscanf(line, "XrayUpperCutoffkeV = %"FSYM, &XrayUpperCutoffkeV);
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

      /* Set particle density. */

      if (SelfGravity) {
	CopyOverlappingParticleMassFields(Temp->GridData, &MetaData, 
					  LevelArray, level);
	if (Temp->GridHierarchyEntry->ParentGrid != NULL)
	  Temp->GridHierarchyEntry->ParentGrid->GridData->DepositParticlePositions(Temp->GridData, Temp->GridHierarchyEntry->ParentGrid->GridData->ReturnTime(), GRAVITATING_MASS_FIELD_PARTICLES);
      }

      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */

      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);

//#ifdef UNUSED

      /* Set boundary (interpolate from parent is good enough). */

      if (level > 0) {
	Temp->GridData->InterpolateBoundaryFromParent
	                   (Temp->GridHierarchyEntry->ParentGrid->GridData);
	CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, level);
      } else
	Temp->GridData->SetExternalBoundaryValues(Exterior);

      /* Set old parent field for children's interpolation (and set the time
	 an arbitrary factor (20%) ahead so interpolation doesn't error). */

      if (level <= ProjectLevel) {
	Temp->GridData->CopyBaryonFieldToOldBaryonField();
	Temp->GridData->SetTime(Temp->GridData->ReturnTime()*1.2);
      }

//#endif /* UNUSED */

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

#ifdef UNUSED
      /* Repair the damage done by zeroing and reprojecting quantities. */

      Temp2 = LevelArray[level+1];
      if (level < ProjectLevel)
	while (Temp2 != NULL) {
	  Temp2->GridData->ProjectSolutionToParentGrid(*Temp->GridData);
	  Temp2 = Temp2->NextGridThisLevel;
	}
#endif /* UNUSED */

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

  /* Write out the projected planes. */

  if (debug) printf("Writing projected planes to %s.\n", ProjectionFileName);

#ifdef USE_HDF5
  char *logname = new char[strlen(ProjectionFileName)+5];
  strcpy(logname, ProjectionFileName);
  strcat(logname, ".log");
  log_fptr = fopen(logname, "a");

  if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", ProjectionFileName);

  file_id = H5Fcreate(ProjectionFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (io_log) fprintf(log_fptr, "H5Fcreate id: %d\n", file_id);
    assert( file_id != h5_error );
#endif /* USE_HDF5 */

  /* Set dimensions (reversed since this is c and we're using f77 order). */

#ifdef USE_HDF4
  int32 OutDims[2];
#endif

  for (dim = 0, i = 0; dim < MetaData.TopGridRank; dim++)
    if (dim != ProjectionDimension)
      OutDims[1 - i++] = ProjectDim[dim];

#ifdef USE_HDF4
  if (DFSDsetdims(2, OutDims) == HDF_FAIL) {
    fprintf(stderr, "Error in DFSDsetdims.\n");
    return FAIL;
  }
#endif

  /* Compute the scales. */

  float32 *float_temp;

  for (dim = 0, i = 0; dim < MetaData.TopGridRank; dim++)
    if (dim != ProjectionDimension) {

      float_temp = new float32[ProjectDim[dim]];

      for (j = 0; j < ProjectDim[dim]; j++)
	float_temp[j] = float32(ProjectLeft[dim] + 
	  (float(j)+0.5) * (ProjectRight[dim]-ProjectLeft[dim]) / 
	    float(ProjectDim[dim]));

#ifdef USE_HDF4
      if (DFSDsetdimscale(2 - i++, ProjectDim[dim], (VOIDP) float_temp) == 
	  HDF_FAIL) {
	fprintf(stderr, "Error in DFSDsetdimscale.\n");
	return FAIL;
      }
#endif

      delete [] float_temp;
    }

  /* Write projected planes to HDF file. */

  float_temp = new float32[size];

  for (i = 0; i < NUMBER_OF_PROJECTED_FIELDS; i++) {

    /* Set name of field. */


#ifdef USE_HDF4
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
      DFSDsetdatastrs("OVII column density", "cm^-2", "f8.3","");
#endif

    /* for i == 3, divide by x-ray luminosity. */

    if (i == 3)
      for (j = 0; j < size; j++)
        ProjectedField[3][j] /= ProjectedField[1][j];

    /* for i == 7, divide by baryon density. */

    if (i == 7)
      for (j = 0; j < size; j++)
        ProjectedField[7][j] /= ProjectedField[0][j];

    /* Copy to float32 space. */

    for (j = 0; j < size; j++)
      float_temp[j] = float32(ProjectedField[i][j]);

    /* Write data. */

#ifdef USE_HDF4
    if (i == 0)
      ret = DFSDputdata(ProjectionFileName, 2, OutDims, (VOIDP) float_temp);
    if (i != 0)
      ret = DFSDadddata(ProjectionFileName, 2, OutDims, (VOIDP) float_temp);

    if (ret == HDF_FAIL) {
      fprintf(stderr, "Error in DFSDput/adddata.\n");
      return FAIL;
    }
#endif

#ifdef USE_HDF5
    file_dsp_id = H5Screate_simple(2, OutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %d\n", file_dsp_id);
      assert( file_dsp_id != h5_error );

    switch(i)
    {
      
      case 0:
        dset_name = "projected_gas_density";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "projected_gas_density", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "Msolar/Mpc^2", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 1:
        dset_name = "projected_x-ray_luminosity_div1e23";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "projected_x-ray_luminosity_div1e23", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "none", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 2:
        dset_name = "projected_dm_density";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "projected_dm_density", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "Msolar/Mpc^2", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 3:
        dset_name = "projected_x-ray_weighted_temperature";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "projected_x-ray_weighted_temperature", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "K", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 4:
        dset_name = "projected_level";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "projected_level", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "none", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "f8.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 5:
        dset_name = "SZ_y_effect";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "SZ y effect", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "none", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "f8.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 6:
        dset_name = "DT_over_T_Doppler_effect";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "DT/T Doppler effect", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "none", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "f8.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 7:
        dset_name = "Metallicity";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "Metallicity", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "none", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "f8.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 8:
        dset_name = "projected_star_density";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "projected_star_density", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "none", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "f8.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      case 9:
        dset_name = "OVII column density";

        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);

        dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %d\n", dset_id);
          assert( dset_id != h5_error );

        HDF5_WriteStringAttr(dset_id, "Label", "OVII column density", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Units", "cm^-2", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Format", "f8.3", log_fptr);
        HDF5_WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        break;

      default:
        fprintf(stderr, "No such case in ProjectToPlane\n");

    }


    h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) float_temp);
      if (io_log) fprintf(log_fptr, "H5Dwrite: %d\n", h5_status);
      assert( h5_status != h5_error );
      JBPERF_COUNT_WRITE(dset_id, float_type_id, H5S_ALL);

    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %d\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %d\n", h5_status);
      assert( h5_status != h5_error );
#endif /* USE_HDF5 */

  }
  delete [] float_temp;


#ifdef USE_HDF5
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log_fptr, "H5Fclose: %d\n", h5_status);
    assert( h5_status != h5_error );

  fclose(log_fptr);
#endif

  return SUCCESS;
}

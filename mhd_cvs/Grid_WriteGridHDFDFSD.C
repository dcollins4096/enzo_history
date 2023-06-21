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
/  GRID CLASS (WRITE OUT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/


//  Write grid to file pointer fptr
//     (we assume that the grid is at an appropriate stopping point,
//      where the Old values aren't required)
//

#ifdef USE_HDF4

#include <string.h>
#include <stdio.h>
#include <df.h>
#include "macros_and_parameters.h"

#ifdef WRITEGRID_HDF_DFSD

#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);

int grid::WriteGridHDFDFSD(FILE *fptr, char *base_name, int grid_id)
{

  /* declarations */

  int i, j, k, dim, field, ret, size, ActiveDim[MAX_DIMENSION];
  int32 OutDims[MAX_DIMENSION];
  float32 *temp;

  char *ParticlePositionLabel[] = 
     {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleVelocityLabel[] = 
     {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
				    "metallicity_fraction", "alpha_fraction"};

  /* initialize */

  char id[10];
  sprintf(id, "%4.4d", grid_id);

  /* make sure quantities defined at least for 3d */

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;

  /* ------------------------------------------------------------------- */
  /* 1) Save general grid class data */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(fptr, "GridRank          = %d\n", GridRank);

    fprintf(fptr, "GridDimension     = ");
    WriteListOfInts(fptr, GridRank, GridDimension);

    fprintf(fptr, "GridStartIndex    = ");
    WriteListOfInts(fptr, GridRank, GridStartIndex);

    fprintf(fptr, "GridEndIndex      = ");
    WriteListOfInts(fptr, GridRank, GridEndIndex);

    fprintf(fptr, "GridLeftEdge      = ");
    WriteListOfFloats(fptr, GridRank, GridLeftEdge);

    fprintf(fptr, "GridRightEdge     = ");
    WriteListOfFloats(fptr, GridRank, GridRightEdge);

    fprintf(fptr, "Time              = %"GOUTSYM"\n", Time);

    fprintf(fptr, "SubgridsAreStatic = %d\n", SubgridsAreStatic);
  
    fprintf(fptr, "NumberOfBaryonFields = %d\n", NumberOfBaryonFields);

  }

  /* ------------------------------------------------------------------- */
  /* 2) save baryon field quantities (including fields). */

  char *name = new char[strlen(base_name)+strlen(id)+1];
  strcpy(name, base_name);
  strcat(name, id);

  if (NumberOfBaryonFields > 0) {

    if (MyProcessorNumber == ROOT_PROCESSOR) {

      fprintf(fptr, "FieldType = ");
      WriteListOfInts(fptr, NumberOfBaryonFields, FieldType);

      fprintf(fptr, "BaryonFileName = %s\n", name);

      fprintf(fptr, "CourantSafetyNumber    = %f\n", CourantSafetyNumber);
      fprintf(fptr, "PPMFlatteningParameter = %d\n", PPMFlatteningParameter);
      fprintf(fptr, "PPMDiffusionParameter  = %d\n", PPMDiffusionParameter);
      fprintf(fptr, "PPMSteepeningParameter = %d\n", PPMSteepeningParameter);

    }

    if (MyProcessorNumber == ProcessorNumber) {

    /* 2a) Set HDF file dimensions (use FORTRAN ordering). */

    for (dim = 0; dim < GridRank; dim++)
      OutDims[GridRank-dim-1] = ActiveDim[dim];

    if (DFSDsetdims(GridRank, OutDims) == HDF_FAIL) {
      fprintf(stderr, "Error in DFSDsetdims.\n");
      return FAIL;
    }

    /* 2b) Write out co-ordinate values.  Use the centre of each cell. */

    size = 1;
    for (dim = GridRank-1; dim >= 0; dim--) {

      /* DFSDsetdimstrs(GridRank-dim, DimLabels[dim], DimUnits[dim],"e10.4");*/
      /* This doesn't work -- why not? */

      /* Compute cell centers and put them in temp. */

      temp = new float32[GridDimension[dim]];
      for (i = 0; i <= GridEndIndex[dim] - GridStartIndex[dim]; i++)
	temp[i] = CellLeftEdge[dim][GridStartIndex[dim] + i] +
	          0.5 * CellWidth[dim][GridStartIndex[dim] + i];
      if (DFSDsetdimscale(GridRank-dim, ActiveDim[dim], (VOIDP) temp) 
	                                                        == HDF_FAIL)
	fprintf(stderr, "WriteGrid: Error in DFSDsetdimscale.\n");
      delete [] temp;

      size *= GridDimension[dim];
    }

    /* create temporary buffer */

    temp = new float32[size];

    /* 2c) Loop over fields, writing each one. */

    for (field = 0; field < NumberOfBaryonFields; field++) {

      /* copy active part of field into grid */

      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           + 
	         (j-GridStartIndex[1])*ActiveDim[0]              + 
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		       float32(
	      BaryonField[field][i + j*GridDimension[0] +
		                     k*GridDimension[0]*GridDimension[1]]
                              );

      /* set datafield name and units, etc. */
      
      if (DFSDsetdatastrs(DataLabel[field], DataUnits[field], "e10.4",
			  "Cartesian") == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDsetdatastrs.\n");
	return FAIL;
      }

      /* create a new file if this is the first field, otherwise, add it
	 to the existing file. */
      
      if (field == 0)
	ret = DFSDputdata(name, GridRank, OutDims, (VOIDP) temp);
      else
	ret = DFSDadddata(name, GridRank, OutDims, (VOIDP) temp);
      if (ret == HDF_FAIL) {
	fprintf(stderr, "Error in DFSD.\n");
	return FAIL;
      }
    }   // end of loop over fields

    /* If this is cosmology, compute the temperature field as well since
       its such a pain to compute after the fact. */

    if (ComovingCoordinates) {

      /* Allocate field and compute temperature. */

      float *temperature = new float[size];
      if (this->ComputeTemperatureField(temperature) == FAIL) {
	fprintf(stderr, "Error in grid->ComputeTemperatureField.\n");
	return FAIL;
      }

      /* Copy active part of field into grid */

      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           + 
	         (j-GridStartIndex[1])*ActiveDim[0]              + 
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     float32(
		   temperature[(k*GridDimension[1] + j)*GridDimension[0] + i]
			     );

      /* set datafield name and units, etc. */
      
      if (DFSDsetdatastrs("Temperature", "K", "e10.4", "Cartesian") 
	  == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDsetdatastrs (temperature).\n");
	return FAIL;
      }

      if (DFSDadddata(name, GridRank, OutDims, (VOIDP) temp) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDadddata (temperature).\n");
	return FAIL;
      }

      delete [] temperature;

    }

    /* Make sure that there is a copy of dark matter field to save
       (and at the right resolution). */

    if (SelfGravity && NumberOfParticles > 0) {
      float SaveGravityResolution = GravityResolution;
      GravityResolution = 1;
      this->InitializeGravitatingMassFieldParticles(RefineBy);
      this->ClearGravitatingMassFieldParticles();
      this->DepositParticlePositions(this, Time, 
				     GRAVITATING_MASS_FIELD_PARTICLES);
      GravityResolution = SaveGravityResolution;
    }

    /* If present, write out the GravitatingMassFieldParticles. */

    if (GravitatingMassFieldParticles != NULL) {

      /* Set dimensions. */

      int StartIndex[MAX_DIMENSION], EndIndex[MAX_DIMENSION];
      for (dim = 0; dim < 3; dim++) {
	StartIndex[dim] = nint((GridLeftEdge[dim] - 
				GravitatingMassFieldParticlesLeftEdge[dim])/
			       GravitatingMassFieldParticlesCellSize);
	EndIndex[dim] = nint((GridRightEdge[dim] - 
			      GravitatingMassFieldParticlesLeftEdge[dim])/
			     GravitatingMassFieldParticlesCellSize) - 1;
      }

      /* Copy active part of field into grid */

      for (k = StartIndex[2]; k <= EndIndex[2]; k++)
	for (j = StartIndex[1]; j <= EndIndex[1]; j++)
	  for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	    temp[(i-StartIndex[0])                           + 
	         (j-StartIndex[1])*ActiveDim[0]              + 
	         (k-StartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     float32(
			     GravitatingMassFieldParticles[ i + 
			       j*GravitatingMassFieldParticlesDimension[0] +
			       k*GravitatingMassFieldParticlesDimension[0]*
			         GravitatingMassFieldParticlesDimension[1]]
			     );

      /* set datafield name and units, etc. */
      
      if (DFSDsetdatastrs("Dark matter density", "", "e10.4", "Cartesian") 
	  == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDsetdatastrs (dark matter).\n");
	return FAIL;
      }

      if (DFSDadddata(name, GridRank, OutDims, (VOIDP) temp) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDadddata (dark matter).\n");
	return FAIL;
      }

      /* Clean up if we modified the resolution. */

      if (SelfGravity && GravityResolution != 1)
	this->DeleteGravitatingMassFieldParticles();

    } // end of (if GravitatingMassFieldParticles != NULL)

    delete [] temp;
    
    /* Write BoundaryFluxes info (why? it's just recreated when the grid
                                  is read in) */

   } // end: if (ProcessorNumber == MyProcessorNumber)
  } // end: if (NumberOfBaryonFields > 0)

  /* ------------------------------------------------------------------- */
  /* 3) Save particle quantities. */
  
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "NumberOfParticles   = %d\n", NumberOfParticles);

  if (NumberOfParticles > 0) {

    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(fptr, "ParticleFileName = %s\n", name);

    if (MyProcessorNumber == ProcessorNumber) {

    /* Sort particles according to their identifier. */

    this->SortParticlesByNumber();
    
    /* Clear any HDF parameters set previously. */

    DFSDclear();

    /* Create a temporary buffer (32 bit). */

    float32 *temp = new float32[NumberOfParticles];
    int32 TempIntArray[1];
    TempIntArray[0] = int32(NumberOfParticles);

    /* Copy particle positions to temp and write them. */

    if (sizeof(float) != sizeof(FLOAT)) {
      fprintf(stderr, "DFSD write grid only supports r4.\n");
      return FAIL;
    }

    for (dim = 0; dim < GridRank; dim++) {

      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticlePosition[dim][i]);

      DFSDsetdims(1, TempIntArray);
      DFSDsetdatastrs(ParticlePositionLabel[dim], "", "", "");
      if (DFSDadddata(name, 1, TempIntArray, (VOIDP) temp) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDadddata (particle position).\n");
	return FAIL;
      }

    }

    /* Copy particle velocities to temp and write them. */

    for (dim = 0; dim < GridRank; dim++) {

      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleVelocity[dim][i]);

      DFSDsetdatastrs(ParticleVelocityLabel[dim], "", "", "");
      if (DFSDadddata(name, 1, TempIntArray, (VOIDP) temp) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDadddata (particle position).\n");
	return FAIL;
      }

    }

    /* Copy mass to temp and write it. */

    for (i = 0; i < NumberOfParticles; i++)
      temp[i] = float32(ParticleMass[i]);

    DFSDsetdatastrs("particle mass", "", "", "");
    if (DFSDadddata(name, 1, TempIntArray, (VOIDP) temp) == HDF_FAIL) {
      fprintf(stderr, "Error in DFSDadddata (particles).\n");
      return FAIL;
    }

    /* Copy number (ID) to temp and write it. */

    int32 *tempint = new int32[NumberOfParticles];
    for (i = 0; i < NumberOfParticles; i++)
      tempint[i] = int32(ParticleNumber[i]);

    DFSDsetNT(DFNT_INT32);
    DFSDsetdatastrs("particle_index", "", "", "");
    if (DFSDadddata(name, 1, TempIntArray, (VOIDP) tempint) == HDF_FAIL) {
      fprintf(stderr, "Error in DFSDadddata (particles).\n");
      return FAIL;
    }
    DFSDsetNT(DFNT_FLOAT32);

    /* Copy particle attributes to temp and write them. */

    for (j = 0; j < NumberOfParticleAttributes; j++) {

      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleAttribute[j][i]);

      DFSDsetdatastrs(ParticleAttributeLabel[j], "", "", "");
      if (DFSDadddata(name, 1, TempIntArray, (VOIDP) temp) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDadddata (particle attribute).\n");
	return FAIL;
      }

    }

    /* clean up */

    delete [] temp;
    delete [] tempint;

  } // end: if (MyProcessorNumber...)
  } // end: if (NumberOfParticles > 0)

  /* 4) Save Gravity info. */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    if (SelfGravity)
      fprintf(fptr, "GravityBoundaryType = %d\n", GravityBoundaryType);

  /* Clean up. */
  
  delete [] name;
  
  return SUCCESS;

}

#endif /* WRITEGRID_HDF_DFSD */
#endif /* USE_HDF4 */


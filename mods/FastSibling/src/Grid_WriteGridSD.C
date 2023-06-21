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

#include <string.h>
#include <stdio.h>
#include "mfhdf.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#ifdef WRITEGRID_HDF_SD

#define NO_WRITE_DIM_SCALES

/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);

int grid::WriteGrid(FILE *fptr, char *base_name, int grid_id)
{

  /* declarations */

  int i, j, k, dim, field, size, ActiveDim[MAX_DIMENSION];
  int32 OutDims[MAX_DIMENSION], sd_id, sds_id, start[] = {0,0,0};
#ifdef WRITE_DIM_SCALES
  int32 dim_id;
#endif /* WRITE_DIM_SCALES */
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

  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }
  for (dim = 0; dim < 3; dim++)
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

  char *name = new char[strlen(base_name)+strlen(id)+1];
  strcpy(name, base_name);
  strcat(name, id);

  /* Open HDF file for writing. */

  if (MyProcessorNumber == ProcessorNumber)
    if ((sd_id = SDstart(name, DFACC_CREATE)) == HDF_FAIL) {
      fprintf(stderr, "Error opening file %s.\n", name);
      return FAIL;
    }

  /* ------------------------------------------------------------------- */
  /* 2) save baryon field quantities (including fields). */

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

    /* 2b) Write out co-ordinate values.  Use the centre of each cell. */

    size = 1;
    float32 *tempdim[MAX_DIMENSION];
    for (dim = GridRank-1; dim >= 0; dim--) {
      
      /* Compute cell centers and put them in temp. */

      tempdim[dim] = new float32[GridDimension[dim]];
      for (i = 0; i <= GridEndIndex[dim] - GridStartIndex[dim]; i++)
	tempdim[dim][i] = CellLeftEdge[dim][GridStartIndex[dim] + i] +
	          0.5 * CellWidth[dim][GridStartIndex[dim] + i];
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

      /* Create a new SDS. */

      if ((sds_id = SDcreate(sd_id, DataLabel[field], DFNT_FLOAT32, GridRank,
			     OutDims)) == HDF_FAIL) {
	fprintf(stderr, "Error in SDCreate for %s\n", DataLabel[field]);
	return FAIL;
      }

      /* set datafield name and units, etc. */
      
      if (SDsetdatastrs(sds_id, DataLabel[field], DataUnits[field], "e10.4",
			"Cartesian") == HDF_FAIL) {
	fprintf(stderr, "Error in SDsetdatastrs.\n");
	return FAIL;
      }

      /* For each exis, add the dimscale. */

#ifdef WRITE_DIM_SCALES
      for (dim = 0; dim < GridRank; dim++) {
	if ((dim_id = SDgetdimid(sds_id, GridRank-dim-1)) == HDF_FAIL) {
	  fprintf(stderr, "Error in SDgetdimid (%d).\n", dim);
	  return FAIL;
	}
	if (SDsetdimscale(dim_id, ActiveDim[dim], DFNT_FLOAT32, 
			  (VOIDP) tempdim[dim]) == HDF_FAIL) {
	  fprintf(stderr, "Error in SDsetdimscale (dim %d)\n", dim);
	  return FAIL;
	}
      }
#endif /* WRITE_DIM_SCALES */

      /* write data and clean up access to this SDS. */
      
      if (SDwritedata(sds_id, start, NULL, OutDims, (VOIDP)temp) == HDF_FAIL) {
	fprintf(stderr, "Error in SDwritedata (field = %d).\n", field);
	return FAIL;
      }

      if (SDendaccess(sds_id) == HDF_FAIL) {
	fprintf(stderr, "Error in SDendaccess.\n");
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

      /* output */

      sds_id = SDcreate(sd_id, "Temperature", DFNT_FLOAT32, GridRank, OutDims);
      SDsetdatastrs(sds_id, "Temperature", "K", "e10.4", "Cartesian");
#ifdef WRITE_DIM_SCALES
      for (dim = 0; dim < GridRank; dim++) {
	dim_id = SDgetdimid(sds_id, GridRank-dim);
	SDsetdimscale(dim_id, ActiveDim[dim], DFNT_FLOAT32, 
		      (VOIDP) tempdim[dim]);
      }
#endif /* WRITE_DIM_SCALES */
      if (SDwritedata(sds_id, start, NULL, OutDims, (VOIDP)temp) == HDF_FAIL) {
	fprintf(stderr, "Error writing temperature\n");
	return FAIL;
      }
      SDendaccess(sds_id);

      delete temperature;

    } // end: if (ComovingCoordinates)

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

      int StartIndex[] = {0,0,0}, EndIndex[] = {0,0,0};
      for (dim = 0; dim < GridRank; dim++) {
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

      /* output */

      sds_id = SDcreate(sd_id, "Dark matter density", DFNT_FLOAT32, GridRank, 
			OutDims);
      SDsetdatastrs(sds_id, "Dark matter density", "", "e10.4", "Cartesian");
#ifdef WRITE_DIM_SCALES
      for (dim = 0; dim < GridRank; dim++) {
	dim_id = SDgetdimid(sds_id, GridRank-dim);
	SDsetdimscale(dim_id, ActiveDim[dim], DFNT_FLOAT32, 
		      (VOIDP) tempdim[dim]);
      }
#endif /* WRITE_DIM_SCALES */
      if (SDwritedata(sds_id, start, NULL, OutDims, (VOIDP)temp) == HDF_FAIL) {
	fprintf(stderr, "Error writing Dark matter density\n");
	return FAIL;
      }
      SDendaccess(sds_id);

      /* Clean up if we modified the resolution. */

      if (SelfGravity && GravityResolution != 1)
	this->DeleteGravitatingMassFieldParticles();

    } // end of (if GravitatingMassFieldParticles != NULL)

    delete temp;
    for (dim = 0; dim < GridRank; dim++)
      delete tempdim[dim];
    
    /* Write BoundaryFluxes info (why? it's just recreated when the grid
                                  is read in) */

   }  // end: if (ProcessorNumber == MyProcessorNumber)
  } // end: if (NumberOfBaryonFields > 0)

  /* ------------------------------------------------------------------- */
  /* 3) Save particle quantities. */
  
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "NumberOfParticles   = %d\n", NumberOfParticles);

  if (NumberOfParticles > 0) {

    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(fptr, "ParticleFileName = %s\n", name); // must be same as above

    if (MyProcessorNumber == ProcessorNumber) {

    /* Sort particles according to their identifier. */

    this->SortParticlesByNumber();

    /* Create a temporary buffer (32 bit). */

    float32 *temp = new float32[NumberOfParticles];
    int32 TempIntArray[1];

    /* Particle positions are not converted to 32 bit first. 
       (128 bit numbers are not supported by HDF so convert to 64). */

    int32 HDFDataType = (sizeof(FLOAT) == 4) ? DFNT_FLOAT32 : DFNT_FLOAT64;
    float64 *temp_pointer = NULL;
    if (sizeof(FLOAT) == 16)
      temp_pointer = new float64[NumberOfParticles];
    TempIntArray[0] = int32(NumberOfParticles);

    for (dim = 0; dim < GridRank; dim++) {

      /* Convert to 64 if 128, either just write out. */

      if (sizeof(FLOAT) == 16)
	for (i = 0; i < NumberOfParticles; i++)
	  temp_pointer[i] = float64(ParticlePosition[dim][i]);
      else
	temp_pointer = (float64*) ParticlePosition[dim];

      sds_id = SDcreate(sd_id, ParticlePositionLabel[dim], HDFDataType,
			1, TempIntArray);
      if (SDwritedata(sds_id, start, NULL, TempIntArray, 
		      (VOIDP) temp_pointer) == HDF_FAIL) {
	fprintf(stderr, "Error in SDwritedata (particle position).\n");
	return FAIL;
      }
      SDendaccess(sds_id);

    }

    if (sizeof(FLOAT) == 16)
      delete [] temp_pointer;  /* clean up if allocated. */

    /* Copy particle velocities to temp and write them. */

    for (dim = 0; dim < GridRank; dim++) {

      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleVelocity[dim][i]);

      sds_id = SDcreate(sd_id, ParticleVelocityLabel[dim], DFNT_FLOAT32, 
			1, TempIntArray);
      if (SDwritedata(sds_id, start, NULL, TempIntArray, (VOIDP) temp) == 
	  HDF_FAIL) {
	fprintf(stderr, "Error in SDwritedata (particle velocity).\n");
	return FAIL;
      }
      SDendaccess(sds_id);

    }

    /* Copy mass to temp and write it. */

    for (i = 0; i < NumberOfParticles; i++)
      temp[i] = float32(ParticleMass[i]);

    sds_id = SDcreate(sd_id, "particle mass", DFNT_FLOAT32, 1, TempIntArray);
    if (SDwritedata(sds_id, start, NULL, TempIntArray, (VOIDP) temp) == 
	HDF_FAIL) {
      fprintf(stderr, "Error in SDwritedata (particle mass).\n");
      return FAIL;
    }
    SDendaccess(sds_id);

    /* Copy number (ID) to temp and write it. */

    int32 *tempint = new int32[NumberOfParticles];
    for (i = 0; i < NumberOfParticles; i++)
      tempint[i] = int32(ParticleNumber[i]);

    sds_id = SDcreate(sd_id, "particle_index", DFNT_INT32, 1, TempIntArray);
    if (SDwritedata(sds_id, start, NULL, TempIntArray, (VOIDP) tempint) == 
	HDF_FAIL) {
      fprintf(stderr, "Error in SDwritedata (particle index).\n");
      return FAIL;
    }
    SDendaccess(sds_id);

    /* Copy type to temp and write it. */

    if (ParticleType == NULL)
      return FAIL;
    for (i = 0; i < NumberOfParticles; i++)
      tempint[i] = int32(ParticleType[i]);

    sds_id = SDcreate(sd_id, "particle_type", DFNT_INT32, 1, TempIntArray);
    if (SDwritedata(sds_id, start, NULL, TempIntArray, (VOIDP) tempint) == 
	HDF_FAIL) {
      fprintf(stderr, "Error in SDwritedata (particle index).\n");
      return FAIL;
    }
    SDendaccess(sds_id);

    /* Copy particle attributes to temp and write them. */

    for (j = 0; j < NumberOfParticleAttributes; j++) {

      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleAttribute[j][i]);

      sds_id = SDcreate(sd_id, ParticleAttributeLabel[j], DFNT_FLOAT32, 
			1, TempIntArray);
      if (SDwritedata(sds_id, start, NULL, TempIntArray, (VOIDP) temp) == 
	  HDF_FAIL) {
	fprintf(stderr, "Error in SDwritedata (particle attribute).\n");
	return FAIL;
      }
      SDendaccess(sds_id);

    }
    /* clean up */

    delete temp;
    delete tempint;

  } // end: if (MyProcessorNumber...)
  } // end: if (NumberOfParticles > 0)

  /* Close HDF file. */

  if (MyProcessorNumber == ProcessorNumber)
    SDend(sd_id);
  
  /* 4) Save Gravity info. */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    if (SelfGravity)
      fprintf(fptr, "GravityBoundaryType = %d\n", GravityBoundaryType);

  /* Clean up. */
  
  delete name;
  
  return SUCCESS;

}
#endif /* WRITEGRID_HDF_SD */

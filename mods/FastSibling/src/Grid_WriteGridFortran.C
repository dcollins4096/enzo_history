/***********************************************************************
/
/  GRID CLASS (WRITE OUT GRID)
/
/  written by: Greg Bryan
/  date:       October, 2000
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
#include "df.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#ifdef WRITEGRID_FORTRAN

/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);

extern "C" void FORTRAN_NAME(fortopen)(char *name, int *unit);
extern "C" void FORTRAN_NAME(fortwrite)(int *unit, float32 *temp, int *rank, 
					int *shape);
extern "C" void FORTRAN_NAME(fortiwrite)(int *unit, int *temp, int *rank, 
					 int *shape);
extern "C" void FORTRAN_NAME(fortclose)(int *unit);

int grid::WriteGrid(FILE *fptr, char *base_name, int grid_id)
{

  /* declarations */

  int i, j, k, dim, field, size, ActiveDim[MAX_DIMENSION];
  float32 *temp;
  int One = 1, FortranUnit = 2;  /* default */

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

  char *name = new char[strlen(base_name)+strlen(id)+6];
  strcpy(name, base_name);
  strcat(name, ".fort");
  strcat(name, id);

  /* Open file for writing. */

  if (MyProcessorNumber == ProcessorNumber)
    FORTRAN_NAME(fortopen)(name, &FortranUnit);

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
      
    /* 2b) Write out co-ordinate values.  Use the centre of each cell. */

    size = 1;
    for (dim = GridRank-1; dim >= 0; dim--)
      size *= GridEndIndex[dim] - GridStartIndex[dim] + 1;

    /* create temporary buffer */

    temp = new float32[size];

    /* Write out four data items indicating dimensions and number
       of particles. */

    i = 3;
    FORTRAN_NAME(fortiwrite)(&FortranUnit, ActiveDim, &One, &i);
    FORTRAN_NAME(fortiwrite)(&FortranUnit, &NumberOfParticles, &One, &One);

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

      /* write data. */

      FORTRAN_NAME(fortwrite)(&FortranUnit, temp, &GridRank, ActiveDim);

    }   // end of loop over fields

    /* If this is cosmology, compute the temperature field as well since
       its such a pain to compute after the fact. */

    if (ComovingCoordinates) {

      /* Allocate field and compute temperature. */

      float *temperature = new float[GridDimension[0]*GridDimension[1]*
				     GridDimension[2]];
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

      FORTRAN_NAME(fortwrite)(&FortranUnit, temperature, &GridRank, ActiveDim);

      delete [] temperature;

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

      FORTRAN_NAME(fortwrite)(&FortranUnit, temp, &GridRank, ActiveDim);

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
      fprintf(fptr, "ParticleFileName = %s\n", name); // must be same as above

    if (MyProcessorNumber == ProcessorNumber) {

    /* Sort particles according to their identifier. */

    this->SortParticlesByNumber();

    /* Create a temporary buffer (32 bit). */

    float32 *temp = new float32[NumberOfParticles];

    /* Copy particle positions to temp and write them. */

    for (dim = 0; dim < GridRank; dim++) {

      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticlePosition[dim][i]);

      FORTRAN_NAME(fortwrite)(&FortranUnit, temp, &One, &NumberOfParticles);

    }

    /* Copy particle velocities to temp and write them. */

    for (dim = 0; dim < GridRank; dim++) {

      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleVelocity[dim][i]);

      FORTRAN_NAME(fortwrite)(&FortranUnit, temp, &One, &NumberOfParticles);

    }

    /* Copy mass to temp and write it. */

    for (i = 0; i < NumberOfParticles; i++)
      temp[i] = float32(ParticleMass[i]);

    FORTRAN_NAME(fortwrite)(&FortranUnit, temp, &One, &NumberOfParticles);

    /* Copy number (ID) to temp and write it. */

    int32 *tempint = new int32[NumberOfParticles];
    for (i = 0; i < NumberOfParticles; i++)
      tempint[i] = int32(ParticleNumber[i]);

    FORTRAN_NAME(fortiwrite)(&FortranUnit, tempint, &One, &NumberOfParticles);

    /* Copy type to temp and write it. */

    int32 *tempint = new int32[NumberOfParticles];
    for (i = 0; i < NumberOfParticles; i++)
      tempint[i] = int32(ParticleNumber[i]);

    FORTRAN_NAME(fortiwrite)(&FortranUnit, tempint, &One, &NumberOfParticles);

    /* Copy particle attributes to temp and write them. */

    for (j = 0; j < NumberOfParticleAttributes; j++) {
      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleAttribute[j][i]);

      FORTRAN_NAME(fortwrite)(&FortranUnit, temp, &One, &NumberOfParticles);
    }

    /* clean up */

    delete [] temp;
    delete [] tempint;

  } // end: if (MyProcessorNumber...)
  } // end: if (NumberOfParticles > 0)

  /* Close file. */

  if (MyProcessorNumber == ProcessorNumber)
    FORTRAN_NAME(fortclose)(&FortranUnit);
  
  /* 4) Save Gravity info. */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    if (SelfGravity)
      fprintf(fptr, "GravityBoundaryType = %d\n", GravityBoundaryType);

  /* Clean up. */
  
  delete [] name;
  
  return SUCCESS;

}
#endif /* WRITEGRID_FORTRAN */

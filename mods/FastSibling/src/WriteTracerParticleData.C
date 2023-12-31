/***********************************************************************
/
/  WRITE OUT GRID DATA AT THE LOCATION OF ANY TRACER PARTICLES
/
/  written by: Greg Bryan
/  date:       March, 2004
/  modified1:
/
/  PURPOSE: This function walks through the grids and writes out,
/     for each grid, the position, density and temperature at the
/     location of tracer particles which are advected at the grid
/     velocity.  These follow stream-lines through the gas.
/     Each processor creates its own file to eliminate communication.
/
************************************************************************/


#include <string.h>
#include <stdio.h>
#include <df.h>
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
#include "CosmologyParameters.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);

static char SummarySuffix[] = ".summary";

int WriteTracerParticleData(char *basename, int dumpnumber, 
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData, 
		   FLOAT WriteTime)
{
  /* Exit if tracer particles not turned on. */

  if (TracerParticleOn == FALSE)
    return SUCCESS;

  /* declarations */

  char id[7], name[MAX_LINE_LENGTH];
  FILE *fptr = NULL;
  int level, Zero = 0;

  /* Compute redshift and units. */

  FLOAT a = 1, dadt, Redshift = 0;
  float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1, 
        TemperatureUnits;
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(WriteTime, &a, &dadt);
    CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		      &TimeUnits, &VelocityUnits, WriteTime);
    Redshift = (1.0+InitialRedshift)/a - 1;
  }

  /* Create output filename */

  strcpy(name, basename);
  sprintf(id, "%-d", MyProcessorNumber);  /* create processor # */
  strcat(name, id);

  if (debug)
    printf("WriteTracerParticleData: writing file %s.\n", name);

  /* Open file. */

  if ((fptr = fopen(name, "ab")) == NULL) {
    fprintf(stderr, "Error opening output file %s\n", name);
    return FAIL;
  }

  /* Output header information to the file :
        0) time in problem units
        1) redshift
        2) time conversion factor
        3) length conversion factor
        4) density conversion factor
        5) velocity conversion factor 
        6) Number of values per tracer particle
           (currently 5: xpos,ypos,zpos,dens,temp) */

  float float_temp = float(WriteTime);  // convert from FLOAT to float
  fwrite((void*) &float_temp, sizeof(float), 1, fptr);
  float_temp = float(Redshift);
  fwrite((void*) &float_temp, sizeof(float), 1, fptr);
  fwrite((void*) &TimeUnits, sizeof(float), 1, fptr);
  fwrite((void*) &LengthUnits, sizeof(float), 1, fptr);
  fwrite((void*) &DensityUnits, sizeof(float), 1, fptr);
  fwrite((void*) &VelocityUnits, sizeof(float), 1, fptr);
  int int_temp = 5;
  fwrite((void*) &int_temp, sizeof(int), 1, fptr);

  /* --------------------------------------------------------------- */
  /* Loop over grids and write grid data to files. */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    /* Loop over all the grids. */

    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {

      /* Write out grid info (also deletes the under subgrid field). */

      if (Temp->GridData->TracerParticleOutputData(fptr, WriteTime) == FAIL) {
	fprintf(stderr, "Error in grid->OutputTracerParticleData.\n");
	return FAIL;
      }

      /* Next grid on this level. */

      Temp = Temp->NextGridThisLevel;

    } // end loop over grids

  } // end loop over levels

  /* Write a 0 to indicate no more data at this time. */

  fwrite((void*) &Zero, sizeof(int), 1, fptr);

  /* close file. */

  fclose(fptr);

  return SUCCESS;
}

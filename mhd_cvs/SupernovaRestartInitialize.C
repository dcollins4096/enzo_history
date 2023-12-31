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
/  INITIALIZE A SUPERNOVA EXPLOSION FROM A RESTART CALCULATION
/
/  written by: Greg Bryan
/  date:       February, 2000
/  modified1:  Brian O'Shea, September 2003.  
/                  added metallicity parameter.
/
/  PURPOSE:  This routine reads in a previously generated output file
/            (presumably from a cosmology calculation), and then
/            initializes a supernova explosion in the center.
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
#include "CosmologyParameters.h"
#include "fortran.def"

/* function prototypes */

void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);



int SupernovaRestartInitialize(FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid, TopGridData &MetaData,
			       ExternalBoundary &Exterior)
{

  /* declarations */

  char line[MAX_LINE_LENGTH];
  int dim, level, ret;

  /* Set default supernova parameters. */

  float SupernovaRestartEjectaMass   = 1.0;   // in solar masses
  float SupernovaRestartMetalMass    = 0.0;   // in solar masses
  float SupernovaRestartEjectaRadius = 1.0;   // in pc
  float SupernovaRestartEjectaEnergy = 1.0;   // in 10^51 erg
  FLOAT SupernovaRestartEjectaCenter[MAX_DIMENSION];
  int   SupernovaRestartColourField   = FALSE;

  char *SupernovaRestartName = NULL;

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    SupernovaRestartEjectaCenter[dim] = FLOAT_UNDEFINED;

  /* Error check. */

  /* Read input from file. */

  char *dummy = new char[MAX_LINE_LENGTH];

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* Read parameters */

    ret += sscanf(line, "SupernovaRestartEjectaMass = %"FSYM, 
		  &SupernovaRestartEjectaMass);
    ret += sscanf(line, "SupernovaRestartMetalMass = %"FSYM, 
		  &SupernovaRestartMetalMass);
    ret += sscanf(line, "SupernovaRestartEjectaRadius = %"FSYM, 
		  &SupernovaRestartEjectaRadius);
    ret += sscanf(line, "SupernovaRestartEjectaEnergy = %"FSYM, 
		  &SupernovaRestartEjectaEnergy);
    ret += sscanf(line, "SupernovaRestartEjectaCenter = %"PSYM" %"PSYM" %"PSYM,
		  SupernovaRestartEjectaCenter,
		  SupernovaRestartEjectaCenter+1,
		  SupernovaRestartEjectaCenter+2);
    ret += sscanf(line, "SupernovaRestartColourField = %d",
		  &SupernovaRestartColourField);
    
    if (sscanf(line, "SupernovaRestartName = %s", dummy) == 1) {
      SupernovaRestartName = dummy;
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "SupernovaRestart") &&
	line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }

  /* More error checking. */

  if (SupernovaRestartName == NULL) {
    fprintf(stderr, "Missing restart file name.\n");
    return FAIL;
  }

  /* -------------------------------------------------------------------- */
  /* Read the restart file. */

  if (debug)
    printf("reading restart parameter file %s\n", SupernovaRestartName);
  if (ReadAllData(SupernovaRestartName, &TopGrid, MetaData, &Exterior) 
      == FAIL) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Error in ParameterFile %s.\n", SupernovaRestartName);
    return FAIL;
  }
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stderr, "Successfully read restart file %s.\n", 
	    SupernovaRestartName);

  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

  /* Convert Mass, Radius and Energy into code units of density, radius
     and thermal energy (by assuming mass and energy is evenly distributed
     within radius). 
     If comoving coordinate is set then assume mass is in solar units,
     radius in pc, and energy in 10^51 erg; otherwise do no conversion.*/


  double MassConversion = 1, LengthConversion = 1, EnergyConversion = 1;
  float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
        TemperatureUnits = 1;

  if (ComovingCoordinates) {

    FLOAT Time = TopGrid.GridData->ReturnTime();

    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }

    LengthConversion = 3.08e18;     // pc
    MassConversion   = 2e33;        // solar masses
    EnergyConversion = 1.0e51;      // 10^51 erg

  }

  float EjectaRadius = SupernovaRestartEjectaRadius * LengthConversion;
  float EjectaDensity = SupernovaRestartEjectaMass * MassConversion/ 
                        (4.0/3.0*3.14159*POW(EjectaRadius, 3));
  float EjectaMetalDensity = SupernovaRestartMetalMass * MassConversion/ 
                        (4.0/3.0*3.14159*POW(EjectaRadius, 3));
  float EjectaThermalEnergy = SupernovaRestartEjectaEnergy * EnergyConversion /
        (SupernovaRestartEjectaMass * MassConversion);

  EjectaRadius        /= LengthUnits;
  EjectaDensity       /= DensityUnits;
  EjectaMetalDensity  /= DensityUnits;
  EjectaThermalEnergy /= VelocityUnits*VelocityUnits;

  if (debug) {
    printf("SupernovaRestart: initial T = %g K\n", 
	   EjectaThermalEnergy*TemperatureUnits*(Gamma-1.0)*0.6);
    printf("SupernovaRestart: r (code units) = %g\n", EjectaRadius);
    printf("SupernovaRestart: density (code units) = %g\n", EjectaDensity);
  }

  /* -------------------------------------------------------------------- */
  /* Loop over all the grids and call the initializer to modify them
     if necessary. */

  LevelHierarchyEntry *Temp, *Temp2;
  int NumberOfCellsSet = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    Temp = LevelArray[level];

    while (Temp != NULL) {

      if (Temp->GridData->SupernovaRestartInitialize(EjectaDensity,
			       EjectaMetalDensity,
			       EjectaRadius, EjectaThermalEnergy,
			       SupernovaRestartEjectaCenter,
			       SupernovaRestartColourField, 
			       &NumberOfCellsSet) == FAIL) {
	fprintf(stderr, "Error in grid->SupernovaRestartInitialize\n");
	return FAIL;
      }
      Temp = Temp->NextGridThisLevel;
    }

  }
  if (debug)
    printf("SupernovaRestart: NumberOfCellsSet = %d\n", NumberOfCellsSet);

  /* -------------------------------------------------------------------- */
  /* Loop over grid and project solution to parent to maintain consistency. */

  for (level = MaximumRefinementLevel; level > 0; level--) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      if (Temp->GridData->ProjectSolutionToParentGrid(
                *(Temp->GridHierarchyEntry->ParentGrid->GridData)) == FAIL) {
	fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	return FAIL;
      }
      Temp2 = Temp->NextGridThisLevel;
      delete Temp;   // clean up as we go along
      Temp = Temp2;
    }
  }

  /* -------------------------------------------------------------------- */
  /* Write parameters to parameter output file */

  
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "SupernovaRestartEjectaMass   = %f\n", 
	    SupernovaRestartEjectaMass);
    fprintf(Outfptr, "SupernovaRestartMetalMass   = %f\n", 
	    SupernovaRestartMetalMass);
    fprintf(Outfptr, "SupernovaRestartEjectaRadius = %f\n", 
	    SupernovaRestartEjectaRadius);
    fprintf(Outfptr, "SupernovaRestartEjectaEnergy = %f\n", 
	    SupernovaRestartEjectaEnergy);
    fprintf(Outfptr, "SupernovaRestartEjectaCenter = ");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, 
		      SupernovaRestartEjectaCenter);
    fprintf(Outfptr, "SupernovaRestartColourField  = %d\n", 
	    SupernovaRestartColourField);
    fprintf(Outfptr, "SupernovaRestartName         = %s\n", 
	    SupernovaRestartName);

  }

  /* Clean up. */

  delete dummy;

  return SUCCESS;
}

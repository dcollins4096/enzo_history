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
/  COMPUTES A SET OF "GLOBAL" STATISTICS
/
/  written by: Greg Bryan
/  date:       May, 2000
/  modified1:  Robert Harkness, Jul 2002 (bugfixes)
/
/  PURPOSE:
/
************************************************************************/

//
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#undef DEFINE_STORAGE

/* function prototypes */

int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int  DepositParticleMassField(HierarchyEntry *Grid, float RequestTime = -1.0);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid, 
				      TopGridData *MetaData, 
				      LevelHierarchyEntry *LevelArray[], 
				      int level);
int CommunicationInitialize(int *argc, char **argv[]);
int CommunicationFinalize();
void my_exit(int status);

main(int argc, char *argv[])
{
  CommunicationInitialize(&argc, &argv);

  /* Main declarations */

  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

  /* Initialize */

  debug                        = FALSE;
  char *ParameterFile          = NULL;
  char *myname                 = argv[0];
  float RegionLeft[MAX_DIMENSION], RegionRight[MAX_DIMENSION];
  int level, dim, i, j;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    RegionLeft[dim] = RegionRight[dim] = FLOAT_UNDEFINED;

  /* --------------------------------------------------------------- */
  /* Interpret command-line arguments. */

  char c;
  while (--argc > 0 && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {

	/* get beginning of region selection (float coordinate). */

      case 'b':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"FSYM, &RegionLeft[dim++]) != 1) {
	    fprintf(stderr, "%s: error reading Begin coordinates.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* debug */

      case 'd':
	debug = TRUE;
	break;

	/* get finish of region selection (float coordinate). */

      case 'f':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"FSYM, &RegionRight[dim++]) != 1) {
	    fprintf(stderr, "%s: error reading Finish coordinates.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* Unknown */

      default:
	fprintf(stderr, "%s: unknown command-line option: -%s.\n", myname, &c);
	my_exit(EXIT_FAILURE);
	
      } // end of switch(c)

  /* Error check for number of parameters, and set parameter file. */

  if (argc != 1) {
    fprintf(stderr, "enzohop [-b #] [-f #] [-d] amr_file\n");
    fprintf(stderr, "  -b)egin region\n");
    fprintf(stderr, "  -f)inish region\n");
    fprintf(stderr, "  -d)ebug\n");
    my_exit(FAIL);
  }
  ParameterFile = argv[0];

  /* Read the saved file. */

  SetDefaultGlobalValues(MetaData);
  if (ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
    fprintf(stderr, "Error in ParameterFile %s.\n", ParameterFile);
    my_exit(EXIT_FAILURE);
  }
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

  /* Set top grid boundary conditions. */

  if (TopGrid.GridData->SetExternalBoundaryValues(&Exterior) == FAIL) {
    fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
    my_exit(EXIT_FAILURE);
  }

  /* Setup temperature distribution. */

  int NumberOfTempBins = 3, NumberOfDensityBins = 100;
  float DensBinStart = -3.0, DensBinEnd = +4.0;
  float *TempBinEdge = new float[NumberOfTempBins+1];
  double *TempBinInfo = new double[(NumberOfTempBins+1)*3];
  double *TempDensInfo = new double[NumberOfTempBins*NumberOfDensityBins*2];

  TempBinEdge[0] = 1.0;
  TempBinEdge[1] = 1.0e5;
  TempBinEdge[2] = 1.0e7;
  TempBinEdge[3] = 1.0e10;

  for (i = 0; i < NumberOfTempBins*3+3; i++)
    TempBinInfo[i] = 0;
  for (i = 0; i < NumberOfTempBins*NumberOfDensityBins*2; i++)
    TempDensInfo[i] = 0;

  /* If undefined, set parameters. */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (RegionLeft[dim] == FLOAT_UNDEFINED)
      RegionLeft[dim] = DomainLeftEdge[dim];
    if (RegionRight[dim] == FLOAT_UNDEFINED)
      RegionRight[dim] = DomainRightEdge[dim];
  }

  if (debug)
    printf("Analysis region: Left = %g %g %g   Right = %g %g %g\n",
	   RegionLeft[0], RegionLeft[1], RegionLeft[2],
	   RegionRight[0], RegionRight[1], RegionRight[2]);

  /* ------------------------------------------------------------ */
  /* Zero the solution (on this grid) which is underneath any subgrid
     (so we get only the high resolution solution from the subgrid). */

  LevelHierarchyEntry *Temp, *Temp2;
  if (debug) printf("->zeroing redundant solution & setting BCs\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
      if (level > 0) {
	Temp->GridHierarchyEntry->ParentGrid->GridData->SetTime(
				Temp->GridData->ReturnTime());
	Temp->GridData->InterpolateBoundaryFromParent
	                   (Temp->GridHierarchyEntry->ParentGrid->GridData);
	CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, level);
      }
      Temp2 = LevelArray[level+1];
      while (Temp2 != NULL) {
	Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						 ZERO_UNDER_SUBGRID_FIELD);
	Temp2 = Temp2->NextGridThisLevel;
      }
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* --------------------------------------------------------------- */
  /* Loop over all the levels and compute stats */

  if (debug) printf("->computing stats\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    /* If SelfGravity, set all the particle mass fields. */

#ifdef UNUSED
    Temp = LevelArray[level];
    if (SelfGravity)
      while (Temp != NULL) {
	DepositParticleMassField(Temp->GridHierarchyEntry);
	Temp = Temp->NextGridThisLevel;
      }
#endif /* UNUSED */

    /* Loop over all the grids. */

    Temp = LevelArray[level];
    if (debug) printf("level = %d\n", level);
    while (Temp != NULL) {

      /* Set particle density. */

#ifdef UNUSED
      if (SelfGravity) {
	CopyOverlappingParticleMassFields(Temp->GridData, &MetaData, 
					  LevelArray, level);
	if (Temp->GridHierarchyEntry->ParentGrid != NULL)
	  Temp->GridHierarchyEntry->ParentGrid->GridData->DepositParticlePositions(Temp->GridData, Temp->GridHierarchyEntry->ParentGrid->GridData->ReturnTime(), GRAVITATING_MASS_FIELD_PARTICLES);
      }
#endif /* UNUSED */

      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */

#ifdef UNUSED
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);

      /* Zero the solution (on this grid) which is underneath any subgrid
	 (so we get only the high resolution solution from the subgrid). */

      Temp2 = LevelArray[level+1];
      while (Temp2 != NULL) {
	Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
        					 ZERO_UNDER_SUBGRID_FIELD);
	Temp2 = Temp2->NextGridThisLevel;
      }
#endif /* UNUSED */

      /* Compute stats. */

      Temp->GridData->ComputeGlobalStatistics(RegionLeft, RegionRight,
			  NumberOfTempBins, TempBinEdge, TempBinInfo,
			  NumberOfDensityBins, DensBinStart, DensBinEnd,
				 TempDensInfo);

      /* Delete grid. */

      delete Temp->GridData;

      /* Next grid on this level. */

      Temp = Temp->NextGridThisLevel;

    } // end loop over grids

  } // end loop over levels

  /* --------------------------------------------------------------- */

  /* Output results. */

  printf("temp range:    f_mass   f_volume   <d^2>/<d>^2\n");
  for (i = 0; i < NumberOfTempBins; i++)
    printf("%g-%g K:  %g   %g  %g\n", TempBinEdge[i], TempBinEdge[i+1],
	   TempBinInfo[i*3+3]/TempBinInfo[0],
	   TempBinInfo[i*3+4]/TempBinInfo[1],
	   TempBinInfo[i*3+5]*TempBinInfo[i*3+4]/POW(TempBinInfo[i*3+3],2) );
  printf("Total mass, Volume = %g  %g\n", TempBinInfo[0], TempBinInfo[1]);

  /* Now open file and output density distribution (big). */

  FILE *fptr;
  if ((fptr = fopen("density_distribution.out", "w")) == NULL) {
    fprintf(stderr, "Error opening output file.\n");
    my_exit(EXIT_FAILURE);
  }

  fprintf(fptr, "# 1st column for each temp range is density weighted\n");
  fprintf(fptr, "# 2nd column for each temp range is volume weighted\n");
  fprintf(fptr, "#density  ");
  for (i = 0; i < NumberOfTempBins; i++)
    fprintf(fptr, "     %g-%g    ", TempBinEdge[i], TempBinEdge[i+1]);
  fprintf(fptr, "\n");

  for (j = 0; j < NumberOfDensityBins; j++) {
    fprintf(fptr, "%g  ", POW(10.0, DensBinStart + 
      (float(j)+0.5)*(DensBinEnd - DensBinStart)/float(NumberOfDensityBins)));
    for (i = 0; i < NumberOfTempBins; i++)
      fprintf(fptr, "  %g %g  ", TempDensInfo[i*NumberOfDensityBins+j],
	      TempDensInfo[(i+NumberOfTempBins)*NumberOfDensityBins+j]);
    fprintf(fptr, "\n");
  }

  fclose(fptr);

  my_exit(EXIT_SUCCESS);
}

void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}

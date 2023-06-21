/***********************************************************************
/
/  FINDS DENSITY PEAKS IN THE DISTRIBUTION
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
 
#define MAX_PARTICLE_NUMBER 2000000
 
/* function prototypes */
 
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
#ifdef UNUSED
int  ComputeAccelerationField(HierarchyEntry &Grid, TopGridData *MetaData,
			      LevelHierarchyEntry *LevelArray[], int level,
			      ExternalBoundary *Exterior = NULL);
#endif
int  DepositParticleMassField(HierarchyEntry *Grid, FLOAT RequestTime = -1.0);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid,
				       TopGridData *MetaData,
				       LevelHierarchyEntry *LevelArray[],
				       int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
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
 
  int level;
 
  /* Initialize */
 
  debug                = TRUE;
  char *myname         = argv[0];
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
 
  int   NumberOfPeaks       = 1;
  float PeakSeparation      = FLOAT_UNDEFINED;
  float MinDensity          = 100.0;
  int   PeakField           = 0;
 
  /* Interpret command-line arguments */
 
  int i, dim;
  char c;
 
  while (--argc > 0 && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {
 
	/* n) Number of Peaks. */
 
      case 'n':
	if (--argc > 0) {
	  if (sscanf((*++argv), "%"ISYM, &NumberOfPeaks) != 1) {
	    fprintf(stderr, "%s: error reading NumberOfPeaks.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  fprintf(stderr, "%s: Need to specify NumberOfPeaks.\n", myname);
	  my_exit(EXIT_FAILURE);
	}
	break;
 
	/* m) min density. */
 
      case 'm':
	if (--argc > 0) {
	  if (sscanf((*++argv), "%"FSYM, &MinDensity) != 1) {
	    fprintf(stderr, "%s: error reading MinDensity.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  fprintf(stderr, "%s: Need to specify MinDensity.\n", myname);
	  my_exit(EXIT_FAILURE);
	}
	break;
 
	/* s) separation. */
 
      case 's':
	if (--argc > 0) {
	  if (sscanf((*++argv), "%"FSYM, &PeakSeparation) != 1) {
	    fprintf(stderr, "%s: error reading PeakSeparation.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  fprintf(stderr, "%s: Need to specify PeakSeparation.\n", myname);
	  my_exit(EXIT_FAILURE);
	}
	break;
 
	/* b) baryon density. */
 
      case 'b':
	PeakField = 0;
	break;
 
	/* d) dark matter density. */
 
      case 'd':
	PeakField = 1;
	break;
 
	/* p) potential field. */
 
      case 'p':
	PeakField = 2;
	break;
 
	/* x) x-ray. */
 
      case 'x':
	PeakField = 3;
	break;
 
	/* Unknown */
 
      default:
	fprintf(stderr, "%s: unknown command-line option: -%s\n", myname, &c);
	my_exit(EXIT_FAILURE);
	
      } // end of switch(c)
 
  /* Error check */
 
  if (argc != 2) {
    fprintf(stderr, "usage: %s [options] amr_file cluster_file\n", myname);
    fprintf(stderr, "  options: -n # (number of peaks)\n");
    fprintf(stderr, "           -m # (minimum density)\n");
    fprintf(stderr, "           -s # (min peak separation, comoving Mpc)\n");
    fprintf(stderr, "           -b   (use baryon density in peak seach, default)\n");
    fprintf(stderr, "           -d   (use dark matter density in peak seach)\n");
    fprintf(stderr, "           -p   (use potential in peak search)\n");
    fprintf(stderr, "           -x   (use bolometric x-ray in peak seach)\n");
    my_exit(EXIT_FAILURE);
  }
 
  /* Read the first hierarchy. */
 
  SetDefaultGlobalValues(MetaData);
  printf("Reading input %s\n", argv[0]);
  if (ReadAllData(argv[0], &TopGrid, MetaData, &Exterior) == FAIL) {
    fprintf(stderr, "Error in ParameterFile %s.\n", argv[0]);
    my_exit(EXIT_FAILURE);
  }
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels
 
  /* Allocate space for the peak list. */
 
  float *Position[MAX_DIMENSION], *Radius = new float[NumberOfPeaks];
  float *PeakValue = new float[NumberOfPeaks];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Position[dim] = new float[NumberOfPeaks];
 
  /* Convert Peak separation from comoving Mpc to box units. */
 
  float BoxSize = 1;
  if (ComovingCoordinates)
    BoxSize = ComovingBoxSize/HubbleConstantNow;
  if (PeakSeparation == FLOAT_UNDEFINED)
    PeakSeparation = 0.01*BoxSize;
  PeakSeparation *= (DomainRightEdge[0]-DomainLeftEdge[0])/BoxSize;
 
  /* Just set the radius to the minimum peak separation. */
 
  for (i = 0; i < NumberOfPeaks; i++)
    Radius[i] = PeakSeparation;
 
  /* ------------------------------------------------------------------- */
  /* Find N largest peaks with a separation greater than that specified. */
 
  /* Loop over levels. */
 
  if (debug) fprintf(stderr, "FindingPeaks");
  int NumberOfPeaksFound = 0;
  PeakValue[0] = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    LevelHierarchyEntry *Temp = LevelArray[level];
 
    if (PeakField == 1 || PeakField == 2)
      while (Temp != NULL) {
	if (DepositParticleMassField(Temp->GridHierarchyEntry) == FAIL) {
	  fprintf(stderr, "Error in DepositParticleMassField.\n");
	  my_exit(EXIT_FAILURE);
	}
	Temp = Temp->NextGridThisLevel;
      }
 
    Temp = LevelArray[level];
 
    /* Loop over the grids on this level. */
 
    while (Temp != NULL) {
      fprintf(stderr, ".");
 
      /* Set boundary. */
 
      if (level > 0) {
	Temp->GridData->InterpolateBoundaryFromParent
	  (Temp->GridHierarchyEntry->ParentGrid->GridData, LevelArray[level]);
	CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, level);
      }
 
      /* Set old parent field for children's interpolation (and set the time
	 an arbitrary factor (1%) ahead so interpolation doesn't error). */
 
      Temp->GridData->CopyBaryonFieldToOldBaryonField();
      Temp->GridData->SetTime(Temp->GridData->ReturnTime()*1.01);
 
      /* Prepare dark-matter field if needed. */
 
      if (PeakField == 1) {
 
	if (Temp->GridData->InitializeGravitatingMassField(RefineBy)
	    == FAIL) {
	  fprintf(stderr, "Error in grid->InitializeGravitatingMassField.\n");
	  my_exit(EXIT_FAILURE);
	}
	Temp->GridData->ClearGravitatingMassField();
 
	if (CopyOverlappingParticleMassFields(Temp->GridData, &MetaData,
					      LevelArray, level) == FAIL) {
	  fprintf(stderr, "Error in CopyOverlappingParticleMassFields.\n");
	  my_exit(EXIT_FAILURE);
	}
 
      } // end: if (PeakField == 1)
 
      /* Prepare potential field if needed. */
 
      if (PeakField == 2) {
	fprintf(stderr, "Potential field not implemented.\n");
	my_exit(EXIT_FAILURE);
#ifdef UNUSED
	ComputePotential = TRUE;
	if (LevelArray[level+1] == NULL)
	  LevelArray[level+1] = LevelArray[0];
	if (ComputeAccelerationField(*(Temp->GridHierarchyEntry), &MetaData,
				     LevelArray, level, &Exterior) == FAIL) {
	  fprintf(stderr, "Error in ComputeAccelerationField.\n");
	  my_exit(EXIT_FAILURE);
	  }
	if (LevelArray[level+1] == LevelArray[0])
	  LevelArray[level+1] = NULL;
#endif /* UNUSED */
      }
 
      /* Find peaks. */
 
      if (Temp->GridData->FindDensityPeaks(NumberOfPeaks, NumberOfPeaksFound,
		   Position, PeakValue, PeakSeparation, MinDensity, PeakField)
	  == FAIL) {
	fprintf(stderr, "Error in grid->FindDensityPeaks.\n");
	exit(EXIT_FAILURE);
      }
 
      /* Clean up. */
 
      if (PeakField == 1)
	Temp->GridData->DeleteGravitatingMassField();
 
      if (PeakField == 2) {
	//	Temp->GridData->DeleteAccelerationFieldForCells();
	Temp->GridData->DeleteParticleAcceleration();
      }
 
      Temp = Temp->NextGridThisLevel;
 
    } // end: loop over grids on this level
 
  } // end: loop over levels
  fprintf(stderr, "\n");
 
  /* ------------------------------------------------------------ */
  /* Write the peak file. */
 
  if (debug)
    printf("NumberOfPeaksFound = %"ISYM"\n", NumberOfPeaksFound);
  FILE *fptr;
  if ((fptr = fopen(argv[1], "w")) == NULL) {
    fprintf(stderr, "Error opening output_file %s\n", argv[1]);
    my_exit(EXIT_FAILURE);
  }
  for (i = 0; i < NumberOfPeaksFound; i++)
    fprintf(fptr, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n", Position[0][i], Position[1][i],
	    Position[2][i], Radius[i], PeakValue[i]);
  fclose(fptr);
 
  my_exit(EXIT_SUCCESS);
}
 
void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}

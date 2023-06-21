/***********************************************************************
/
/  OUTPUT STAR PARTICLES TO THEIR OWN FILE
/
/  written by: Greg Bryan
/  date:       May, 1999
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
#include <string.h>
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
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int CommunicationInitialize(int *argc, char **argv[]);
int CommunicationFinalize();
void my_exit(int status);

static char DefaultDumpName[] = "DumpGridData";

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

  debug                = FALSE;
  char *myname         = argv[0];
  char *ParameterFile = NULL, *BaseName  = NULL;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;

  /* general declarations. */

  int WriteGrid = FALSE, WriteDM = FALSE, WriteStars = FALSE;
  FILE *gridfptr = NULL, *dmfptr = NULL, *starfptr = NULL;
  char GridName[MAX_LINE_LENGTH], DMName[MAX_LINE_LENGTH], 
    StarName[MAX_LINE_LENGTH];
  LevelHierarchyEntry *Temp2, *Temp;

  /* --------------------------------------------------------------- */
  /* Interpret command-line arguments. */

  char c;
  while (--argc > 0 && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {

	/* debug */

      case 'd':
	debug = TRUE;
	break;

	/* grid data */

      case 'g':
	WriteGrid = 1;
	break;

	/* dm particle data */

      case 'p':
	WriteDM = 1;
	break;

	/* star data */

      case 's':
	WriteStars = 1;
	break;

	/* Unknown */

      default:
	fprintf(stderr, "%s: unknown command-line option: -%s.\n", myname, &c);
	my_exit(EXIT_FAILURE);
	
      } // end of switch(c)

  /* Error check for number of parameters, and set parameter file. */

  if (argc < 1 || argc > 2) {
    fprintf(stderr, "usage: %s [-d] [-g] [-p] [-s] amr_saved_filename [dump_filename]\n", myname);
    fprintf(stderr, "  -d)ebug\n");
    fprintf(stderr, "  -g)rid data output\n");
    fprintf(stderr, "  -p)article data output\n");
    fprintf(stderr, "  -s)tar data output\n");
    my_exit(FAIL);
  }
  ParameterFile = argv[0];
  if (argc == 2)
    BaseName = argv[1];
  else
    BaseName = DefaultDumpName;

  /* If none selected, set default. */

  if (WriteStars+WriteDM+WriteGrid == 0)
    WriteGrid = 1;

  /* Read the saved file. */

  SetDefaultGlobalValues(MetaData);
  printf("Reading saved file %s\n", ParameterFile);
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

  /* ------------------------------------------------------------ */
  /* Zero the solution (on this grid) which is underneath any subgrid
     (so we get only the high resolution solution from the subgrid). */

  if (debug) printf("->zeroing redundant solution & setting BCs\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
      if (level > 0) {
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

  /* ------------------------------------------------------------ */
  /* Output the grid data, dark matter particles and star particles. */

  /* Generate names and open files for output */

  strcpy(GridName, BaseName);
  strcat(GridName, ".grid");
  strcpy(DMName, BaseName);
  strcat(DMName, ".dm");
  strcpy(StarName, BaseName);
  strcat(StarName, ".star");

  if (WriteGrid == TRUE) {
    printf("Writing grid file %s\n", GridName);
    if ((gridfptr = fopen(GridName, "w")) == NULL) {
      fprintf(stderr, "error opening %s\n", GridName);
      my_exit(EXIT_FAILURE);
    }
  }

  if (WriteDM == TRUE) {
    printf("Writing dark matter file %s\n", DMName);
    if ((dmfptr = fopen(DMName, "w")) == NULL) {
      fprintf(stderr, "error opening %s\n", DMName);
      my_exit(EXIT_FAILURE);
    }
  }

  if (StarParticleCreation > 0 && WriteStars == TRUE) {
    printf("Writing star particle file %s\n", StarName);
    if ((starfptr = fopen(StarName, "w")) == NULL) {
      fprintf(stderr, "error opening %s\n", StarName);
      my_exit(EXIT_FAILURE);
    }
  }

  /* Loop over levels. */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    LevelHierarchyEntry *Temp = LevelArray[level];

    while (Temp != NULL) {
      
      /* Prepare the dark matter field (includes sub-grids, but not
	 ghost regions). */

      if (SelfGravity)
	if (DepositParticleMassField(Temp->GridHierarchyEntry,
				     Temp->GridData->ReturnTime()) == FAIL) {
	  fprintf(stderr, "Error in DepositParticleMassField.\n");
	  my_exit(EXIT_FAILURE);
	}  

      /* Dump grid data */

      if (Temp->GridData->DumpGridData(gridfptr, dmfptr, starfptr) == FAIL) {
	fprintf(stderr, "Error in grid->DumpGridData.\n");
	my_exit(EXIT_FAILURE);
      }

      Temp = Temp->NextGridThisLevel;

    } // end: loop over grids on this level

  } // end: loop over levels

  /* close files. */

  if (gridfptr != NULL)
    fclose(gridfptr);
  if (dmfptr != NULL)
    fclose(dmfptr);
  if (starfptr != NULL)
    fclose(starfptr);

  my_exit(EXIT_SUCCESS);
}

void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}

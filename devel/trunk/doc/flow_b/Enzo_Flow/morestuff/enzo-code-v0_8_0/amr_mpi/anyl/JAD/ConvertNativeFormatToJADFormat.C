/***********************************************************************
/
/  CONVERT FROM OUR NATIVE FORMAT TO JOHN'S AMR DATA (JAD) FORMAT
/
/  written by: Greg Bryan
/  date:       April, 1996
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
#include "IO.hh"
#include "IEEEIO.hh"
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

  debug                = TRUE;
  char *myname         = argv[0];
  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;

  /* Error check */

  if (argc != 3) {
    fprintf(stderr, "usage: %s amr_saved_filename jad_filename\n", myname);
    my_exit(EXIT_FAILURE);
  }

  /* Read the saved file. */

  SetDefaultGlobalValues(MetaData);
  printf("Reading saved file %s\n", argv[1]);
  if (ReadAllData(argv[1], &TopGrid, MetaData, &Exterior) == FAIL) {
    fprintf(stderr, "Error in ParameterFile %s.\n", argv[1]);
    my_exit(EXIT_FAILURE);
  }
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

  /* Find the lowest level. */

  int FinestLevel;
  for (FinestLevel = 0; FinestLevel < MAX_DEPTH_OF_HIERARCHY-1; FinestLevel++)
    if (LevelArray[FinestLevel+1] == NULL)
      break;

  /* ------------------------------------------------------------ */

  /* Write the grids in JAD format. */

  printf("Writing JAD file %s\n", argv[2]);

  IObase *filehandle = new IEEEIO(argv[2], IObase::Create);
  if (!filehandle->isValid()) {
    fprintf(stderr, "error opening %s\n", argv[2]);
    my_exit(EXIT_FAILURE);
  }

  float TopGridWidth = (DomainRightEdge[0]-DomainLeftEdge[0])/
    float(MetaData.TopGridDims[0]);

  /* Loop over levels. */

  int GridID = 1, Resolution = 1;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    LevelHierarchyEntry *Temp = LevelArray[level];

    while (Temp != NULL) {

      if (Temp->GridData->WriteGridJADFormat(filehandle, TopGridWidth,
			  level, GridID++, Resolution, FinestLevel) == FAIL) {
	fprintf(stderr, "Error in grid->WriteGridJADFormat.\n");
	my_exit(EXIT_FAILURE);
      }

      Temp = Temp->NextGridThisLevel;

    } // end: loop over grids on this level

    Resolution *= RefineBy;

  } // end: loop over levels

  /* Close the SDS file. */

  delete filehandle;

  my_exit(EXIT_SUCCESS);
}

void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}

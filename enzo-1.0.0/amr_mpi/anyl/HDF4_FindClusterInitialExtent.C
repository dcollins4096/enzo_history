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
/  FINDS AN OBJECT'S INITIAL PARTICLE POSITIONS
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
#include "df.h"
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
int SetDefaultGlobalValues(TopGridData &MetaData);
void QuickSortAndDrag(int List[], int left, int right,
                      int NumberToDrag, float *DragList1[],
                      int NumberToDrag2, FLOAT *DragList2[]);
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

  /* Error check */

  if (argc < 3 || argc > 4) {
    fprintf(stderr, "usage: %s amr_final_output amr_initial_output cluster_file\n", myname);
    fprintf(stderr, "   or: %s amr_output cluster_file\n", myname);
    fprintf(stderr, "  (cluster_file format: x y z r   in box units)\n");
    my_exit(EXIT_FAILURE);
  }

  /* Read the cluster file. */

  FILE *fptr;
  float Position[MAX_DIMENSION], Radius;
  if ((fptr = fopen(argv[argc-1], "r")) == NULL) {
    fprintf(stderr, "Error opening cluster_file %s\n", argv[argc-1]);
    my_exit(EXIT_FAILURE);
  }
  if (fscanf(fptr, "%"FSYM" %"FSYM" %"FSYM" %"FSYM, Position, Position+1, 
	     Position+2, &Radius) != 4) {
    fprintf(stderr, "Error reading cluster_file %s\n", argv[argc-1]);
    my_exit(EXIT_FAILURE);
  }
  fclose(fptr);

  /* Read the first hierarchy. */

  SetDefaultGlobalValues(MetaData);
  printf("Reading input1 %s\n", argv[1]);
  if (ReadAllData(argv[1], &TopGrid, MetaData, &Exterior) == FAIL) {
    fprintf(stderr, "Error in ParameterFile %s.\n", argv[1]);
    my_exit(EXIT_FAILURE);
  }
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

  /* ------------------------------------------------------------ */
  /* Find the particles within the specified radius of the given position. */

  /* Allocate a big buffer. */

  int *ParticleNumberList = new int[MAX_PARTICLE_NUMBER];
  int ParticlesFound = 0;
  LevelHierarchyEntry *Temp, *Temp2;

  /* Loop over levels. */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    Temp = LevelArray[level];

    while (Temp != NULL) {

      /* Find & put particles in ParticleNumberList, seting ParticlesFound. */

      if (Temp->GridData->FindParticlesWithinSphere(Position, Radius,
			       &ParticlesFound, ParticleNumberList) == FAIL) {
	fprintf(stderr, "Error in grid->FindParticlesWithinSphere.\n");
	my_exit(EXIT_FAILURE);
      }

      if (ParticlesFound > MAX_PARTICLE_NUMBER) {
	fprintf(stderr, "Increase MAX_PARTICLE_NUMBER\n");
	my_exit(EXIT_FAILURE);
      }

      /* Delete grid after checking. */

      delete Temp->GridData;
      Temp2 = Temp->NextGridThisLevel;
      delete Temp;

      Temp = Temp2;

    } // end: loop over grids on this level

    LevelArray[level] = NULL; // clean up

  } // end: loop over levels
  printf("NumberOfParticlesInRegion = %d\n", ParticlesFound);

  /* ------------------------------------------------------------ */
  /* If there is no second hierarchy, then just output the
     particle indices. */

  if (argc == 3) {

    int32 TempInt = ParticlesFound;

    DFSDsetdims(1, &TempInt);
    DFSDsetdatastrs("particle_index", "", "", "");
    if (DFSDputdata("findinit.hdf", 1, &TempInt, (VOIDP) ParticleNumberList) ==
	HDF_FAIL) {
      fprintf(stderr, "Error writint findinit.hdf\n");
      my_exit(EXIT_FAILURE);
    }

    my_exit(EXIT_SUCCESS);
  }

  /* ------------------------------------------------------------ */
  /* Find those particles in the second hierarchy. */

  int dim;
  float LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    LeftEdge[dim] = DomainRightEdge[dim];
    RightEdge[dim] = DomainLeftEdge[dim];
  }

  /* Open output file. */

  if ((fptr = fopen("FindClusterInitialExtent.out", "w")) == NULL) {
    fprintf(stderr, "Error opening output file.\n");
    my_exit(EXIT_FAILURE);
  }

  /* Sort list of particles. */

  QuickSortAndDrag(ParticleNumberList, 0, ParticlesFound-1, 0, NULL, 0, NULL);

  /* Read the second hierarchy. */

  printf("Reading input2 %s\n", argv[2]);
  if (ReadAllData(argv[2], &TopGrid, MetaData, &Exterior) == FAIL) {
    fprintf(stderr, "Error in ParameterFile %s.\n", argv[1]);
    my_exit(EXIT_FAILURE);
  }
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

  /* Loop over levels. */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    Temp = LevelArray[level];

    while (Temp != NULL) {

      /* Find particles which match the numbers found above, set LeftEdge
	 and RightEdge and output particle positions to open file. */

      if (Temp->GridData->FindMatchingParticles(ParticlesFound, 
		 ParticleNumberList, LeftEdge, RightEdge, fptr) == FAIL) {
	fprintf(stderr, "Error in grid->FindMatchingParticles.\n");
	my_exit(EXIT_FAILURE);
      }

      Temp = Temp->NextGridThisLevel;

    } // end: loop over grids on this level

  } // end: loop over levels

  /* Output results. */

  printf("LeftEdge = %f %f %f\n", LeftEdge[0], LeftEdge[1], LeftEdge[2]);
  printf("RightEdge = %f %f %f\n", RightEdge[0], RightEdge[1], RightEdge[2]);

  fclose(fptr);

  my_exit(EXIT_SUCCESS);
}

void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}

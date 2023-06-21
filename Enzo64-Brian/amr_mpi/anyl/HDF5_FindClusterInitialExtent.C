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
 
#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
 
#include "hdf4.h"
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
 
// function prototypes
 
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
void QuickSortAndDrag(int List[], int left, int right,
                      int NumberToDrag1, float *DragList1[],
                      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, int   *DragList3[]);
int CommunicationInitialize(int *argc, char **argv[]);
int CommunicationFinalize();
void my_exit(int status);
 
// HDF5 function prototypes
 
#include "extern_hdf5.h"
 
 
 
 
main(int argc, char *argv[])
{
 
  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, attr_dsp_id;
  hid_t       file_type_id;
  hid_t       int_type_id;
 
  hsize_t     attr_count;
  hsize_t     part_count;
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  FILE        *log_fptr;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
 
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
  printf("NumberOfParticlesInRegion = %"ISYM"\n", ParticlesFound);
 
  /* ------------------------------------------------------------ */
  /* If there is no second hierarchy, then just output the
     particle indices. */
 
  if (argc == 3) {
 
//    int32 TempInt = ParticlesFound;
 
    if (io_log) log_fptr = fopen("FindInit_Log", "a");
 
    int_type_id = HDF5_I4;
    file_type_id = HDF5_FILE_I4;
 
    attr_count = 1;
    part_count = ParticlesFound;
 
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    file_dsp_id = H5Screate_simple(1, &part_count, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate_simple: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error );
 
    file_id = H5Fcreate("findinit.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
      assert( file_id != h5_error );
 
    dset_id = H5Dcreate(file_id, "particle_index", file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "ParticlesFound", HDF5_FILE_I4, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Acreate id: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id, int_type_id, &ParticlesFound);
      if (io_log) fprintf(log_fptr, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Dwrite(dset_id, int_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticleNumberList);
      if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Fclose(file_id);
      if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    if (io_log) fclose(log_fptr);
 
/*
    DFSDsetdims(1, &TempInt);
    DFSDsetdatastrs("particle_index", "", "", "");
    if (DFSDputdata("findinit.hdf", 1, &TempInt, (VOIDP) ParticleNumberList) ==
	HDF_FAIL) {
      fprintf(stderr, "Error writint findinit.hdf\n");
      my_exit(EXIT_FAILURE);
    }
*/
 
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
 
  QuickSortAndDrag(ParticleNumberList, 0, ParticlesFound-1, 0, NULL, 0, NULL, 0, NULL);
 
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
 
  printf("LeftEdge = %"FSYM" %"FSYM" %"FSYM"\n", LeftEdge[0], LeftEdge[1], LeftEdge[2]);
  printf("RightEdge = %"FSYM" %"FSYM" %"FSYM"\n", RightEdge[0], RightEdge[1], RightEdge[2]);
 
  fclose(fptr);
 
  my_exit(EXIT_SUCCESS);
}
 
void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}

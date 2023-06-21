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
/  EXTERNAL BOUNDARY CLASS (READ THE EXTERNAL BOUNDARY VALUES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  James Bordner,   June 2003   added USE_HDF5

/  PURPOSE:
/
************************************************************************/

#ifdef USE_HDF5

#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "hdf4.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

// This routine reads the external boundary from the provided file pointer

// HDF5 function prototypes

#include "extern_hdf5.h"

// function prototypes

int ReadListOfInts(FILE *fptr, int N, int nums[]);
int ReadListOfFloats(FILE *fptr, int N, float floats[]);




int ExternalBoundary::ReadExternalBoundaryHDF5(FILE *fptr)
{

  int Dims[MAX_DIMENSION], index, size, i;
  int BoundaryValuePresent[2*MAX_DIMENSION];
  int dim, field, TempInt, j;

  float32 *buffer;

  char hdfname[MAX_LINE_LENGTH];

  FILE *log_fptr;

  hid_t       file_id;
  hid_t       dset_id1, dset_id2;
  hid_t       float_type_id;
  hid_t       file_type_id;
  hid_t       file_dsp_id;
  hid_t       mem_dsp_id;

  hsize_t     mem_stride, mem_count, file_stride, file_count;

  hssize_t    mem_offset, file_offset;

  herr_t      h5_status;
  herr_t      h5_error = -1;

  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     mem_size, file_size;

  const char *dname_type = "BoundaryDimensionType";
  const char *dname_value = "BoundaryDimensionValue";

#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  float_type_id = (sizeof(float32)==8) ? HDF5_R8      : HDF5_R4;
  file_type_id  = (sizeof(float32)==8) ? HDF5_FILE_R8 : HDF5_FILE_R4;


  /* read general class data */

  if (fscanf(fptr, "BoundaryRank = %d\n", &BoundaryRank) != 1) {
    fprintf(stderr, "Error reading BoundaryRank.\n");
    return FAIL;
  }

  fscanf(fptr, "BoundaryDimension =");

  if (ReadListOfInts(fptr, BoundaryRank, BoundaryDimension) == FAIL) {
    fprintf(stderr, "Error reading BoundaryDimension.\n");
    return FAIL;
  }

  /* read baryon field quantities */

  if (fscanf(fptr, "NumberOfBaryonFields = %d\n", 
	     &NumberOfBaryonFields) != 1) {
    fprintf(stderr, "Error reading NumberOfBaryonFields.\n");
    return FAIL;
  }

  /* Read particle boundary type. */

  if (fscanf(fptr, "ParticleBoundaryType = %d\n",&ParticleBoundaryType) != 1) {
    fprintf(stderr, "Error reading ParticleBoundaryType.\n");
    return FAIL;
  }

  if (NumberOfBaryonFields > 0) {

    /* read field types */

    fscanf(fptr, "BoundaryFieldType = ");

    if (ReadListOfInts(fptr, NumberOfBaryonFields, BoundaryFieldType) 
        == FAIL) {
      fprintf(stderr, "Error reading BoundaryFieldType.\n");
      return FAIL;
    }

    /* read hdf file name */

    if (fscanf(fptr, "BaryonFileName = %s\n", hdfname) != 1) {
      fprintf(stderr, "Error reading BaryonFileName.\n");
      return FAIL;
    }    

    /* read BoundaryValue present line */

    fscanf(fptr, "BoundaryValuePresent = ");

    if (ReadListOfInts(fptr, BoundaryRank*2, BoundaryValuePresent) == FAIL) {
      fprintf(stderr, "Error reading BoundaryValuePresent.\n");
      return FAIL;
    }

    /* Read HDF files */

    char *logname = new char[strlen(hdfname)+6];
    strcpy(logname, hdfname);
    strcat(logname, ".log2");
    if (io_log) log_fptr = fopen(logname, "a");

    if (io_log) fprintf(log_fptr, "ReadEB start\n");
    if (io_log) fprintf(log_fptr, "  NumberOfBaryonFields %d\n", NumberOfBaryonFields);
    if (io_log) fprintf(log_fptr, "  BoundaryRank %d\n", BoundaryRank);

    for (dim = 0; dim < BoundaryRank; dim++)
    {
       if (io_log) fprintf(log_fptr, "    BoundaryDimension[%d] %d\n", dim, BoundaryDimension[dim]);
    }

    if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", hdfname);

    file_id = H5Fopen(hdfname, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fopen id: %d\n", file_id);
      assert( file_id != h5_error );

    /* loop over faces, reading each */

    for (dim = 0; dim < BoundaryRank; dim++)
      if (BoundaryDimension[dim] > 1) {

	/* calculate size and dims of flux plane */
	
	index = 0;
	size  = 1;
	Dims[0] = 1;

	for (i = 0; i < BoundaryRank; i++)
	  if (i != dim) {
	    Dims[index++] = BoundaryDimension[i];
	    size *= BoundaryDimension[i];
	  }

	index = max(BoundaryRank-1, 1);   // make index at least 1

        char *nfile = new char[2];
        char *dname1 = new char[strlen(dname_type)+3];
        char *dname2 = new char[strlen(dname_value)+3];

        nfile[0] = '\0';
        dname1[0] = '\0';
        dname2[0] = '\0';

        sprintf(nfile,"%d",dim);
        strcat(strcat(strcat(dname1,dname_type),"."),nfile);
        strcat(strcat(strcat(dname2,dname_value),"."),nfile);


        mem_size = size;
        file_size  = mem_size * 2 * NumberOfBaryonFields;

        file_dsp_id = H5Screate_simple(1, &file_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %d\n", file_dsp_id);
          assert( file_dsp_id != h5_error );

        if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", dname1);

        dset_id1 =  H5Dopen(file_id, dname1);
          if (io_log) fprintf(log_fptr, "H5Dopen id: %d\n", dset_id1);
          assert( dset_id1 != h5_error );

        if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", dname2);

        dset_id2 =  H5Dopen(file_id, dname2);
          if (io_log) fprintf(log_fptr, "H5Dopen id: %d\n", dset_id2);
          assert( dset_id2 != h5_error );

        file_offset = 0;

        mem_dsp_id = H5Screate_simple(1, &mem_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %d\n", mem_dsp_id);
          assert( mem_dsp_id != h5_error );

        file_dsp_id = H5Screate_simple(1, &file_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %d\n", file_dsp_id);
          assert( file_dsp_id != h5_error );

	/* Read HDF dims */

        /* Read attributes for BoundaryType and BoundaryValue

           NumberOfBaryonFields
           BoundaryRank
             BoundaryDimension[dim]
             Index
             Size
             OutDims[]
        */

//RH
//        if (io_log) fprintf(log_fptr, "REB hdf file %s\n", hdfname);
//        if (io_log) fprintf(log_fptr, "REB hdf rank %d\n", TempInt);
//        if (io_log) fprintf(log_fptr, "REB max buff %d\n", BoundaryRank);
//        for (i=0; i < TempInt; i++)
//        {
//          if (io_log) fprintf(log_fptr, "%d  %d\n", i, TempIntArray[i]);
//        }
//RH

	/* Check rank and dimensions (dims are stored backwards for us). */
/*
	if (TempInt != index) {
	  fprintf(stderr, "HDF file rank does not match BoundaryRank.\n");
	  return FAIL;
	}

	for (i = 0; i < index; i++)
	  if (TempIntArray[index-i-1] != Dims[i]) {
	    fprintf(stderr, "HDF file dims do not match BoundaryDims.\n");
	    fprintf(stderr, " Dims[%d] = %d   HDF Dims[%d] = %d\n", i, Dims[i],
		    index-i-1, TempIntArray[index-i-1]);
	    return FAIL;
	  }
*/
	/* Allocate temporary space. */
	
	buffer = new float32[size];

	/* loop over fields, reading each */

	for (field = 0; field < NumberOfBaryonFields; field++)
	  for (i = 0; i < 2; i++) {

          if (io_log) fprintf(log_fptr, "        dim %d : field %d : i %d\n", dim, field, i);

	    /* allocate room for BoundaryType */

	    BoundaryType[field][dim][i] = new boundary_type[size];

	    /* read BoundaryType (then convert to int) */

            mem_offset = 0;
            mem_stride = 1;
            mem_count = size;

            h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %d\n", h5_status);
              assert( h5_status != h5_error );

            file_stride = 1;
            file_count = size;

            h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &file_offset, &file_stride, &file_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Sselect file slab: %d\n", h5_status);
              assert( h5_status != h5_error );

            file_offset = file_offset + size;

            h5_status = H5Dread(dset_id1, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
              if (io_log) fprintf(log_fptr, "H5Dread boundary type: %d\n", h5_status);
              assert( h5_status != h5_error );
	    JBPERF_COUNT_READ(dset_id1, float_type_id, mem_dsp_id);

//RH
//            if (io_log) fprintf(log_fptr, "REB read BoundaryType\n");
//            if (io_log) fprintf(log_fptr, "REB getdata from %s\n", hdfname);
//            if (io_log) fprintf(log_fptr, "REB getdata rank %d\n", TempInt);
//            for (i = 0; i < TempInt; i++)
//            {
//              if (io_log) fprintf(log_fptr, "%d  %d\n", i, TempIntArray[i]);
//            }
//RH

	    for (j = 0; j < size; j++)
	      BoundaryType[field][dim][i][j] = (boundary_type) nint(buffer[j]);

	    /* read BoundaryValue */

	    if (BoundaryValuePresent[2*dim+i]) {
	      BoundaryValue[field][dim][i] = new float[size];

              h5_status = H5Dread(dset_id2, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
                if (io_log) fprintf(log_fptr, "H5Dread boundary value: %d\n", h5_status);
                assert( h5_status != h5_error );
	    JBPERF_COUNT_READ(dset_id2, float_type_id, mem_dsp_id);

//RH
//            if (io_log) fprintf(log_fptr, "REB read BoundaryValue\n");
//            if (io_log) fprintf(log_fptr, "REB getdata from %s\n", hdfname);
//            if (io_log) fprintf(log_fptr, "REB getdata rank %d\n", TempInt);
//            for (i = 0; i < TempInt; i++)
//            {
//              if (io_log) fprintf(log_fptr, "%d  %d\n", i, TempIntArray[i]);
//            }
//RH

	      for (j = 0; j < size; j++)
		BoundaryValue[field][dim][i][j] = float(buffer[j]);

	    }

	  }  // end of loop over fields

	delete buffer;

        delete nfile;
        delete dname1;
        delete dname2;

        h5_status = H5Dclose(dset_id1);
          if (io_log) fprintf(log_fptr,"H5Dclose 1: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Dclose(dset_id2);
          if (io_log) fprintf(log_fptr,"H5Dclose 2: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(mem_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

        h5_status = H5Sclose(file_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %d\n", h5_status);
          assert( h5_status != h5_error );

      }   // end of loop over dims

      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %d\n", h5_status);
        assert( h5_status != h5_error );

      if (io_log) fclose(log_fptr);

  }

  return SUCCESS;

}
#else

// HDF5 is not used, so HDF5ReadExternalBoundary should not be called

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "ExternalBoundary.h"
#include "message.h"

int ExternalBoundary::ReadExternalBoundaryHDF5(FILE *fptr)
{
  ERROR_MESSAGE;
  return FAIL;
}
#endif

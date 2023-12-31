/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (READ THE EXTERNAL BOUNDARY VALUES)
/
/  written by: Greg Bryan / Robert Harkness
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Robert Harkness, November 2005
/              Out-of-core handling for the boundary
/
/  PURPOSE:
/
************************************************************************/
 
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
 
#include "hdf4.h"
 
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
 
int WRITE_BT(boundary_type *bt_buffer, int field, int dim, int face, int slabsize, int BoundaryDimension[], int BoundaryRank, int Nfields);
int WRITE_BV(float         *bv_buffer, int field, int dim, int face, int slabsize, int BoundaryDimension[], int BoundaryRank, int Nfields);
 
 
 
int ExternalBoundary::ReadExternalBoundary(FILE *fptr)
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

#ifdef OOC_BOUNDARY
  boundary_type *bt_buffer;
  float *bv_buffer;

  int face;
  int slabsize;
#endif
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
 
  int ii = sizeof(float32);
 
  switch(ii)
  {
 
    case 4:
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
      break;
 
    case 8:
      float_type_id = HDF5_R8;
      file_type_id = HDF5_FILE_R8;
      break;
 
    default:
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
 
  }
 
  /* read general class data */
 
  if (fscanf(fptr, "BoundaryRank = %"ISYM"\n", &BoundaryRank) != 1) {
    fprintf(stderr, "Error reading BoundaryRank.\n");
    return FAIL;
  }
 
  fscanf(fptr, "BoundaryDimension =");
 
  if (ReadListOfInts(fptr, BoundaryRank, BoundaryDimension) == FAIL) {
    fprintf(stderr, "Error reading BoundaryDimension.\n");
    return FAIL;
  }
 
  /* read baryon field quantities */
 
  if (fscanf(fptr, "NumberOfBaryonFields = %"ISYM"\n",
	     &NumberOfBaryonFields) != 1) {
    fprintf(stderr, "Error reading NumberOfBaryonFields.\n");
    return FAIL;
  }
 
  /* Read particle boundary type. */
 
  if (fscanf(fptr, "ParticleBoundaryType = %"ISYM"\n",&ParticleBoundaryType) != 1) {
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
    if (io_log) fprintf(log_fptr, "  NumberOfBaryonFields %"ISYM"\n", NumberOfBaryonFields);
    if (io_log) fprintf(log_fptr, "  BoundaryRank %"ISYM"\n", BoundaryRank);
 
    for (dim = 0; dim < BoundaryRank; dim++)
    {
       if (io_log) fprintf(log_fptr, "    BoundaryDimension[%"ISYM"] %"ISYM"\n", dim, BoundaryDimension[dim]);
    }
 
    if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", hdfname);
 
    file_id = H5Fopen(hdfname, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
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
 
        sprintf(nfile,"%"ISYM,dim);
        strcat(strcat(strcat(dname1,dname_type),"."),nfile);
        strcat(strcat(strcat(dname2,dname_value),"."),nfile);
 
 
        mem_size = size;
        file_size  = mem_size * 2 * NumberOfBaryonFields;
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &file_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
          assert( file_dsp_id != h5_error );
 
        if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", dname1);
 
        dset_id1 =  H5Dopen(file_id, dname1);
          if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id1);
          assert( dset_id1 != h5_error );
 
        if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", dname2);
 
        dset_id2 =  H5Dopen(file_id, dname2);
          if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id2);
          assert( dset_id2 != h5_error );
 
        file_offset = 0;
 
        mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
          assert( mem_dsp_id != h5_error );
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &file_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
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
//        if (io_log) fprintf(log_fptr, "REB hdf rank %"ISYM"\n", TempInt);
//        if (io_log) fprintf(log_fptr, "REB max buff %"ISYM"\n", BoundaryRank);
//        for (i=0; i < TempInt; i++)
//        {
//          if (io_log) fprintf(log_fptr, "%"ISYM"  %"ISYM"\n", i, TempIntArray[i]);
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
	    fprintf(stderr, " Dims[%"ISYM"] = %"ISYM"   HDF Dims[%"ISYM"] = %"ISYM"\n", i, Dims[i],
		    index-i-1, TempIntArray[index-i-1]);
	    return FAIL;
	  }
*/
	/* Allocate temporary space. */
	
	buffer = new float32[size];

#ifdef OOC_BOUNDARY
        bt_buffer = new boundary_type[size];
        bv_buffer = new float[size];
#endif

/*
    bt_buffer = new
    bv_buffer = new

    bt_buffer[j] = buffer[j];
    bv_buffer[j] = buffer[j];
    WRITE_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
    WRITE_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
*/
 
	/* loop over fields, reading each */
 
	for (field = 0; field < NumberOfBaryonFields; field++)
	  for (i = 0; i < 2; i++) {
 
            if (io_log) fprintf(log_fptr, "        dim %"ISYM" : field %"ISYM" : i %"ISYM"\n", dim, field, i);
 
	    /* read BoundaryType (then convert to int) */
 
            mem_offset = 0;
            mem_stride = 1;
            mem_count = size;
 
            h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
              assert( h5_status != h5_error );
 
            file_stride = 1;
            file_count = size;
 
            h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &file_offset, &file_stride, &file_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
              assert( h5_status != h5_error );
 
            file_offset = file_offset + size;
 
            h5_status = H5Dread(dset_id1, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
              if (io_log) fprintf(log_fptr, "H5Dread boundary type: %"ISYM"\n", h5_status);
              assert( h5_status != h5_error );
 
#ifdef OOC_BOUNDARY

            for (j = 0; j < size; j++)
              bt_buffer[j] = (boundary_type) nint(buffer[j]);

            face = i;
            slabsize = size;

            WRITE_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields); 

#else

            BoundaryType[field][dim][i] = new boundary_type[size];

	    for (j = 0; j < size; j++)
	      BoundaryType[field][dim][i][j] = (boundary_type) nint(buffer[j]);

#endif 

	    /* read BoundaryValue */

#ifdef OOC_BOUNDARY

            if (ExternalBoundaryValueIO) {

              h5_status = H5Dread(dset_id2, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
                if (io_log) fprintf(log_fptr, "H5Dread boundary value: %"ISYM"\n", h5_status);
                assert( h5_status != h5_error );

              for (j = 0; j < size; j++)
                bv_buffer[j] = float(buffer[j]);

              face = i;
              slabsize = size;

              WRITE_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

            }

#else

            if (BoundaryValuePresent[2*dim+i]) {

              BoundaryValue[field][dim][i] = new float[size];

              h5_status = H5Dread(dset_id2, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
                if (io_log) fprintf(log_fptr, "H5Dread boundary value: %"ISYM"\n", h5_status);
                assert( h5_status != h5_error );
 
	      for (j = 0; j < size; j++)
		BoundaryValue[field][dim][i][j] = float(buffer[j]);

            }

#endif
 
	  }  // end of loop over fields

#ifdef OOC_BOUNDARY
        delete [] bt_buffer;
        delete [] bv_buffer;
#endif
 
	delete buffer;
        delete nfile;
        delete dname1;
        delete dname2;
 
        h5_status = H5Dclose(dset_id1);
          if (io_log) fprintf(log_fptr,"H5Dclose 1: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Dclose(dset_id2);
          if (io_log) fprintf(log_fptr,"H5Dclose 2: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Sclose(mem_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Sclose(file_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
      }   // end of loop over dims
 
      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      if (io_log) fclose(log_fptr);
 
  }
 
  return SUCCESS;
 
}

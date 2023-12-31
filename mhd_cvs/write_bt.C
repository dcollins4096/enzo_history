#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "OOC.h"
#include "macros_and_parameters.h"

#include "typedefs.h"
#include "global_data.h"

// function prototypes

int WRITE_BT(boundary_type *bt_buffer,
             int field, int dim, int face, int slabsize,
             int BoundaryDimension[], int BoundaryRank,
             int NumberOfBaryonFields)
{

  FILE *log;
  int io_log = 0;
  int i;

  char *Name = "BoundaryType";
  char *FieldName = "BaryonName";

  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id;
  hid_t       int_file_type_id, int_mem_type_id;

  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[4], file_count[4], file_block[4];
  hsize_t     bt_dims[4];

  hsize_t     bsize;

  hssize_t    mem_offset;
  hssize_t    file_offset[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

/*
  int BoundaryType[field][dim][face][index]
  float BoundaryValue[field][dim][face][index]
*/

    bsize = slabsize;
    int cubesize = 1;
    int facesize = 1;
    for ( i = 0; i < BoundaryRank; i++ ) {
      cubesize = cubesize * BoundaryDimension[i];
    }
    int MaxFaceSize = 1;
    for ( i = 0; i < BoundaryRank; i++ ) {
      MaxFaceSize = max(MaxFaceSize, cubesize/BoundaryDimension[i]);
    }

    facesize = cubesize/BoundaryDimension[dim];

    bt_dims[3] = MaxFaceSize;
    bt_dims[2] = 2;
    bt_dims[1] = BoundaryRank;
    bt_dims[0] = NumberOfBaryonFields; 

    file_type_id = HDF5_FILE_INT;
    mem_type_id = HDF5_INT;


    if ( MyProcessorNumber == 0 ) {

/*
  fprintf(stderr, "wbt field = %"ISYM"\n", field);
  fprintf(stderr, "wbt dim = %"ISYM"\n", dim);
  fprintf(stderr, "wbt face = %"ISYM"\n", face);
  fprintf(stderr, "wbt slabsize = %"ISYM"\n", slabsize);
  fprintf(stderr, "wbt BDims = %"ISYM" %"ISYM" %"ISYM"\n", BoundaryDimension[0], BoundaryDimension[1], BoundaryDimension[2]);
  fprintf(stderr, "wbt BRank = %"ISYM"\n", BoundaryRank);
  fprintf(stderr, "wbt NBF = %"ISYM"\n", NumberOfBaryonFields);
  fprintf(stderr, "wbt cubesize = %"ISYM"\n", cubesize);
  fprintf(stderr, "wbt facesize = %"ISYM"\n", facesize);
  fprintf(stderr, "wbt MaxFaceSize = %"ISYM"\n", MaxFaceSize);
*/

// 1D array for memory buffer
    mem_dsp_id = H5Screate_simple(1, &bsize, NULL);
      if (io_log) fprintf(stderr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
      assert( mem_dsp_id != h5_error );

// 4D array for file, or 2?
    file_dsp_id = H5Screate_simple(4, bt_dims, NULL);
      if (io_log) fprintf(stderr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error );

    if (io_log) fprintf(stderr, "Calling H5Fopen with Name = %s\n", Name);

    file_id = H5Fopen(Name, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(stderr, "H5Fopen id: %"ISYM"\n", file_id);

    if ( file_id < 0 ) {

    if (io_log) fprintf(stderr, "Calling H5Fcreate with Name = %s\n", Name);

    file_id = H5Fcreate(Name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(stderr, "H5Fcreate id: %"ISYM"\n", file_id);
      assert( file_id != h5_error );

    dset_id =  H5Dcreate(file_id, Name, file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(stderr, "H5Dcreate id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );

    } else {

    if (io_log) fprintf(stderr, "Calling H5Dopen with Name = %s\n", Name);

    dset_id =  H5Dopen(file_id, Name);
      if (io_log) fprintf(stderr, "H5Dopen id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );

    }

  mem_stride = 1;         // contiguous elements
  mem_count = facesize;   // number of elements on boundary face
  mem_offset = 0;         // zero offset in buffer
  mem_block = 1;          // single element blocks

  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(stderr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  file_stride[0] = 1;                   // contiguous elements
  file_count[0] = 1;                    // one component per call
  file_offset[0] = field;               // field component 0 to N-1
  file_block[0] = 1;                    // single element blocks

  file_stride[1] = 1;                   // contiguous elements
  file_count[1] = 1;                    // one component per call
  file_offset[1] = dim;                 // dimension 0, 1 or 2
  file_block[1] = 1;                    // single element blocks

  file_stride[2] = 1;                   // contiguous elements
  file_count[2] = 1;                    // one component per call
  file_offset[2] = face;                // face, 0 or 1
  file_block[2] = 1;                    // single element blocks

  file_stride[3] = 1;                   // contiguous elements
  file_count[3] = facesize;             // data for one plane
  file_offset[3] = 0;                   // complete field, no offset
  file_block[3] = 1;                    // single element blocks

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
    if (io_log) fprintf(stderr, "H5Sselect file slab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );


  h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, bt_buffer);
    if (io_log) fprintf(stderr, "H5Dwrite: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );


  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(stderr, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(stderr, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(stderr, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(stderr, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  }

  return SUCCESS;

}

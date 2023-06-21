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
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))


#include <hdf5.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

void pcol(float *x, int n, int m, FILE *log_fptr);


int main(int argc, char *argv[])
{

  hid_t       file_id, dset_id;
  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       mem_type_id;
  hid_t       dsp_id;
  hid_t       typ_id;

  hsize_t     size;
  hsize_t     dims[4];

  hsize_t     xdims[4];
  hsize_t     maxdims[4];

  hsize_t     ydims[3];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  hsize_t     mem_stride, mem_count;
  hsize_t     file_stride[3], file_count[3];

  hssize_t    mem_offset;
  hssize_t    file_offset[3];

  int i;
  int idim, jdim, slice;
  int ndims;
  int *ibuff;
  float *buff;

  char *file1;
  char *datax;
  char *plane;

  file1 = argv[1];
  datax = argv[2];
  plane = argv[3];

  sscanf(plane, "%d", &slice);
  printf("Slice = %d\n", slice);

  file_id = H5Fopen(file1, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert( file_id != h5_error );
  dset_id = H5Dopen(file_id, datax);
    assert( dset_id != h5_error );

  dsp_id = H5Dget_space(dset_id);
    assert( dsp_id != h5_error );
  typ_id = H5Dget_type(dset_id);
    assert( typ_id != h5_error );
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  size = 1;
  printf("Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
  {
    dims[i] = xdims[i];
    size = size * dims[i];
    printf(" Dim %d\n", (int) xdims[i]);
  }

  printf("Size %d\n", (int) size);

  jdim = xdims[1];
  idim = xdims[0];

  file_offset[0] = slice;
  file_stride[0] = 1;
  file_count[0] = 1;

  file_offset[1] = 0;
  file_stride[1] = 1;
  file_count[1] = jdim;

  file_offset[2] = 0;
  file_stride[2] = 1;
  file_count[2] = idim;

  mem_offset = 0;
  mem_stride = 1;
  mem_count = idim*jdim;

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
    assert( file_dsp_id != h5_error );
  mem_dsp_id = H5Screate_simple(1, &mem_count, NULL);
    assert( mem_dsp_id != h5_error );

  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    assert( h5_status != h5_error );
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
    assert( h5_status != h5_error );

  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
  {
    buff = new float[(int) idim*jdim];
    mem_type_id = H5T_NATIVE_FLOAT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, buff);
      printf("float read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal( typ_id, H5T_STD_I32BE ) )
  {
    ibuff = new int[(int) idim*jdim];
    mem_type_id = H5T_NATIVE_INT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, ibuff);
      printf("int read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  h5_status = H5Sclose(dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(file_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Dclose(dset_id);
    assert( h5_status != h5_error );
  h5_status = H5Fclose(file_id);
    assert( h5_status != h5_error );

  ydims[0] = jdim;
  ydims[1] = idim;

  file_dsp_id = H5Screate_simple(ndims-1, ydims, NULL);
    assert( file_dsp_id != h5_error );
  file_id = H5Fcreate( "Slice", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert( file_id != h5_error );
  dset_id = H5Dcreate( file_id, "Slice", mem_type_id, file_dsp_id, H5P_DEFAULT);
    assert( dset_id != h5_error );

  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
  {
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, buff);
      printf("float write status %d\n", (int) h5_status); 
      assert( h5_status != h5_error );
  }

  if ( H5Tequal( typ_id, H5T_STD_I32BE ) )
  {
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, ibuff);
      printf("int write status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  h5_status = H5Tclose(typ_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(mem_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(file_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Dclose(dset_id);
    assert( h5_status != h5_error );
  h5_status = H5Fclose(file_id);
    assert( h5_status != h5_error );

}



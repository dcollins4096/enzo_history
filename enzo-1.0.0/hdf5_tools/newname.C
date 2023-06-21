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
#include <math.h>

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

  herr_t      h5_status;
  herr_t      h5_error = -1;


  int i;
  int ndims;
  int *ibuff;
  float *rbuff;
  double *dbuff;

  char *filename;
  char *newfile;
  char *dset_in;
  char *dset_out;

  filename = argv[1];
  newfile = argv[2];
  dset_in = argv[3];
  dset_out = argv[4];

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert( file_id != h5_error);
  dset_id = H5Dopen(file_id, dset_in);
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

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
    assert( file_dsp_id != h5_error );
  mem_dsp_id = H5Screate_simple(1, &size, NULL);
    assert( mem_dsp_id != h5_error );

  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
  {
    rbuff = new float[(int) size];
    mem_type_id = H5T_NATIVE_FLOAT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, rbuff);
      printf("float read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )
  {
    dbuff = new double[(int) size];
    mem_type_id = H5T_NATIVE_DOUBLE;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, dbuff );
      printf("double read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal( typ_id, H5T_STD_I32BE ) )
  {
    ibuff = new int[(int) size];
    mem_type_id = H5T_NATIVE_FLOAT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, ibuff);
      printf("int read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  h5_status = H5Dclose(dset_id);
  h5_status = H5Fclose(file_id);

  file_id = H5Fcreate(newfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert( file_id != h5_error );
  dset_id = H5Dcreate(file_id, dset_out, mem_type_id, file_dsp_id, H5P_DEFAULT);
    assert( dset_id != h5_error );

  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
  {
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, rbuff);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )
  {
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, dbuff);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal( typ_id, H5T_STD_I32BE ) )
  {
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, ibuff);
      assert( h5_status != h5_error );
  }

  h5_status = H5Sclose(file_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(mem_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Tclose(typ_id);
    assert( h5_status != h5_error );
  h5_status = H5Dclose(dset_id);
    assert( h5_status != h5_error );
  h5_status = H5Fclose(file_id);
    assert( h5_status != h5_error );
}

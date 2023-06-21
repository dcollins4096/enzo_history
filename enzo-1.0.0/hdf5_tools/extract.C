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

typedef void *VOIDP;

int main(int argc, char *argv[])
{

  hid_t       file_id, dset_id;
  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       mem_type_id, file_type_id;
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
  int dtype;

  char *file1;
  char *datax;
  char *file2;

  file1 = argv[1];
  datax = argv[2];
  file2 = argv[3];

  dtype = 0;

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

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
    assert( file_dsp_id != h5_error );
  mem_dsp_id = H5Screate_simple(1, &size, NULL);
    assert( mem_dsp_id != h5_error );

  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
  {
    file_type_id = H5T_IEEE_F32BE;
    mem_type_id = H5T_NATIVE_FLOAT;
    dtype = 1;
  }

  if ( H5Tequal( typ_id, H5T_IEEE_F64BE ) )
  {
    file_type_id = H5T_IEEE_F64BE;
    mem_type_id = H5T_NATIVE_DOUBLE;
    dtype = 2;
  }

  if ( H5Tequal( typ_id, H5T_STD_I32BE ) )
  {
    file_type_id = H5T_STD_I32BE;
    mem_type_id = H5T_NATIVE_INT;
    dtype = 1;
  }

  int *buff = new int[size*dtype];

  h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) buff);
    printf("read status %d\n", (int) h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Dclose(dset_id);
    assert( h5_status != h5_error );
  h5_status = H5Fclose(file_id);
    assert( h5_status != h5_error );

  file_id = H5Fcreate(file2, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert( file_id != h5_error );
  dset_id = H5Dcreate(file_id, datax, file_type_id, file_dsp_id, H5P_DEFAULT);
    assert( dset_id != h5_error );
  h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) buff);
    printf("write status %d\n", (int) h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Dclose(dset_id);
    assert( h5_status != h5_error );
  h5_status = H5Fclose(file_id);
    assert( h5_status != h5_error );

  h5_status = H5Sclose(dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Tclose(typ_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(mem_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(file_dsp_id);
    assert( h5_status != h5_error );

}

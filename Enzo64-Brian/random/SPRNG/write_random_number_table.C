#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "macros_and_parameters.h"

#if defined(IRIS4) || defined(SUN)   || defined(COMPAQ) || \
    defined(IA64) || defined(LINUX) || defined(NEC) || defined(CRAYX1)
extern "C" void write_random_number_table_(int *nx, int *slice, double tab[])
#endif
#if defined(SP2) || defined(HP)
extern "C" void write_random_number_table(int *nx, int *slice, double tab[])
#endif
{

  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id;

  hsize_t     slab_dims[2];
  hsize_t     slab_rank;
  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[2], file_count[2], file_block[2];
  hsize_t     attr_count;
  hsize_t     bufsize;

  hssize_t    mem_offset;
  hssize_t    file_offset[2];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int i;
  int Nblock, Nslice;

  Nslice = *nx/2;
  Nblock = 2 * 5 * ( 2*(*nx/2) + 1 ) * ( 2*(*nx/2) + 1 );

  slab_rank = 2;
  slab_dims[0] = Nslice;
  slab_dims[1] = Nblock;
  bufsize = Nblock;

  mem_type_id = HDF5_R8;
  file_type_id = HDF5_FILE_R8;

  mem_dsp_id = H5Screate_simple(1, &bufsize, NULL);
    assert( mem_dsp_id != h5_error );

  file_dsp_id = H5Screate_simple(((Eint32) slab_rank), slab_dims, NULL);
    assert( file_dsp_id != h5_error );

  if ( *slice == 0 ) {

    file_id = H5Fcreate("RandomField", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      assert( file_id != h5_error );

    dset_id =  H5Dcreate(file_id, "RandomField", file_type_id, file_dsp_id, H5P_DEFAULT);
      assert( dset_id != h5_error );

  } else {

    file_id = H5Fopen("RandomField", H5F_ACC_RDWR, H5P_DEFAULT);
      assert( file_id != h5_error );

    dset_id =  H5Dopen(file_id, "RandomField");
      assert( dset_id != h5_error );
  }


// Data in memory is considered 1D, stride 1, with zero offset

  mem_stride = 1;      // contiguous elements
  mem_count = Nblock;  // number of elements in chunk
  mem_offset = 0;      // zero offset in buffer
  mem_block = 1;       // single element blocks

  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    assert( h5_status != h5_error );


  file_stride[0] = 1;      // contiguous elements
  file_count[0] = 1;       // one component per call
  file_offset[0] = *slice; // component Part of Npart
  file_block[0] = 1;       // single element blocks

  file_stride[1] = 1;      // contiguous elements
  file_count[1] = Nblock;  // size of chunk
  file_offset[1] = 0;      // complete field, no offset
  file_block[1] = 1;       // single element blocks

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
    assert( h5_status != h5_error );


  h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, tab);
    assert( h5_status != h5_error );

  h5_status = H5Dclose(dset_id);
    assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
    assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
    assert( h5_status != h5_error );

  h5_status = H5Fclose(file_id);
    assert( h5_status != h5_error );

}

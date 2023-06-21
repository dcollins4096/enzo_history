/***********************************************************************
/
/  OUTPUT THE FIELD TO A HDF5 FILE
/
/  written by: Robert Harkness
/  date:       April 2002
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "macros_and_parameters.h"

// function prototypes

void fcol(float *x, int n, int m, FILE *log_fptr);

// HDF5 function prototypes

#include "extern_hdf5.h"




int WriteFieldHDF5(int Rank, int Dims[3], float *Field, char *Name, int Part, int Npart)
{


  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id, attr_type_id;

  hsize_t     OutDims[3];
  hsize_t     Dimm;
  hsize_t     Slab_Dims[4];
  hsize_t     Slab_Rank;
  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[4], file_count[4], file_block[4];
  hsize_t     attr_count;

  hssize_t    mem_offset;
  hssize_t    file_offset[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int         dim;

  int         component_rank_attr;
  int         component_size_attr;
  int         field_rank_attr;
  int         field_dims_attr[3];

  FILE        *dumpfile;

#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  if (io_log) dumpfile = fopen("DumpWF","a");

  if (io_log) printf("On entry to WriteFieldHDF5\n");
  if (io_log) printf("  Rank %d\n",Rank);
  if (io_log) printf("  Dims %d  %d  %d\n",Dims[0],Dims[1],Dims[2]);
  if (io_log) printf("  Name %s\n",Name);
  if (io_log) printf("  Part %d of %d\n",Part,Npart);

  /* GB: Reverse dim ordering since we are using fortran array ordering
         (actually, we don't have to do this here). */

  for ( dim = 0; dim < Rank; dim++ )
  {
    OutDims[dim] = Dims[dim];
  }

  /* OutDims[Rank-dim-1] = Dims[dim]; */


  Slab_Rank = Rank+1;

  Slab_Dims[0] = Npart;

  for ( dim = 1; dim < Slab_Rank; dim++ )
  {
    Slab_Dims[dim] = OutDims[dim-1];
  }

  if (io_log) printf("  Extended Rank %d\n", (int) Slab_Rank);

  for ( dim = 0; dim < Slab_Rank; dim++ )
  {
    if (io_log) printf("    %d:  %d\n",dim, (int) Slab_Dims[dim]);
  }

  Dimm = 1;

  for ( dim = 0; dim < Rank; dim++ )
  {
    Dimm = Dimm * Dims[dim];
  }

  if (io_log) printf("  Grid Elements %d\n", (int) Dimm);

  if (io_log) fcol(Field, (int) Dimm, 8, dumpfile);

  component_rank_attr = Npart;
  component_size_attr = Dimm;

  field_rank_attr = Rank;

  for ( dim = 0; dim < Rank; dim++ )
  {
    field_dims_attr[dim] = Dims[dim];
  }

  /* HDF5
   
     The HDF4 DFSD interface does not have a direct analogue in HDF5.
     In particular, there is no "append" on open and there are no
     specific features to store rank, dimensions, scales etc.
     These are replaced by dataset attributes, which can be named.
     To eliminate redundancy and to keep the data structure close
     to the original, each ENZO file will contain a single dataset
     of the same name (as opposed to separate datasets for each 
     component).

     The dataspace uses the number of components for each of the
     the dimensions of the dataset as an additional "dimension".
     For example, a 3D scalar field of 4x4x4 points will have a 
     dataspace of {1,4,4,4}, while a 3D vector field of, say,
     {Vx,Vy,Vz} and 4x4x4 points will have a dataspace of {3,4,4,4}.
     
     Slab I/O is used in preparation for parallel I/O.

     It is ASSUMED that this routine is called with
     Part = {0,1,2,...,Npart-1} to write Npart components of equal length
     (the product of the field dimensions {Dims[i] for i=0,Rank-1}.
     Each component is offset in the file by {0,1,2,...,Npart-1} * the field size.

     The rank and dimensions of the field and the number of field
     components are HDF5 attributes.

  */

  int ii = sizeof(float);

  switch(ii)
  {

    case 4:
      mem_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
      break;

    case 8:
      mem_type_id = HDF5_R8;
      file_type_id = HDF5_FILE_R8;
      break;

    default:
      mem_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;

  }

// Data in memory is considered 1D, stride 1, with zero offset

  mem_stride = 1;      // contiguous elements
  mem_count = Dimm;    // number of elements in field
  mem_offset = 0;      // zero offset in buffer
  mem_block = 1;       // single element blocks

// 1D memory model

  mem_dsp_id = H5Screate_simple(1, &Dimm, NULL);

  if (io_log) printf("H5Screate mem_dsp_id %d\n",mem_dsp_id);
  assert( mem_dsp_id != h5_error );

  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);

  if (io_log) printf("H5Sselect mem slab: %d\n",h5_status);
  assert( h5_status != h5_error );


// Data in the file is (1+Rank)D with Npart components per grid point.
// Offset[0] is the component Part of Npart components.  Data for each
// Part are contiguous in the file, so stride = 1.

  file_stride[0] = 1;      // contiguous elements
  file_count[0] = 1;       // one component per call
  file_offset[0] = Part;   // component Part of Npart
  file_block[0] = 1;       // single element blocks

  for ( dim = 1; dim < Slab_Rank; dim++ )
  {
    file_stride[dim] = 1;                   // contiguous elements
    file_count[dim] = OutDims[dim-1];       // field dimensions
    file_offset[dim] = 0;                   // complete field, no offset
    file_block[dim] = 1;                    // single element blocks
  }

  file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);

  if (io_log) printf("H5Screate file_dsp_id %d\n",file_dsp_id);
  assert( file_dsp_id != h5_error );

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);

  if (io_log) printf("H5Sselect file slab: %d\n",h5_status);
  assert( h5_status != h5_error );

  /* If Part is zero, create an HDF5 file with a single dataset of 
     the same name and attach the dataset attributes, otherwise
     open an existing file and dataset. */

  if ( Part == 0 )
  {
    if (io_log) printf("Calling H5Fcreate with Name = %s\n",Name);

    file_id = H5Fcreate(Name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (io_log) printf("H5Fcreate File id: %d\n",file_id);
    assert( file_id != h5_error );

    if (io_log) printf("Calling H5Dcreate with Name = %s\n",Name);

    dset_id =  H5Dcreate(file_id, Name, file_type_id, file_dsp_id, H5P_DEFAULT);

    if (io_log) printf("H5Dcreate Dataset id: %d\n",dset_id);
    assert( dset_id != h5_error );

//

    attr_count = 1;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);

    assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Component_Rank",  HDF5_FILE_I4, attr_dsp_id, H5P_DEFAULT);

    assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  HDF5_I4, &component_rank_attr);

    assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);

    assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);

    assert( h5_status != h5_error );

//

    attr_count = 1;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);

    assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Component_Size",  HDF5_FILE_I4, attr_dsp_id, H5P_DEFAULT);

    assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  HDF5_I4, &component_size_attr);

    assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);

    assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);

    assert( h5_status != h5_error );

//

    attr_count = 1;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);

    assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Rank", HDF5_FILE_I4, attr_dsp_id, H5P_DEFAULT);

    assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  HDF5_I4, &field_rank_attr);

    assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);

    assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);

    assert( h5_status != h5_error );

//

    attr_count = Rank;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);

    assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Dimensions", HDF5_FILE_I4, attr_dsp_id, H5P_DEFAULT);

    assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  HDF5_I4, field_dims_attr);

    assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);

    assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);

    assert( h5_status != h5_error );

  }

  else

  {
    if (io_log) printf("Calling H5Fopen with Name = %s\n",Name);

    file_id = H5Fopen(Name, H5F_ACC_RDWR, H5P_DEFAULT);

    if (io_log) printf("H5Fopen File id: %d\n",file_id);
    assert( file_id != h5_error );

    if (io_log) printf("Calling H5Dopen with Name = %s\n",Name);

    dset_id =  H5Dopen(file_id, Name);

    if (io_log) printf("H5Dopen Dataset id: %d\n",dset_id);
    assert( dset_id != h5_error );

  }


  h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, Field);

  if (io_log) printf("H5Dwrite %d\n",h5_status);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(dset_id);

  if (io_log) printf("H5Dclose %d\n",h5_status);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);

  if (io_log) printf("H5Sclose %d\n",h5_status);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);

  if (io_log) printf("H5Sclose %d\n",h5_status);
  assert( h5_status != h5_error );

  h5_status = H5Fclose(file_id);

  if (io_log) printf("H5Fclose %d\n",h5_status);
  assert( h5_status != h5_error );

  if (io_log) printf("Exit WriteFieldHDF5\n");

  if (io_log) fclose(dumpfile);

  return SUCCESS;
}

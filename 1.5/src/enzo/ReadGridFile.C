/***********************************************************************
/
/  READ AN HDF5 FILE
/
/  written by: Robert Harkness
/  date:       February 2007
/              April 2008
/              May 2008
/
/  PURPOSE: Read pre-decomposed grid tiles
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

 


 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
void my_exit(int status);
 
// function prototypes
 
void fcol(float *x, int n, int m, FILE *log_fptr);
 
// HDF5 function prototypes
 

 
 
 
 
int ReadGridFile(char *name, int Rank, int Dim[], int StartIndex[],
                      int EndIndex[], int BufferOffset[], float *buffer,
                      inits_type **tempbuffer, int Part, int Npart)
{
 
  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id, attr_type_id;
 
  hsize_t     xfer_size;
  hsize_t     Slab_Dims[4];
  int         Slab_Rank;
  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[4], file_count[4], file_block[4];
  hsize_t     slab_stride[4], slab_count[4], slab_block[4];
  hsize_t     attr_count;
 
  hssize_t    mem_offset;
  hssize_t    file_offset[4];
  hssize_t    slab_offset[4];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int component_rank_attr;
  int component_size_attr;
  int field_rank_attr;
  int field_dims_attr[3];
 
  int dim, i, j, k;
 
  int TempInt;
  int TempIntArray[MAX_DIMENSION];
 
  FILE *log_fptr;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
  char pid[MAX_TASK_TAG_SIZE];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  char *logname = new char[MAX_NAME_LENGTH];
  strcpy(logname, "GRlog.");
  strcat(logname, pid);
 
  if (io_log) log_fptr = fopen(logname, "a");
 
  delete logname;
 
  if (io_log) fprintf(log_fptr, "\n");
  if (io_log) fprintf(log_fptr, "GR file %s\n", name);
  if (io_log) fprintf(log_fptr, "GR rank %"ISYM"\n", Rank);
 
  for ( i = 0; i < Rank; i++)
  {
    if (io_log) fprintf(log_fptr,"%"ISYM"  %"ISYM"  %"ISYM"  %"ISYM"\n", Dim[i], StartIndex[i], EndIndex[i], BufferOffset[i]);
  }
 
// This routine reads only data from Inits: inits_type is 32- or 64-bit
 
  int ii = sizeof(inits_type);
 
  if (io_log) fprintf(log_fptr, "GR size of inits_type is %"ISYM"\n", ii);
 
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
 
  // Open the HDF5 file and dataset
 
  if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", name);
 
  file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
  fprintf(stderr, "GR H5Fopen %s on CPU %"ISYM"\n", name, MyProcessorNumber);
    if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
    if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", name);
 
  dset_id =  H5Dopen(file_id, name);
    if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
    if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}

 
  // Read the HDF5 attributes 
 
 
  if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");
 
  attr_id = H5Aopen_name(dset_id, "Component_Rank");
    if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aread(attr_id, HDF5_INT, &component_rank_attr);
    if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  if (io_log) fprintf(log_fptr, "COMPONENT_RANK %"ISYM"\n", component_rank_attr);
 
 
 
 
  if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");
 
  attr_id = H5Aopen_name(dset_id, "Component_Size");
    if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aread(attr_id, HDF5_INT, &component_size_attr);
    if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  if (io_log) fprintf(log_fptr, "COMPONENT_SIZE %"ISYM"\n", component_size_attr);
 
 
 
 
  if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");
 
  attr_id = H5Aopen_name(dset_id, "Rank");
    if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aread(attr_id, HDF5_INT, &field_rank_attr);
    if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  if (io_log) fprintf(log_fptr, "RANK %"ISYM"\n", field_rank_attr);
 
 
 
 
  if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");
 
  attr_id = H5Aopen_name(dset_id, "Dimensions");
    if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  attr_count = field_rank_attr;
 
  attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
    if (io_log) fprintf(log_fptr, "H5Screate_simple %"ISYM"\n", attr_dsp_id);
    if( attr_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aread(attr_id, HDF5_INT, field_dims_attr);
    if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  for (dim = 0; dim < field_rank_attr; dim++)
  {
    if (io_log) fprintf(log_fptr, "DIMS %"ISYM":  %"ISYM"\n", dim, (int) field_dims_attr[dim]);
  }
 
 
 
 
  xfer_size = 1;
 
  for ( dim = 0; dim < field_rank_attr; dim++ )
  {
    xfer_size = xfer_size * field_dims_attr[dim];
  }
 
  if (io_log) fprintf(log_fptr, "  Grid Elements %"ISYM"\n", (int) xfer_size);
 
  // Size of ENTIRE array (should be the same)
 
  TempInt = field_rank_attr;
 
  for (dim = 0; dim < Rank; dim++)
  {
    TempIntArray[dim] = field_dims_attr[dim];
  }
 
  Slab_Rank = field_rank_attr+1;
 
  Slab_Dims[0] = component_rank_attr;
 
  for ( dim = 1; dim < Slab_Rank; dim++ )
  {
    Slab_Dims[dim] = field_dims_attr[dim-1];
  }
 
  if (io_log) fprintf(log_fptr, "  Extended Rank %"ISYM"\n", (int) Slab_Rank);
  for ( dim = 0; dim < Slab_Rank; dim++ )
  {
    if (io_log) fprintf(log_fptr, "    %"ISYM":  %"ISYM"\n", dim, (int) Slab_Dims[dim]);
  }
 
  // Error check
 
  if (Rank < 1 || Rank > 3) {
    fprintf(stderr, "Rank %"ISYM" not supported.\n", Rank);
    return FAIL;
  }
 
  if (Npart != component_rank_attr) {
    fprintf(stderr, "Npart and Component_Rank do not agree!\n");
    return FAIL;
  }
 
  if (TempInt != Rank) {
    fprintf(stderr, "Rank mismatch in %s.\n", name);
    return FAIL;
  }
 
  // Compute size of the expected HDF5 field
 
  int size = 1;
 
  for (dim = 0; dim < Rank; dim++)
    size *= (EndIndex[dim]-StartIndex[dim]+1);

  // Check that size = xfer_size

  if (size != xfer_size ) {
    fprintf(stderr, "Expected size and HDF5 grid file size do not agree.\n");
    return FAIL;
  }

  for (dim = 0; dim < Rank; dim++) {
    if ( field_dims_attr[dim] != (EndIndex[dim]-StartIndex[dim]+1) )
      fprintf(stderr, "Size mismatch: %s : %"ISYM"  %"ISYM"\n", name, field_dims_attr[dim], (EndIndex[dim]-StartIndex[dim]+1) );
  } 

  // Allocate space for temp buffer
 
  if (io_log) fprintf(log_fptr, "Allocate %"ISYM" inits_types for *tempbuffer\n", size);
 
  (*tempbuffer) = new inits_type[size];

  // Read the whole data set

  // Data in memory is considered 1D, stride 1, with zero offset
 
  mem_stride = 1;           // contiguous elements
  mem_count = xfer_size;    // number of elements in field
  mem_offset = 0;           // zero offset in buffer
  mem_block = 1;            // single element blocks
 
  // 1D memory model
 
  mem_dsp_id = H5Screate_simple((Eint32) 1, &xfer_size, NULL);
    if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
    if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
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
    file_count[dim] = TempIntArray[dim-1];  // field dimensions
    file_offset[dim] = 0;                   // complete field, no offset
    file_block[dim] = 1;                    // single element blocks
  }
 
  file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
    if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
    if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
    if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) (*tempbuffer));
    if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

 
  if (io_log) fprintf(log_fptr, "Dump *tempbuffer xfer_size = %"ISYM"\n", (int) xfer_size);
 
  // If buffer is not defined, then just return w/o clearing (*tempbuffer)
 
  if (buffer == NULL)
  {
    if (io_log) fprintf(log_fptr, "NULL buffer, return data in *tempbuffer\n");
    if (io_log) fclose(log_fptr);
    return SUCCESS;
  }
 
  // Clear buffer (primarily to prevent errors in unused area)
 
  if (io_log) fprintf(log_fptr, "Non-NULL buffer, copy *tempbuffer to buffer\n");
 
  size = 1;
 
  for (dim = 0; dim < Rank; dim++)
    size *= Dim[dim];
 
  if (io_log) fprintf(log_fptr, "Buffer = 0, size = %"ISYM"\n", size);
 
  for (i = 0; i < size; i++)
    buffer[i] = 0;
 
  // Copy field into real array

  if (Rank == 1)
    for (i = StartIndex[0]; i <= EndIndex[0]; i++)
      buffer[i] =
	(float) ((*tempbuffer)[i-StartIndex[0]]);
 
  if (Rank == 2)
    for (j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	buffer[j*Dim[0] + i] =
	  (float) ((*tempbuffer)[(j-StartIndex[1])*TempIntArray[0] +
		                 (i-StartIndex[0])]);

  if (Rank == 3)
    for (k = StartIndex[2]; k <= EndIndex[2]; k++)
      for (j = StartIndex[1]; j <= EndIndex[1]; j++)
	for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	  buffer[k*Dim[0]*Dim[1] + j*Dim[0] + i] =
	    (float) ((*tempbuffer)[(k-StartIndex[2])*TempIntArray[0]*TempIntArray[1] +
	                           (j-StartIndex[1])*TempIntArray[0] +
	                           (i-StartIndex[0])]);
 
  if (Rank > 0)
    if (io_log) fprintf(log_fptr, "Dim[0] = %"ISYM", TempIntArray[0] = %"ISYM", StartIndex[0] = %"ISYM", EndIndex[0] = %"ISYM"\n", Dim[0], TempIntArray[0], StartIndex[0], EndIndex[0]);
  if (Rank > 1)
    if (io_log) fprintf(log_fptr, "Dim[1] = %"ISYM", TempIntArray[1] = %"ISYM", StartIndex[1] = %"ISYM", EndIndex[1] = %"ISYM"\n", Dim[1], TempIntArray[1], StartIndex[1], EndIndex[1]);
  if (Rank > 2)
    if (io_log) fprintf(log_fptr, "Dim[2] = %"ISYM", TempIntArray[2] = %"ISYM", StartIndex[2] = %"ISYM", EndIndex[2] = %"ISYM"\n", Dim[2], TempIntArray[2], StartIndex[2], EndIndex[2]);

  // clean up
 
  if (io_log) fprintf(log_fptr, "De-Allocate (*tempbuffer)\n");
 
  delete [] (*tempbuffer);
 
  if (io_log) fclose(log_fptr);
 
  return SUCCESS;
 
}

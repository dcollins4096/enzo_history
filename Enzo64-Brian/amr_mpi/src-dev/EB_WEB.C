/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (WRITE THE EXTERNAL BOUNDARY VALUES)
/
/  written by: Greg Bryan / Robert Harkness
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       February, 2004
/              Direct or indirect SRB driver
/  modified2:  Robert Harkness
/  date:       November, 2005
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
 
// This routine writes the external boundary to the provided file pointer
 
// HDF5 function prototypes
 
#include "extern_hdf5.h"
 
// function prototypes
 
void WriteListOfInts(FILE *fptr, int N, int nums[]);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
 
int FileMover(char *file_to_move);
 
int READ_BT(boundary_type *bt_buffer, int field, int dim, int face, int slabsize, int BoundaryDimension[], int BoundaryRank, int Nfields);
int READ_BV(float         *bv_buffer, int field, int dim, int face, int slabsize, int BoundaryDimension[], int BoundaryRank, int Nfields);
 
 
 
int ExternalBoundary::WriteExternalBoundary(FILE *fptr, char *hdfname)
{
 
  int dim, field, i, j, index, ret, size;
  int BoundaryValuePresent[MAX_DIMENSION*2], Temp[MAX_DIMENSION];
  int file_status;
  float32 *buffer;
 
  FILE *log_fptr;
 
  hid_t       file_id;
  hid_t       dset_id1, dset_id2;
  hid_t       float_type_id;
  hid_t       file_type_id;
  hid_t       attr_id;
  hid_t       attr_dsp_id;
  hid_t       file_dsp_id;
  hid_t       mem_dsp_id;
 
  hsize_t     mem_stride, mem_count, file_stride, file_count;
 
  hssize_t    mem_offset, file_offset;
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     mem_size, file_size;
  hsize_t     n_attr;
 
  int         Dims[MAX_DIMENSION];
 
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
 
  // Save general class data
 
  fprintf(fptr, "BoundaryRank         = %"ISYM"\n", BoundaryRank);
  fprintf(fptr, "BoundaryDimension    = ");
 
  WriteListOfInts(fptr, BoundaryRank, BoundaryDimension);
 
  // Save baryon field quantities
 
  fprintf(fptr, "NumberOfBaryonFields = %"ISYM"\n", NumberOfBaryonFields);
 
  // Save particle boundary type
 
  fprintf(fptr, "ParticleBoundaryType = %"ISYM"\n", ParticleBoundaryType);
 
  if (NumberOfBaryonFields > 0) {
 
    fprintf(fptr, "BoundaryFieldType    = ");
 
    WriteListOfInts(fptr, NumberOfBaryonFields, BoundaryFieldType);
 
    fprintf(fptr, "BaryonFileName       = %s\n", hdfname);
 
    // Write out information about the BoundaryValue fields
 
    for (dim = 0; dim < BoundaryRank; dim++)
      for (i = 0; i < 2; i++)
	if (BoundaryValue[0][dim][i] == NULL)
	  BoundaryValuePresent[2*dim+i] = FALSE;
	else
	  BoundaryValuePresent[2*dim+i] = TRUE;
 
    fprintf(fptr, "BoundaryValuePresent = ");
 
    WriteListOfInts(fptr, BoundaryRank*2, BoundaryValuePresent);
 
    char *logname = new char[strlen(hdfname)+5];
    strcpy(logname, hdfname);
    strcat(logname, ".log");
    if (io_log) log_fptr = fopen(logname, "a");
    delete logname;
 
    if (io_log) fprintf(log_fptr, "WriteEB start\n");
    if (io_log) fprintf(log_fptr, "  NumberOfBaryonFields %"ISYM"\n", NumberOfBaryonFields);
    if (io_log) fprintf(log_fptr, "  BoundaryRank %"ISYM"\n", BoundaryRank);
 
    for (dim = 0; dim < BoundaryRank; dim++)
    {
       if (io_log) fprintf(log_fptr, "    BoundaryDimension[%"ISYM"] %"ISYM"\n", dim, BoundaryDimension[dim]);
    }
 
    if (io_log) fprintf(log_fptr, "H5Fcreate with Name = %s\n", hdfname);
 
    file_id = H5Fcreate(hdfname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
      assert( file_id != h5_error );
 
    for (dim = 0; dim < BoundaryRank; dim++)
      if (BoundaryDimension[dim] > 1) {
 
	// Calculate size and dims of flux plane
	
	index   = 0;
	size    = 1;
	Temp[0] = 1;
 
	for (i = 0; i < BoundaryRank; i++)
	  if (i != dim) {
	    Temp[index++] = BoundaryDimension[i];
	    size *= BoundaryDimension[i];
	  }
 
	index = max(BoundaryRank-1, 1);   // make index at least 1
 
	// Reverse outdims (for HDF)
 
	for (i = 0; i < index; i++)
	  OutDims[index-i-1] = Temp[i];
 
        for (i = 0; i < index; i++)
        {
          Dims[i] = BoundaryDimension[i];
        }
 
        if (io_log) fprintf(log_fptr, "    Index %"ISYM"\n", index);
        for (i = 0; i < index; i++)
        {
          if (io_log) fprintf(log_fptr, "      OutDims[%"ISYM"] %"ISYM"\n", i, (int) OutDims[i]);
        }
        if (io_log) fprintf(log_fptr, "    Size %"ISYM"\n", size);
 
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
 
        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dname1);
 
        dset_id1 =  H5Dcreate(file_id, dname1, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id1);
          assert( dset_id1 != h5_error );
 
        if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dname2);
 
        dset_id2 =  H5Dcreate(file_id, dname2, file_type_id, file_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id2);
          assert( dset_id2 != h5_error );
 
        file_offset = 0;
 
        mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
          assert( mem_dsp_id != h5_error );
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &file_size, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
          assert( file_dsp_id != h5_error );
 
        // Attributes for BoundaryType and BoundaryValue
 
        n_attr = 1;
 
        attr_dsp_id = H5Screate_simple((Eint32) 1, &n_attr, NULL);
         if (io_log) fprintf(log_fptr, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
         assert( attr_dsp_id != h5_error );
 
        attr_id = H5Acreate(dset_id1, "NumberOfBaryonFields", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Acreate: %"ISYM"\n", attr_id);
          assert( attr_id != h5_error );
 
        h5_status = H5Awrite(attr_id,  HDF5_INT, &NumberOfBaryonFields);
          if (io_log) fprintf(log_fptr, "H5Awrite: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        attr_id = H5Acreate(dset_id1, "BoundaryRank", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Acreate: %"ISYM"\n", attr_id);
          assert( attr_id != h5_error );
 
        h5_status = H5Awrite(attr_id,  HDF5_INT, &BoundaryRank);
          if (io_log) fprintf(log_fptr, "H5Awrite: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        attr_id = H5Acreate(dset_id1, "Index", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Acreate: %"ISYM"\n", attr_id);
          assert( attr_id != h5_error );
 
        h5_status = H5Awrite(attr_id,  HDF5_INT, &index);
          if (io_log) fprintf(log_fptr, "H5Awrite: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        attr_id = H5Acreate(dset_id1, "Size", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Acreate: %"ISYM"\n", attr_id);
          assert( attr_id != h5_error );
 
        h5_status = H5Awrite(attr_id,  HDF5_INT, &size);
          if (io_log) fprintf(log_fptr, "H5Awrite: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Sclose(attr_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        n_attr = BoundaryRank;
 
        attr_dsp_id = H5Screate_simple((Eint32) 1, &n_attr, NULL);
          if (io_log) fprintf(log_fptr, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
          assert( attr_dsp_id != h5_error );
 
        attr_id = H5Acreate(dset_id1, "BoundaryDimension", HDF5_FILE_INT, attr_dsp_id, H5P_DEFAULT);
          if (io_log) fprintf(log_fptr, "H5Acreate: %"ISYM"\n", attr_id);
          assert( attr_id != h5_error );
 
        h5_status = H5Awrite(attr_id,  HDF5_INT, Dims);
          if (io_log) fprintf(log_fptr, "H5Awrite: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Aclose(attr_id);
          if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Sclose(attr_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
	// Allocate temporary space

/*
    bt_buffer = new
    bv_buffer = new
    READ_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
    READ_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

    buffer[j] = float32(bt_buffer[j]);
    buffer[j] = float32(bv_buffer[j]);
*/

	buffer = new float32[size];

#ifdef OOC_BOUNDARY
        bt_buffer = new boundary_type[size];
        bv_buffer = new float[size];
#endif
 
	for (field = 0; field < NumberOfBaryonFields; field++)
	  for (i = 0; i < 2; i++) {
 
          if (io_log) fprintf(log_fptr, "        dim %"ISYM" : field %"ISYM" : i %"ISYM"\n", dim, field, i);
 
	    // Write out BoundaryType (convert to float first)

#ifdef OOC_BOUNDARY

            face = i;
            slabsize = size;

            if (ExternalBoundaryTypeIO)
            READ_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

            for (j = 0; j < size; j++)
              buffer[j] = float32(bt_buffer[j]);

#else
 
	    for (j = 0; j < size; j++)
	      buffer[j] = float32(BoundaryType[field][dim][i][j]);
#endif
 
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
 
            h5_status = H5Dwrite(dset_id1, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
              if (io_log) fprintf(log_fptr, "H5Dwrite boundary type: %"ISYM"\n", h5_status);
              assert( h5_status != h5_error );
 
	    /* write out BoundaryValue */

#ifdef OOC_BOUNDARY

            if (ExternalBoundaryValueIO) {

              face = i;
              slabsize = size;

              READ_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields); 

              for (j = 0; j < size; j++)
                buffer[j] = float32(bv_buffer[j]);

              h5_status = H5Dwrite(dset_id2, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
                if (io_log) fprintf(log_fptr, "H5Dwrite boundary value: %"ISYM"\n", h5_status);
                assert( h5_status != h5_error );

            }

#else

            if (BoundaryValue[field][dim][i] != NULL) {

	      for (j = 0; j < size; j++)
		buffer[j] = float32(BoundaryValue[field][dim][i][j]);

 
              h5_status = H5Dwrite(dset_id2, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
                if (io_log) fprintf(log_fptr, "H5Dwrite boundary value: %"ISYM"\n", h5_status);
                assert( h5_status != h5_error );
 
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
          if (io_log) fprintf(log_fptr, "H5Dclose 1: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Dclose(dset_id2);
          if (io_log) fprintf(log_fptr, "H5Dclose 2: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Sclose(mem_dsp_id);
          if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
        h5_status = H5Sclose(file_dsp_id);
          if (io_log) fprintf(log_fptr,"H5Sclose file_dsp: %"ISYM"\n", h5_status);
          assert( h5_status != h5_error );
 
      }  // end of loop over dims
 
      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      file_status = FileMover(hdfname);
        assert( file_status == 0 );
 
      if (io_log) fclose(log_fptr);
 
  }
 
  return SUCCESS;
 
}

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
#include <hdf5.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define GRID 512
#define NGRIDS 512
#define FAIL -1

#include "macros_and_parameters.h"

// HDF5 prototypes

// #include "extern_hdf5.h"

int CrackHierarchyFile(FILE *fptr, int *Grid, float GridLeftEdge[], float GridRightEdge[]);
void pcol(float *x, int n, int m, FILE *log_fptr);


int main(int argc, char **argv)
{

  hid_t       in_file_id, in_dset_id, in_attr_id;
  hid_t       in_mem_dsp_id, in_file_dsp_id;

  hid_t       file_id, dset_id;
  hid_t       file_type_id, mem_type_id, attr_type;
  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       attr_id, attr_dsp_id;

  hsize_t     max_dims[3];
  hsize_t     input_dims[3];
  hsize_t     output_dims[3];

  hsize_t     dbuff_size;
  hsize_t     gbuff_size;

  herr_t      h5_status;
  herr_t      h5_error = -1;

/*
  hsize_t     mem_stride, mem_count;
  hsize_t     file_stride, file_count;

  hssize_t    mem_offset;
  hssize_t    file_offset;
*/

  FILE *log;
  char pid[5];

  int i, j, k, m, n;
  int ijk;

  int rank, ngrids, gridcounter;
  int ndims;

  float Left[3];
  float Right[3];
  float eps = 1.0e-06;

  int StartIndex[3], EndIndex[3];
  int Dim[3],BigDim[3];

  FILE *fptr;
  int Grid;
  float GridLeftEdge[3];
  float GridRightEdge[3];

  char dset_label[25];
  char dset_units[25];
  char dset_format[25];
  char dset_geometry[25];

  ngrids = NGRIDS;

  char *logname = new char[5+1];
  strcpy(logname, "GGlog");
  log = fopen(logname, "w");

  char *GridName = new char[strlen(argv[1])+1];
  char *FieldName = new char[80+1];

  GridName = argv[1];
  FieldName = argv[2];

  char *GlueFile = new char[strlen(FieldName)+1];
  strcpy(GlueFile,FieldName);

  char *Hierarchy = new char[strlen(GridName)+10+1];
  char *DumpName = new char[strlen(GridName)+5+4+1];

  strcpy(Hierarchy,GridName);
  strcat(Hierarchy,".hierarchy");

  fptr = fopen(Hierarchy, "r");

  printf("Open %s\n",Hierarchy);

  BigDim[0] = GRID;
  BigDim[1] = GRID;
  BigDim[2] = GRID;

  output_dims[0] = GRID;
  output_dims[1] = GRID;
  output_dims[2] = GRID;
  gbuff_size = output_dims[0] * output_dims[1] * output_dims[2];

  float *output_buffer = new float[gbuff_size];

  mem_type_id = HDF5_R4;

  attr_type = H5Tcopy(H5T_C_S1);
              H5Tset_size(attr_type, 80);

  rank = 3;

  fprintf(log, "NumberOfGrids = %d\n", ngrids);

  for(gridcounter = 0; gridcounter < ngrids; gridcounter++)
  {

    n=gridcounter;

    Grid = n;

    ijk = CrackHierarchyFile(fptr, &Grid, GridLeftEdge, GridRightEdge);

    Left[0] = GridLeftEdge[0];
    Left[1] = GridLeftEdge[1];
    Left[2] = GridLeftEdge[2];

    Right[0] = GridRightEdge[0];
    Right[1] = GridRightEdge[1];
    Right[2] = GridRightEdge[2];

    StartIndex[0] = (int)(GRID * Left[0] + eps);
    EndIndex[0] = (int)(GRID * Right[0] + eps) - 1;
    StartIndex[1] = (int)(GRID * Left[1] + eps); 
    EndIndex[1] = (int)(GRID * Right[1] + eps) - 1;
    StartIndex[2] = (int)(GRID * Left[2] + eps); 
    EndIndex[2] = (int)(GRID * Right[2] + eps) - 1;

    for (m = 0; m < rank; m++)
      fprintf(log, "Grid %d    Left   %8.2f   Right  %8.2f\n", gridcounter, Left[m], Right[m]);

    for (m = 0; m < rank; m++)
      fprintf(log, "Grid %d    [%d,%d]\n", gridcounter, StartIndex[m], EndIndex[m]);

    sprintf(pid, "%4.4d", n+1);
    strcpy(DumpName, GridName);
    strcat(DumpName, ".grid");
    strcat(DumpName, pid);
    fprintf(log, "Input file %s\n", DumpName);

    dbuff_size = 1;
    file_type_id = HDF5_FILE_R4;

    in_file_id = H5Fopen(DumpName, H5F_ACC_RDONLY, H5P_DEFAULT);
      assert( in_file_id != h5_error );
    in_dset_id = H5Dopen(in_file_id, FieldName);
      assert( in_dset_id != h5_error );
    in_file_dsp_id = H5Dget_space(in_dset_id);
      assert( in_file_dsp_id != h5_error );
    ndims = H5Sget_simple_extent_dims(in_file_dsp_id, input_dims, max_dims);
      assert( ndims != h5_error );

    for (i = 0; i < ndims; i++)
    {
      dbuff_size = dbuff_size * input_dims[i];
    }
    fprintf(log, "dbuff_size %d [%Ld,%Ld,%Ld]\n", (int) dbuff_size, (int) input_dims[0], (int) input_dims[1], (int) input_dims[2]);

    in_mem_dsp_id = H5Screate_simple(1, &dbuff_size, NULL);
      assert ( in_mem_dsp_id != h5_error );

    float *input_buffer = new float[dbuff_size];

    h5_status = H5Dread(in_dset_id, mem_type_id, in_mem_dsp_id, in_file_dsp_id, H5P_DEFAULT, input_buffer);
      assert( h5_status != h5_error );


    if ( n == 0 )
    {

      in_attr_id = H5Aopen_name(in_dset_id,"Label");
        assert( in_attr_id != h5_error );
      h5_status = H5Aread(in_attr_id, attr_type, dset_label);
        assert( h5_status != h5_error );
      h5_status = H5Aclose(in_attr_id);
        assert( h5_status != h5_error );

      fprintf(log, "Dset_label %s\n", dset_label);

      in_attr_id = H5Aopen_name(in_dset_id,"Units");
        assert( in_attr_id != h5_error );
      h5_status = H5Aread(in_attr_id, attr_type, dset_units);
        assert( h5_status != h5_error );
      h5_status = H5Aclose(in_attr_id);
        assert( h5_status != h5_error );

      fprintf(log, "Dset_units %s\n", dset_units);

      in_attr_id = H5Aopen_name(in_dset_id,"Format");
        assert( in_attr_id != h5_error );
      h5_status = H5Aread(in_attr_id, attr_type, dset_format);
        assert( h5_status != h5_error );
      h5_status = H5Aclose(in_attr_id);
        assert( h5_status != h5_error );

      fprintf(log, "Dset_format %s\n", dset_format);

      in_attr_id = H5Aopen_name(in_dset_id,"Geometry");
        assert( in_attr_id != h5_error );
      h5_status = H5Aread(in_attr_id, attr_type, dset_geometry);
        assert( h5_status != h5_error );
      h5_status = H5Aclose(in_attr_id);
        assert( h5_status != h5_error );

      fprintf(log, "Dset_geometry %s\n", dset_geometry);

    }


    h5_status = H5Sclose(in_mem_dsp_id);
      assert( h5_status != h5_error );
    h5_status = H5Sclose(in_file_dsp_id);
      assert( h5_status != h5_error );
    h5_status = H5Dclose(in_dset_id);
      assert( h5_status != h5_error );
    h5_status = H5Fclose(in_file_id);
      assert( h5_status != h5_error );

//  pcol(input_buffer, ((int) dbuff_size), 10, log);

    Dim[0] = (EndIndex[0]-StartIndex[0]) + 1;
    Dim[1] = (EndIndex[1]-StartIndex[1]) + 1;
    Dim[2] = (EndIndex[2]-StartIndex[2]) + 1;

    for (k = StartIndex[2]; k <= EndIndex[2]; k++)
      for (j = StartIndex[1]; j <= EndIndex[1]; j++)
        for (i = StartIndex[0]; i <= EndIndex[0]; i++)
          output_buffer[k*BigDim[0]*BigDim[1] + j*BigDim[0] + i] =
            input_buffer[(k-StartIndex[2])*Dim[0]*Dim[1] +
                         (j-StartIndex[1])*Dim[0] +
                         (i-StartIndex[0])];

/*
    mem_offset = 0;
    mem_stride = 1;
    mem_count = dbuff_size;
    mem_dsp_id = H5Screate_simple(1, &dbuff_size, NULL);
      assert( mem_dsp_id != h5_error );

    file_stride = 1;
    file_count = dbuff_size;

    h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
      assert( h5_status != h5_error );
    h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &file_offset, &file_stride, &file_count, NULL);
      assert( h5_status != h5_error );
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, input_buffer);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(mem_dsp_id);
      assert( h5_status != h5_error );

    file_offset = file_offset + dbuff_size;
*/

    delete [] input_buffer;

  } // End of loop over hierarchy

// Output grid

  mem_dsp_id = H5Screate_simple(1, &gbuff_size, NULL);
    assert( mem_dsp_id != h5_error );
  file_dsp_id = H5Screate_simple(3, output_dims, NULL);
    assert( file_dsp_id != h5_error );
  file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert( file_id != h5_error );
  dset_id = H5Dcreate(file_id, GlueFile, file_type_id, file_dsp_id, H5P_DEFAULT);
    assert( dset_id != h5_error );
  h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, output_buffer);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(file_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(mem_dsp_id);
    assert( h5_status != h5_error );


  attr_dsp_id = H5Screate(H5S_SCALAR);
    assert( attr_dsp_id != h5_error );
  attr_id = H5Acreate(dset_id, "Label", attr_type,  attr_dsp_id, H5P_DEFAULT);
    assert( attr_id != h5_error );
  h5_status = H5Awrite(attr_id, attr_type, dset_label);
    assert( h5_status != h5_error );
  h5_status = H5Aclose(attr_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(attr_dsp_id);
    assert( h5_status != h5_error );

  attr_dsp_id = H5Screate(H5S_SCALAR);
    assert( attr_dsp_id != h5_error );
  attr_id = H5Acreate(dset_id, "Units", attr_type,  attr_dsp_id, H5P_DEFAULT);
    assert( attr_id != h5_error );
  h5_status = H5Awrite(attr_id, attr_type, dset_units);
    assert( h5_status != h5_error );
  h5_status = H5Aclose(attr_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(attr_dsp_id);
    assert( h5_status != h5_error );

  attr_dsp_id = H5Screate(H5S_SCALAR);
    assert( attr_dsp_id != h5_error );
  attr_id = H5Acreate(dset_id, "Format", attr_type,  attr_dsp_id, H5P_DEFAULT);
    assert( attr_id != h5_error );
  h5_status = H5Awrite(attr_id, attr_type, dset_format);
    assert( h5_status != h5_error );
  h5_status = H5Aclose(attr_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(attr_dsp_id);
    assert( h5_status != h5_error );

  attr_dsp_id = H5Screate(H5S_SCALAR);
    assert( attr_dsp_id != h5_error );
  attr_id = H5Acreate(dset_id, "Geometry", attr_type,  attr_dsp_id, H5P_DEFAULT);
    assert( attr_id != h5_error );
  h5_status = H5Awrite(attr_id, attr_type, dset_geometry);
    assert( h5_status != h5_error );
  h5_status = H5Aclose(attr_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(attr_dsp_id);
    assert( h5_status != h5_error );

  h5_status = H5Dclose(dset_id);
    assert( h5_status != h5_error );
  h5_status = H5Fclose(file_id);
    assert( h5_status != h5_error );

  fclose(log);

}





void pcol(float *x, int n, int m, FILE *log_fptr)
{

  int nrow,mrow;
  int i,j;

#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  if (io_log)
  {
    nrow = n/m;
    mrow = n - nrow * m;
    if( mrow > 0 )
    {
      nrow = nrow+1;
    }

    fprintf(log_fptr, "\n");

    for(j=0;j<n;j=j+m)
    {
      for(i=j;i<min(j+m,n);i++)
      {
        fprintf(log_fptr, "%12.4e", x[i]);
      }
      fprintf(log_fptr, "\n");
    }
  }

}

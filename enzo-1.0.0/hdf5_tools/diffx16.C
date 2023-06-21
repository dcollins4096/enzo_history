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
#include <math.h>
#include <assert.h>

void pcol(float *x, int n, int m, FILE *log_fptr);


int main(int argc, char *argv[])
{

  hid_t       file_id, dset_id, attr_id;
  hid_t       mem_dsp_id, file_dsp_id, attr_dsp_id;
  hid_t       mem_type_id;

  hid_t       dsp_id;
  hid_t       typ_id;

  hsize_t     size;
  hsize_t     dims[4];

  hsize_t     xdims[4];
  hsize_t     maxdims[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

//htri_t      h5_equal;
  

  int i,j,k,m;
  int ndims;
  int *ibuff_1;
  int *ibuff_2;
  float *buff_1;
  float *buff_2;

  int *idiff;
  float *fdiff;

  char *file1;
  char *file2;
  char *datax;

  if ( argc != 4)
  {
    printf("Usage: differ file1 file2 datasetname\n");
    return( -1 );
  }

  file1 = argv[1];
  file2 = argv[2];
  datax = argv[3];

  file_id = H5Fopen(file1, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_id = H5Dopen(file_id, datax);

  dsp_id = H5Dget_space(dset_id);
  typ_id = H5Dget_type(dset_id);
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  if( H5Tequal( typ_id, H5T_IEEE_F32BE) )
  {
    printf("type equality ok\n");
    printf("data type is H5T_IEEE_F32BE\n");
  }

  if( H5Tequal( typ_id, H5T_STD_I32BE) )
  {
    printf("type equality ok\n");
    printf("data type is H5T_STD_I32BE\n");
  } 

  printf("Datum size %d\n", (int) H5Tget_size( typ_id ) );

  H5T_order_t byte_order;

  byte_order = H5Tget_order( typ_id );

  if( byte_order == H5T_ORDER_LE )
    printf("Little-Endian\n");

  if( byte_order == H5T_ORDER_BE )
    printf("Big-Endian\n");

  if( byte_order == H5T_ORDER_NONE )
    printf("No specific byte order (strings)\n");

  if( byte_order == H5T_ORDER_ERROR )
    printf("H5Tget_order FAILED\n");

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
  mem_dsp_id = H5Screate_simple(1, &size, NULL);

  if ( H5Tequal(typ_id, H5T_IEEE_F32BE) )
  {
    buff_1 = new float[(int) size];
    mem_type_id = H5T_NATIVE_FLOAT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, buff_1);
    printf("float read status %d\n", (int) h5_status);
  }


  if ( H5Tequal(typ_id, H5T_STD_I32BE) )
  {
    ibuff_1 = new int[(int) size];
    mem_type_id = H5T_NATIVE_INT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, ibuff_1);
    printf("int read status %d\n", (int) h5_status);
  }


  h5_status = H5Sclose(dsp_id);
  h5_status = H5Tclose(typ_id);
  h5_status = H5Sclose(mem_dsp_id);
  h5_status = H5Sclose(file_dsp_id);
  h5_status = H5Dclose(dset_id);
  h5_status = H5Fclose(file_id);


  file_id = H5Fopen(file2, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_id = H5Dopen(file_id, datax);

  dsp_id = H5Dget_space(dset_id);
  typ_id = H5Dget_type(dset_id);
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
  mem_dsp_id = H5Screate_simple(1, &size, NULL);

  if ( H5Tequal(typ_id, H5T_IEEE_F32BE) )
  {
    buff_2 = new float[(int) size];
    mem_type_id = H5T_NATIVE_FLOAT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, buff_2);
    printf("float read status %d\n", (int) h5_status);
  }


  if ( H5Tequal(typ_id, H5T_STD_I32BE) )
  {
    ibuff_2 = new int[(int) size];
    mem_type_id = H5T_NATIVE_INT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, ibuff_2);
    printf("int read status %d\n", (int) h5_status);
  }


  h5_status = H5Sclose(dsp_id);
//  h5_status = H5Tclose(typ_id);
  h5_status = H5Sclose(mem_dsp_id);
  h5_status = H5Sclose(file_dsp_id);
  h5_status = H5Dclose(dset_id);
  h5_status = H5Fclose(file_id);


  if ( H5Tequal(typ_id, H5T_IEEE_F32BE) )
  {

    fdiff = new float[(int) size];
    float eps = 1.0e-15;

    for (i = 0; i < size; i++)
    {
      fdiff[i] = fabs((buff_1[i] - buff_2[i])/(buff_1[i] + buff_2[i] + eps));
    }


    int nfail = 0;
    int nchecks = 0;
    int nprec = 0;

    for (i = 0; i < size; i++)
    {
      if( buff_1[i] != buff_2[i] )
      {
        nfail++;
      }
      nchecks++;
    }

    printf("N checks %d\n",nchecks);
    printf("N fails %d\n",nfail);

    if ( nfail != 0 )
    {

      FILE *out;
      out = fopen("dump","w");

      for (i = 0; i < size; i++)
      {
        if( fdiff[i] > 0.00001 )
        {
          fprintf(out,"  %8d  %16.4e  %16.8f  %16.8f\n", i, fdiff[i],buff_1[i],buff_2[i]);
          nprec++;
        }
      }
      fclose(out);
    }
    printf("Failures above tolerance %d\n",nprec);
    H5Tclose(typ_id);

  }  // float only for now

  m = 0;

  for(k=0; k<16; k++)
  {
    printf("K=%d\n",m/256);
    printf("    \n");
    for(j=0; j<16; j++)
    {
      printf("%2d: ",j);
      for(i=0; i<16; i++)
      {
        if( fdiff[m] > 0.00001 )
        {
          printf("x");
        }
        else
        {
          printf(" ");
        }
        m++;
      }
      printf(":\n");
    }
    printf("    \n");
    printf("\n");
    printf("\n");
  }
}



void pcol(float *x, int n, int m, FILE *log_fptr)
{

int nrow,mrow;
int i,j;

  int         io_log = 1;

if (io_log)
{
  nrow = n/m;
  mrow = n - nrow * m;
  if( mrow > 0 )
  {
    nrow = nrow+1;
  }

  fprintf(log_fptr,"\n");

  for(j=0;j<n;j=j+m)
  {
    for(i=j;i<min(j+m,n);i++)
    {
      fprintf(log_fptr, "%12.3e", x[i]);
    }
    fprintf(log_fptr,"\n");
  }
}

}

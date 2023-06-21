#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))


#include <hdf5.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// HDF5 prototypes

// #include "extern_hdf5.h"

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

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int i;
  int ndims;
  int *ibuff_1;
  int *ibuff_2;
  float *buff_1;
  float *buff_2;
  double *dbuff_1;
  double *dbuff_2;

  int *idiff;
  float *fdiff;
  double *ddiff;

  char *file1;
  char *file2;
  char *datax;
  char *datay;

  if ( argc != 5)
  {
    printf("Usage: differ file1 file2 dataset1 dataset2\n");
    return( -1 );
  }

  file1 = argv[1];
  file2 = argv[2];
  datax = argv[3];
  datay = argv[4];

  file_id = H5Fopen(file1, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert( file_id != h5_error );
  dset_id = H5Dopen(file_id, datax);
    assert( dset_id != h5_error );

  dsp_id = H5Dget_space(dset_id);
    assert( dsp_id != h5_error );
  typ_id = H5Dget_type(dset_id);
    assert( typ_id != h5_error );
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

  if( H5Tequal( typ_id, H5T_IEEE_F64BE) )
  {
    printf("type equality ok\n");
    printf("data type is H5T_IEEE_F64BE\n");
  }

  printf("Datum size %d\n", (int) H5Tget_size( typ_id ) );

  H5T_order_t byte_order;

  byte_order = H5Tget_order( typ_id );
    assert( byte_order != h5_error );

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
    assert( file_dsp_id != h5_error );
  mem_dsp_id = H5Screate_simple(1, &size, NULL);
    assert( mem_dsp_id != h5_error );

  if ( H5Tequal(typ_id, H5T_IEEE_F32BE) )
  {
    buff_1 = new float[(int) size];
    mem_type_id = H5T_NATIVE_FLOAT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, buff_1);
      printf("float read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal(typ_id, H5T_IEEE_F64BE) )
  {
    dbuff_1 = new double[(int) size];
    mem_type_id = H5T_NATIVE_DOUBLE;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, dbuff_1);
      printf("double read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal(typ_id, H5T_STD_I32BE) )
  {
    ibuff_1 = new int[(int) size];
    mem_type_id = H5T_NATIVE_INT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, ibuff_1);
      printf("int read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }


  h5_status = H5Sclose(dsp_id);
    assert( h5_status != h5_error );
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


  file_id = H5Fopen(file2, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert( file_id != h5_error );
  dset_id = H5Dopen(file_id, datay);
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

  if ( H5Tequal(typ_id, H5T_IEEE_F32BE) )
  {
    buff_2 = new float[(int) size];
    mem_type_id = H5T_NATIVE_FLOAT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, buff_2);
      printf("float read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal(typ_id, H5T_IEEE_F64BE) )
  {
    dbuff_2 = new double[(int) size];
    mem_type_id = H5T_NATIVE_DOUBLE;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, dbuff_2);
      printf("double read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }

  if ( H5Tequal(typ_id, H5T_STD_I32BE) )
  {
    ibuff_2 = new int[(int) size];
    mem_type_id = H5T_NATIVE_INT;
    h5_status = H5Dread(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, ibuff_2);
      printf("int read status %d\n", (int) h5_status);
      assert( h5_status != h5_error );
  }


  h5_status = H5Sclose(dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(mem_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Sclose(file_dsp_id);
    assert( h5_status != h5_error );
  h5_status = H5Dclose(dset_id);
    assert( h5_status != h5_error );
  h5_status = H5Fclose(file_id);
    assert( h5_status != h5_error );


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

  }  // float only for now


  if ( H5Tequal(typ_id, H5T_IEEE_F64BE) )
  {

    ddiff = new double[(int) size];
    double eps = 1.0e-15;

    for (i = 0; i < size; i++)
    {
      ddiff[i] = fabs((dbuff_1[i] - dbuff_2[i])/(dbuff_1[i] + dbuff_2[i] + eps));
    }


    int nfail = 0;
    int nchecks = 0;
    int nprec = 0;

    for (i = 0; i < size; i++)
    {
      if( dbuff_1[i] != dbuff_2[i] )
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
        if( ddiff[i] > 0.0000000001 )
        {
          fprintf(out,"  %8d  %16.4e  %16.8e  %16.8e\n", i, ddiff[i],dbuff_1[i],dbuff_2[i]);
          nprec++;
        }
      }
      fclose(out);
    }
    printf("Failures above tolerance %d\n",nprec);

  }


  if ( H5Tequal(typ_id, H5T_STD_I32BE) )
  {

    idiff = new int[(int) size];

    for (i = 0; i < size; i++)
    {
      idiff[i] = ibuff_1[i] - ibuff_2[i];
    }


    int nfail = 0;
    int nchecks = 0;
    int nprec = 0;

    for (i = 0; i < size; i++)
    {
      if( ibuff_1[i] != ibuff_2[i] )
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
        if( idiff[i] != 0 )
        {
          fprintf(out,"  %8d  %20d  %20d  %20d\n", i, idiff[i],ibuff_1[i],ibuff_2[i]);
          nprec++;
        }
      }
      fclose(out);
    }
    printf("Failures above tolerance %d\n",nprec);

  }

  h5_status = H5Tclose(typ_id);
    assert( h5_status != h5_error );

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

#include <mpi.h>
#include <hdf5.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define MAXCPU 1024
#define GRID 16

#include "macros_and_parameters.h"

// HDF5 prototypes

#include "extern_hdf5.h"

void pcol(float *x, int n, int m, FILE *log_fptr);


int main(int argc, char **argv)
{

  hid_t       in_file_id, in_dset_id, in_type_id;
  hid_t       in_mem_dsp_id, in_file_dsp_id;

  hid_t       file_id, dset_id, type_id;
  hid_t       mem_dsp_id, file_dsp_id;

  hsize_t     max_dims[3];
  hsize_t     input_dims[3];
  hsize_t     output_dims[3];

  hsize_t     dbuff_size;
  hsize_t     gbuff_size;

  herr_t      h5_status;
  herr_t      h5_error = -1;

  hsize_t     mem_stride, mem_count;
  hsize_t     file_stride, file_count;

  hssize_t    mem_offset;
  hssize_t    file_offset;

  FILE *log;
  char pid[5];

  int i, j, k, m, n;
  int a, b, c;
  int ii, jj, kk, mx, ijk;

  int dim, rank, ngrids, gridcounter;
  int ndims;
  int ntiles;
  int mpi_size, mpi_rank;
  int mpi_layout[3] = {0,0,0};
  int enzo_layout[3] = {0,0,0};

  float Left[3];
  float Right[3];
  float dx, dy, dz;
  float eps = 1.0e-06;

  int StartIndex[3], EndIndex[3];
  int Dim[3],BigDim[3];
  int *GridDims[3], *Index[3];
  int TempDims[3];
  float LeftEdge[3];
  float RightEdge[3];



  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

  ntiles = 4;

  char *pgname = new char[strlen(pid)+6+1];
  strcpy(pgname, "GGlog");
  log = fopen(pgname, "w");

  char *GlueFile = "Glued";
  char *GridName = "E_38.4_256cube0001";
  char *FieldName = "Temperature";

  char *DumpName = new char[strlen(GridName)+5+4+1];

  BigDim[0] = GRID;
  BigDim[1] = GRID;
  BigDim[2] = GRID;

  output_dims[0] = GRID;
  output_dims[1] = GRID;
  output_dims[2] = GRID;
  gbuff_size = output_dims[0] * output_dims[1] * output_dims[2];

  float *output_buffer = new float[gbuff_size];

  type_id = HDF5_R4;

  rank = 3;
  ngrids = 1;

  MPI_Dims_create(ntiles, rank, mpi_layout);

  if (rank == 3 && ntiles == 64)
  {
    for (dim = 0; dim < rank; dim++)
    {
      mpi_layout[dim] = 4;
    }
    printf("Ntiles = 64 ==> coerced to 4**3\n");
  }

  for (dim = 0; dim < rank; dim++)
  {
    enzo_layout[dim] = mpi_layout[rank-1-dim];
    ngrids *= enzo_layout[dim];
  }


  for (dim = 0; dim < rank; dim++)
  {
    GridDims[dim] = new int[enzo_layout[dim]];
    Index[dim] = new int[enzo_layout[dim]];

    float ExactDims = float(GRID)/float(enzo_layout[dim]);
    float ExactCount = 0;
    int Disp = 0;

    for (i = 0; i < enzo_layout[dim]; i++)
    {
      ExactCount += ExactDims;
      if (dim == 0)
        GridDims[dim][i] = nint(ExactCount*0.5)*2 - Disp;
      else
        GridDims[dim][i] = nint(ExactCount) - Disp;
      Index[dim][i] = Disp;
      Disp += GridDims[dim][i];
    }
  }

  for (dim = 0; dim < rank; dim++)
  {
    printf("GridDims[%d]: ",dim);
    for (i = 0; i < enzo_layout[dim]; i++)
    {
      printf(" %d",GridDims[dim][i]);
    }
    printf("\n");
  }

  for (dim = 0; dim < rank; dim++)
  {
    printf("Index[%d]: ",dim);
    for (i = 0; i < enzo_layout[dim]; i++)
    {
      printf(" %d",Index[dim][i]);
    }
    printf("\n");
  }



  fprintf(log,"NumberOfGrids = %d\n",ngrids);
  fprintf(log,"ENZO_layout %d %d %d\n",enzo_layout[0],enzo_layout[1],enzo_layout[2]);

  dx = 1.0 / ((float) enzo_layout[0]);
  dy = 1.0 / ((float) enzo_layout[1]);
  dz = 1.0 / ((float) enzo_layout[2]);

  a=enzo_layout[0];
  b=enzo_layout[1];
  c=enzo_layout[2];

  gridcounter = 0;

  file_offset = 0;

  for (kk = 0; kk < enzo_layout[2]; kk++)
    for (jj = 0; jj < enzo_layout[1]; jj++)
      for (ii = 0; ii < enzo_layout[0]; ii++)
      {

        n=gridcounter;

        for (dim = 0; dim < rank; dim++)
        {
          ijk = (dim == 0) ? ii : ((dim == 1) ? jj : kk);
          TempDims[dim] = GridDims[dim][ijk];
          LeftEdge[dim] = float(Index[dim][ijk])/float(GRID);
          RightEdge[dim] = float(Index[dim][ijk]+TempDims[dim])/float(GRID);
        }

        fprintf(log,"Grid %d Edges %10.4f %10.4f\n",n,LeftEdge[0],RightEdge[0]);
        fprintf(log,"Grid %d Edges %10.4f %10.4f\n",n,LeftEdge[1],RightEdge[1]);
        fprintf(log,"Grid %d Edges %10.4f %10.4f\n",n,LeftEdge[2],RightEdge[2]);

      // rank to coordinate
        m = n;
        i = m/(b*c);
        m = m%(b*c);
        j = m/c;
        m = m%c;
        k = m;
      // coordinate to rank check
        m = ((i*b*c) + j*c) + k;
        mx = ((ii*b*c) + jj*c) + kk;

        fprintf(log,"Grid %d  {%d %d %d} m=%d mx=%d\n",n,i,j,k,m,mx);

        Left[0] =  dx * (float) ii;
        Right[0] = dx * (float) (ii+1);
        Left[1] =  dy * (float) jj;
        Right[1] = dy * (float) (jj+1);
        Left[2] =  dz * (float) kk;
        Right[2] = dz * (float) (kk+1);

        StartIndex[0] = (int)(GRID * Left[0] + eps);
        EndIndex[0] = (int)(GRID * Right[0] + eps) - 1;
        StartIndex[1] = (int)(GRID * Left[1] + eps); 
        EndIndex[1] = (int)(GRID * Right[1] + eps) - 1;
        StartIndex[2] = (int)(GRID * Left[2] + eps); 
        EndIndex[2] = (int)(GRID * Right[2] + eps) - 1;

        for (m = 0; m < rank; m++)
          fprintf(log,"Grid %d    Left   %8.2f   Right  %8.2f\n",gridcounter,Left[m],Right[m]);

        for (m = 0; m < rank; m++)
          fprintf(log,"Grid %d    [%d,%d]\n",gridcounter,StartIndex[m],EndIndex[m]);

/*
        in_mem_offset = 0;
        in_mem_stride = 1;
        in_mem_count = dbuff_size;
        in_mem_dsp_id = H5Screate_simple(1, &dbuff_size, NULL);
*/

        sprintf(pid, "%4.4d", n+1);
        strcpy(DumpName, GridName);
        strcat(DumpName, ".grid");
        strcat(DumpName, pid);
        fprintf(log,"Input file %s\n", DumpName);

        dbuff_size = 1;
        in_type_id = HDF5_R4;
        in_file_id = H5Fopen(DumpName, H5F_ACC_RDONLY, H5P_DEFAULT);
        in_dset_id = H5Dopen(in_file_id, FieldName);
        in_file_dsp_id = H5Dget_space(in_dset_id);
        ndims = H5Sget_simple_extent_dims(in_file_dsp_id, input_dims, max_dims);
        for (i = 0; i < ndims; i++)
        {
          dbuff_size = dbuff_size * input_dims[i];
        }
        fprintf(log,"dbuff_size %d [%Ld,%Ld,%Ld]\n",(int) dbuff_size,(int) input_dims[0],(int) input_dims[1],(int) input_dims[2]);
/*
        Dim[2] = input_dims[0];
        Dim[1] = input_dims[1];
        Dim[0] = input_dims[2];
*/
        in_mem_dsp_id = H5Screate_simple(1, &dbuff_size, NULL);

        float *input_buffer = new float[dbuff_size];

        h5_status = H5Dread(in_dset_id, in_type_id, in_mem_dsp_id, in_file_dsp_id, H5P_DEFAULT, input_buffer);

        assert( h5_status != h5_error );

        h5_status = H5Sclose(in_mem_dsp_id);
        h5_status = H5Sclose(in_file_dsp_id);
        h5_status = H5Dclose(in_dset_id);
        h5_status = H5Fclose(in_file_id);

//        pcol(input_buffer, ((int) dbuff_size),10,log);


     if ( n == 0 )
     {
        StartIndex[0] = 0;
        StartIndex[1] = 0;
        StartIndex[2] =0;
        EndIndex[0] = 15;
        EndIndex[1] = 7;
        EndIndex[2] = 7;
     }
     if ( n == 1 )
     {
        StartIndex[0] = 0;
        StartIndex[1] = 8;
        StartIndex[2] = 8;
        EndIndex[0] = 15;
        EndIndex[1] = 15;
        EndIndex[2] = 15;
     }
     if ( n == 2 )
     {
        StartIndex[0] = 0;
        StartIndex[1] = 0;
        StartIndex[2] = 8;
        EndIndex[0] = 15;
        EndIndex[1] = 7;
        EndIndex[2] = 15;
     }
     if ( n == 3 )
     {
        StartIndex[0] = 0;
        StartIndex[1] = 8;
        StartIndex[2] = 0;
        EndIndex[0] = 15;
        EndIndex[1] = 15;
        EndIndex[2] = 7;
     }

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

        file_stride = 1;
        file_count = dbuff_size;

        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &file_offset, &file_stride, &file_count, NULL);
        h5_status = H5Dwrite(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, input_buffer);

        h5_status = H5Sclose(mem_dsp_id);

        file_offset = file_offset + dbuff_size;
*/

        delete [] input_buffer;

        gridcounter++;

      }

  mem_dsp_id = H5Screate_simple(1, &gbuff_size, NULL);
  file_dsp_id = H5Screate_simple(3, output_dims, NULL);
  file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dset_id = H5Dcreate(file_id, FieldName, type_id, file_dsp_id, H5P_DEFAULT);
  h5_status = H5Dwrite(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, output_buffer);
  h5_status = H5Sclose(file_dsp_id);
  h5_status = H5Sclose(mem_dsp_id);
  h5_status = H5Dclose(dset_id);
  h5_status = H5Fclose(file_id);

  fclose(log);

  MPI_Finalize();
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

  fprintf(log_fptr,"\n");

  for(j=0;j<n;j=j+m)
  {
    for(i=j;i<min(j+m,n);i++)
    {
      fprintf(log_fptr, "%12.4e", x[i]);
    }
    fprintf(log_fptr,"\n");
  }
}

}

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
#include <mpi.h>
#include <hdf5.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXCPU 1024

#define XDIM 16
#define YDIM 16
#define ZDIM 16

// HDF5 prototypes

#include "extern_hdf5.h"

int main(int argc, char **argv)
{

  hid_t       file_id, dset_id, type_id, attr_id;
  hid_t       mem_dsp_id, file_dsp_id, attr_dsp_id;
  hid_t       file_acc_template;
  hid_t       xfer_prop_list;

  hsize_t     dims[3];
  herr_t      h5_status;

  hsize_t     mem_stride, mem_count, attr_count;
  hsize_t     slab_stride[3], slab_count[3];
  hsize_t     Slab_Rank;
  hsize_t     Slab_Dims[2];

  hssize_t    mem_offset;
  hssize_t    slab_offset[3];

  int ncpu, jcpu;

  int i, j, k, m, n;
  int a, b, c;
  int ii, jj, kk;
  int ic, ipc, jpc, ppc, iblk;

  int thisnode, prevnode, nextnode, ltype, rtype, stype, ier;
  
  int dim, rank, ngrids, gridcounter;

  int mpi_size, mpi_rank;

  int mpi_layout[3] = {0,0,0};
  int enzo_layout[3] = {0,0,0};
  float Left[3];
  float Right[3];
  float dx, dy, dz;

  float GridLeft[MAXCPU][3];
  float GridRight[MAXCPU][3];

  int dbuff_size;
  int TotalParticleCount;
  int NumberOfParticles;
  int ParticleRank;

  double Start_Wall_Time, End_Wall_Time, Total_Wall_Time;

  MPI_Request req1[100];
  MPI_Request req2[100];
  MPI_Status  stat[100];

  int io_log = 0;

// 

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

  ncpu = mpi_size;
  jcpu = mpi_rank;


  FILE *log;
  char pid[5];
  sprintf(pid, "%4.4d", jcpu);

  char *pgname = new char[strlen(pid)+6+1];
  strcpy(pgname, "PGlog.");
  strcat(pgname, pid);

  log = fopen(pgname, "a");

  char *PPos = new char[strlen(pid)+5+1];
  strcpy(PPos, "PPos.");
  strcat(PPos, pid);

  char *PVel = new char[strlen(pid)+5+1];
  strcpy(PVel, "PVel.");
  strcat(PVel, pid);

  char *PPro = new char[strlen(pid)+5+1];
  strcpy(PPro, "PPro.");
  strcat(PPro, pid);


//  TotalParticleCount = XDIM*YDIM*ZDIM;

  file_id = H5Fopen("ParticlePositions", H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_id = H5Dopen(file_id, "ParticlePositions");
  attr_id = H5Aopen_name(dset_id, "Component_Rank");
  h5_status = H5Aread(attr_id, H5T_NATIVE_INT, &ParticleRank);
  h5_status = H5Aclose(attr_id);
  attr_id = H5Aopen_name(dset_id, "Component_Size");
  h5_status = H5Aread(attr_id, H5T_NATIVE_INT, &TotalParticleCount);
  h5_status = H5Aclose(attr_id);
  h5_status = H5Dclose(dset_id);
  h5_status = H5Fclose(file_id);

  dbuff_size = TotalParticleCount/ncpu;

  float *buff[2];
  buff[0] = new float[dbuff_size];
  buff[1] = new float[dbuff_size];

  rank = 3;
  ngrids = 1;

  MPI_Dims_create(ncpu, rank, mpi_layout);

  for (dim = 0; dim < rank; dim++)
  {
    enzo_layout[dim] = mpi_layout[rank-1-dim];
    ngrids *= enzo_layout[dim];
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


  for (kk = 0; kk < enzo_layout[2]; kk++)
    for (jj = 0; jj < enzo_layout[1]; jj++)
      for (ii = 0; ii < enzo_layout[0]; ii++)
      {
        n=gridcounter;

        if ( n == jcpu )
        {
      // rank to coordinate
        m = n;
        i = m/(b*c);
        m = m%(b*c);
        j = m/c;
        m = m%c;
        k = m;
      // coordinate to rank check
        m = ((i*b*c) + j*c) + k;
        fprintf(log,"Grid %d  {%d %d %d}  %d\n",n,i,j,k,m);

        Left[0] =  dx * (float) ii;
        Right[0] = dx * (float) (ii+1);
        Left[1] =  dy * (float) jj;
        Right[1] = dy * (float) (jj+1);
        Left[2] =  dz * (float) kk;
        Right[2] = dz * (float) (kk+1);

        GridLeft[jcpu][0] = Left[0];
        GridLeft[jcpu][1] = Left[1];
        GridLeft[jcpu][2] = Left[2];

        GridRight[jcpu][0] = Right[0];
        GridRight[jcpu][1] = Right[1];
        GridRight[jcpu][2] = Right[2];

        for (m = 0; m < rank; m++)
          fprintf(log,"Grid %d    Left   %8.2f   Right  %8.2f\n",gridcounter,Left[m],Right[m]);

        }

        gridcounter++;

      }

    unsigned int *BitArray = NULL;
    unsigned int *BitMask = NULL;
    unsigned int BitsPerInt;
    unsigned int BitMaskSize;
    unsigned int BitMaskTrue;
    unsigned int TestBit;
    unsigned int MaskAddr;
    unsigned int WordAddr;
    unsigned int BitAddr;


    int XMask;


    BitsPerInt = 8 * sizeof(BitsPerInt);
    BitMaskSize = (TotalParticleCount/BitsPerInt)+1;

    BitMask = new unsigned int[BitMaskSize];
    BitArray = new unsigned int[BitsPerInt];

    for (i = 0; i < BitMaskSize; i++)
      BitMask[i] = 0;

    BitArray[BitsPerInt-1] = 1;
    BitMaskTrue = 1;

    for (i=BitsPerInt-1; i>0; i--)
    {
      BitArray[i-1] = 2 * BitArray[i];
      BitMaskTrue = (BitMaskTrue | BitArray[i-1]);
    }

    for (i = 0; i < BitMaskSize; i++)
      BitMask[i] = BitMaskTrue;

/*
    int *Mask = NULL; 
    Mask = new int[TotalParticleCount];
    for (i = 0; i < TotalParticleCount; i++)
      Mask[i] = 1;
*/

Start_Wall_Time = MPI_Wtime();

for ( dim = 0; dim < rank; dim++)
{

  dims[0] = 3;
  dims[1] = ncpu;
  dims[2] = TotalParticleCount/ncpu;

  type_id = H5T_NATIVE_FLOAT;

//  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);

//  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);

  file_acc_template = H5P_DEFAULT;
  file_id = H5Fopen("ParticlePositions", H5F_ACC_RDONLY, file_acc_template);
  printf("File id: %d\n",file_id);
  dset_id = H5Dopen(file_id, "ParticlePositions");
  printf("Data id: %d\n",dset_id);
  mem_dsp_id = H5Screate_simple(1, &dims[2], NULL);
  printf("Mem dsp: %d\n", mem_dsp_id);

  mem_offset = 0;
  mem_count = dims[2];
  mem_stride = 1;

  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
  printf("Mem hs: %d\n",h5_status);
  file_dsp_id = H5Screate_simple(3, dims, NULL);
  printf("File dsp: %d\n",file_dsp_id);

  slab_offset[0] = dim;   // x,y,z
  slab_stride[0] = 1;
  slab_count[0] = 1;

  slab_offset[1] = jcpu;
  slab_stride[1] = 1;
  slab_count[1] = 1;

  slab_offset[2] = 0;
  slab_stride[2] = 1;
  slab_count[2] = dims[2];

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
  printf("File hs: %d\n",h5_status);

//  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);

//  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);

  xfer_prop_list = H5P_DEFAULT;
  h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, buff[0]);
  printf("Read: %d\n",h5_status);

//  h5_status = H5Pclose(xfer_prop_list);
//  h5_status = H5Pclose(file_acc_template);
  h5_status = H5Sclose(file_dsp_id);
  h5_status = H5Sclose(mem_dsp_id);
  h5_status = H5Dclose(dset_id);
  h5_status = H5Fclose(file_id);

  if (io_log)
  {
        fprintf(log,"Proc %d Dim %d\n",jcpu,dim);
        for ( i = 0; i < dims[2]; i++ )
        {
          fprintf(log,"%8.4f",buff[0][i]);
          if ( ((i+1) % 16) == 0 )
            fprintf(log,"\n");
        }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  thisnode = jcpu;
  nextnode = (thisnode+1) % ncpu;
  prevnode = (thisnode-1+ncpu) % ncpu;

  fprintf(log,"P:  %d %d %d\n",prevnode,thisnode,nextnode);

  stype = 30000+thisnode;
  rtype = 30000+prevnode;
  ltype = 30000+nextnode;

  for (k = 0; k < ncpu; k++)
  {

  i = (k  ) % 2;
  j = (k+1) % 2;

  MPI_Barrier(MPI_COMM_WORLD);

/* Transmit right
  ier = MPI_Irecv( buff[j], dims[2], MPI_FLOAT, prevnode, rtype, MPI_COMM_WORLD, req1 );
  if (io_log) fprintf(log,"IRECV proc %d buffer[%d] err %d\n",thisnode,j,ier);

  ier = MPI_Isend( buff[i], dims[2], MPI_FLOAT, nextnode, stype, MPI_COMM_WORLD, req2 );
  if (io_log) fprintf(log,"SEND proc %d buffer[%d] err %d\n",thisnode,i,ier);
*/

/* Transmit left */
  ier = MPI_Irecv( buff[j], dims[2], MPI_FLOAT, nextnode, ltype, MPI_COMM_WORLD, req1 );
  if (io_log) fprintf(log,"IRECV proc %d buffer[%d] err %d\n",thisnode,j,ier);

  ier = MPI_Isend( buff[i], dims[2], MPI_FLOAT, prevnode, stype, MPI_COMM_WORLD, req2 );
  if (io_log) fprintf(log,"SEND proc %d buffer[%d] err %d\n",thisnode,i,ier);

  // work on buff[i]

  if (io_log) fprintf(log,"WORK on buffer[%d]\n",i);

  if (io_log)
  {
        fprintf(log,"Proc %d Dim %d K %d Buffer%d\n",jcpu,dim,k,i);
        for (int idim = 0; idim < dims[2]; idim++ )
        {
          fprintf(log,"%8.4f",buff[i][idim]);
          if ( ((idim+1) % 16) == 0 )
            fprintf(log,"\n");
        }
  }

/* Transmit right
  ipc = (jcpu-k+ncpu) % ncpu;
*/

/* Transmit left */
  ipc = (jcpu+k+ncpu) % ncpu;

  ipc = ipc * dims[2];
  jpc = ipc + dims[2] - 1;

  if (io_log) fprintf(log,"ipc = %d to %d\n",ipc, jpc);
  if (io_log) fprintf(log,"left %4.2f, right %4.2f\n", Left[dim], Right[dim]);

  for (ic = 0; ic < dims[2]; ic++)
  {
    if ( buff[i][ic] < Left[dim] || buff[i][ic] > Right[dim] )
    {
//    Mask[ipc] = 0;
      MaskAddr = ipc;
      WordAddr = MaskAddr/BitsPerInt;
      BitAddr  = MaskAddr%BitsPerInt;
      TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
      if ( TestBit != 0 )
        BitMask[WordAddr] = (BitMask[WordAddr] ^ BitArray[BitAddr]);
    }
    ipc++;
  }

/*
  if (io_log)
  {
    for (ic = 0; ic < TotalParticleCount; ic++)
    {
      fprintf(log,"%1d", Mask[ic]);
      if ( ((ic+1) % 128) == 0 )
        fprintf(log,"\n");
    }
  }
*/

  ier = MPI_Wait( req1, stat );
  if (io_log) fprintf(log,"WAIT proc %d err %d\n",thisnode,ier);

  MPI_Barrier(MPI_COMM_WORLD);

  } // End loop over ring

} // End loop over {x,y,z}

End_Wall_Time = MPI_Wtime();
Total_Wall_Time = End_Wall_Time - Start_Wall_Time;
printf("%d: %16.8f %16.8f %16.8f\n",jcpu,Start_Wall_Time,End_Wall_Time,Total_Wall_Time);

// CHOP if(0)
// CHOP {

  NumberOfParticles = 0;

  for (i = 0; i < TotalParticleCount; i++)
  {
    XMask = 1;
    MaskAddr = i;
    WordAddr = MaskAddr/BitsPerInt;
    BitAddr  = MaskAddr%BitsPerInt;
    TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
    if ( TestBit == 0 )
      XMask = 0;
    if (io_log)
    {
      fprintf(log, "%1d", XMask);
      if ( ((i+1) % 128) == 0 ) 
      {
        fprintf(log, "\n");
      }
    }
    NumberOfParticles += XMask;
  }

  fprintf(log,"Grid %d:  NumberOfParticles %d\n",jcpu,NumberOfParticles);

/*
  if (io_log)
  {
    for (i = 0; i < TotalParticleCount; i++)
    {
      fprintf(log,"%1d", Mask[i]);
      if ( ((i+1) % 128) == 0 )
        fprintf(log,"\n");
    }
  }
*/

//  Read particle data under bit mask

  float *outbuff = NULL;

  outbuff = new float[NumberOfParticles];

  for (dim = 0; dim < rank; dim++)
  {

  ppc = 0; // output particle counter

  dims[0] = 3;
  dims[1] = ncpu;
  dims[2] = TotalParticleCount/ncpu;

  type_id = H5T_NATIVE_FLOAT;

//  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);

//  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_acc_template = H5P_DEFAULT;
  file_id = H5Fopen("ParticlePositions", H5F_ACC_RDONLY, file_acc_template);

  dset_id = H5Dopen(file_id, "ParticlePositions");

  mem_dsp_id = H5Screate_simple(1, &dims[2], NULL);

  mem_offset = 0;
  mem_count = dims[2];
  mem_stride = 1;

  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);

  file_dsp_id = H5Screate_simple(3, dims, NULL);

  slab_offset[0] = dim;   // x,y,z
  slab_stride[0] = 1;
  slab_count[0] = 1;

  slab_offset[1] = jcpu;
  slab_stride[1] = 1;
  slab_count[1] = 1;

  slab_offset[2] = 0;
  slab_stride[2] = 1;
  slab_count[2] = dims[2];

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);

//  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);

//  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);

  xfer_prop_list = H5P_DEFAULT;
  h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, buff[0]);

//  h5_status = H5Pclose(xfer_prop_list);
//  h5_status = H5Pclose(file_acc_template);
  h5_status = H5Sclose(file_dsp_id);
  h5_status = H5Sclose(mem_dsp_id);
  h5_status = H5Dclose(dset_id);
  h5_status = H5Fclose(file_id);

  MPI_Barrier(MPI_COMM_WORLD);

  thisnode = jcpu;
  nextnode = (thisnode+1) % ncpu;
  prevnode = (thisnode-1+ncpu) % ncpu;

  stype = 30000+thisnode;
  rtype = 30000+prevnode;
  ltype = 30000+nextnode;

  iblk = - ( (ncpu - jcpu) % ncpu );

  for (k = 0; k < 2*ncpu; k++) // two passes 
  {

  i = (k  ) % 2;
  j = (k+1) % 2;

  MPI_Barrier(MPI_COMM_WORLD);

  // Transmit left

  ier = MPI_Irecv( buff[j], dims[2], MPI_FLOAT, nextnode, ltype, MPI_COMM_WORLD, req1 );
  if (io_log) fprintf(log,"IRECV proc %d buffer[%d] err %d\n",thisnode,j,ier);

  ier = MPI_Isend( buff[i], dims[2], MPI_FLOAT, prevnode, stype, MPI_COMM_WORLD, req2 );
  if (io_log) fprintf(log,"SEND proc %d buffer[%d] err %d\n",thisnode,i,ier);

  // work on buff[i]

  if (io_log) fprintf(log,"K = %d, IBLK = %d\n",k,iblk);

  if ( iblk > -1 && iblk < ncpu )
  {

    ipc = (jcpu+k+ncpu) % ncpu;

    ipc = ipc * dims[2];
    jpc = ipc + dims[2] - 1;

    if (io_log) fprintf(log,"  ipc = %d to %d\n",ipc,jpc);

    for (ic = 0; ic < dims[2]; ic++)
    {
      MaskAddr = ipc;
      WordAddr = MaskAddr/BitsPerInt;
      BitAddr  = MaskAddr%BitsPerInt;
      TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
      if ( TestBit == 0 )
        XMask = 0;
      if ( TestBit != 0 )
      {
        XMask = 1;
        outbuff[ppc] = buff[i][ic];
        ppc++;
      }
      ipc++;
    }

  }

  iblk++;

  ier = MPI_Wait( req1, stat );

  if (io_log) fprintf(log,"WAIT proc %d err %d\n",thisnode,ier);

  } // end loop over rings

  fprintf(log,"L %8.4f  R %8.4f\n",GridLeft[jcpu][dim],GridRight[jcpu][dim]);

  if (io_log)
  {
    for (i = 0; i < NumberOfParticles; i++)
    {
      fprintf(log,"%8.4f",outbuff[i]);
      if ( ((i+1) % 16) == 0 )
        fprintf(log,"\n");
    }
  }

  Slab_Rank = 2;
  Slab_Dims[0] = 3;
  Slab_Dims[1] = NumberOfParticles;

/*
  component_rank_attr = 3;
  component_size_attr = NumberOfParticles;
  field_rank_attr = 1;
  field_dims_attr = NumberOfParticles;
*/

// Data in memory is considered 1D, stride 1, with zero offset

  mem_stride = 1;                   // contiguous elements
  mem_count = NumberOfParticles;    // number of elements in field
  mem_offset = 0;                   // zero offset in buffer

// 1D memory model

  mem_dsp_id = H5Screate_simple(1, &mem_count, NULL);

  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);

// Data in the file is (1+Rank)D with Npart components per grid point.
// Offset[0] is the component Part of Npart components.  Data for each
// Part are contiguous in the file, so stride = 1.

  slab_stride[0] = 1;      // contiguous elements
  slab_count[0] = 1;       // one component per call
  slab_offset[0] = dim;    // component Part of Npart

  slab_stride[1] = 1;                   // contiguous elements
  slab_count[1] = NumberOfParticles;    // field dimensions
  slab_offset[1] = 0;                   // complete field, no offset

  file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);

  if ( dim == 0 )
  {
    file_id = H5Fcreate(PPos, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dset_id =  H5Dcreate(file_id, PPos, type_id, file_dsp_id, H5P_DEFAULT);

    attr_count = 1;
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
    attr_id = H5Acreate(dset_id, "NumberOfParticles", H5T_NATIVE_INT, attr_dsp_id, H5P_DEFAULT);
    h5_status = H5Awrite(attr_id, H5T_NATIVE_INT, &NumberOfParticles);
    h5_status = H5Aclose(attr_id);
    attr_id = H5Acreate(dset_id, "TotalParticleCount", H5T_NATIVE_INT, attr_dsp_id, H5P_DEFAULT);
    h5_status = H5Awrite(attr_id, H5T_NATIVE_INT, &TotalParticleCount);
    h5_status = H5Aclose(attr_id);
    h5_status = H5Sclose(attr_dsp_id);

  }
  else
  {
    file_id = H5Fopen(PPos, H5F_ACC_RDWR, H5P_DEFAULT);
    dset_id =  H5Dopen(file_id, PPos);
  }

  h5_status = H5Dwrite(dset_id, type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, outbuff);

  h5_status = H5Dclose(dset_id);
  h5_status = H5Sclose(mem_dsp_id);
  h5_status = H5Sclose(file_dsp_id);
  h5_status = H5Fclose(file_id);

  } // end of loop over dims


//  Read particle velocities under bit mask

  for (dim = 0; dim < rank; dim++)
  {

  ppc = 0; // output particle counter

  dims[0] = 3;
  dims[1] = ncpu;
  dims[2] = TotalParticleCount/ncpu;

  type_id = H5T_NATIVE_FLOAT;

//  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);

//  h5_status = H5Pset_fapl_mpio(file_acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_acc_template = H5P_DEFAULT;
  file_id = H5Fopen("ParticleVelocities", H5F_ACC_RDONLY, file_acc_template);

  dset_id = H5Dopen(file_id, "ParticleVelocities");

  mem_dsp_id = H5Screate_simple(1, &dims[2], NULL);

  mem_offset = 0;
  mem_count = dims[2];
  mem_stride = 1;

  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);

  file_dsp_id = H5Screate_simple(3, dims, NULL);

  slab_offset[0] = dim;   // x,y,z
  slab_stride[0] = 1;
  slab_count[0] = 1;

  slab_offset[1] = jcpu;
  slab_stride[1] = 1;
  slab_count[1] = 1;

  slab_offset[2] = 0;
  slab_stride[2] = 1;
  slab_count[2] = dims[2];

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);

//  xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);

//  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);

  xfer_prop_list = H5P_DEFAULT;
  h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, buff[0]);

//  h5_status = H5Pclose(xfer_prop_list);
//  h5_status = H5Pclose(file_acc_template);
  h5_status = H5Sclose(file_dsp_id);
  h5_status = H5Sclose(mem_dsp_id);
  h5_status = H5Dclose(dset_id);
  h5_status = H5Fclose(file_id);

  MPI_Barrier(MPI_COMM_WORLD);

  thisnode = jcpu;
  nextnode = (thisnode+1) % ncpu;
  prevnode = (thisnode-1+ncpu) % ncpu;

  stype = 30000+thisnode;
  rtype = 30000+prevnode;
  ltype = 30000+nextnode;

  iblk = - ( (ncpu - jcpu) % ncpu );

  for (k = 0; k < 2*ncpu; k++) // two passes
  {

  i = (k  ) % 2;
  j = (k+1) % 2;

  MPI_Barrier(MPI_COMM_WORLD);

  // Transmit left

  ier = MPI_Irecv( buff[j], dims[2], MPI_FLOAT, nextnode, ltype, MPI_COMM_WORLD, req1 );
  if (io_log) fprintf(log,"IRECV proc %d buffer[%d] err %d\n",thisnode,j,ier);

  ier = MPI_Isend( buff[i], dims[2], MPI_FLOAT, prevnode, stype, MPI_COMM_WORLD, req2 );
  if (io_log) fprintf(log,"SEND proc %d buffer[%d] err %d\n",thisnode,i,ier);

  // work on buff[i]

  if (io_log) fprintf(log,"K = %d, IBLK = %d\n",k,iblk);

  if ( iblk > -1 && iblk < ncpu )
  {

    ipc = (jcpu+k+ncpu) % ncpu;

    ipc = ipc * dims[2];
    jpc = ipc + dims[2] - 1;

    if (io_log) fprintf(log,"  ipc = %d to %d\n",ipc,jpc);

    for (ic = 0; ic < dims[2]; ic++)
    {
      MaskAddr = ipc;
      WordAddr = MaskAddr/BitsPerInt;
      BitAddr  = MaskAddr%BitsPerInt;
      TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
      if ( TestBit == 0 )
        XMask = 0;
      if ( TestBit != 0 )
      {
        XMask = 1;
        outbuff[ppc] = buff[i][ic];
        ppc++;
      }
      ipc++;
    }

  }

  iblk++;

  ier = MPI_Wait( req1, stat );

  if (io_log) fprintf(log,"WAIT proc %d err %d\n",thisnode,ier);

  } // end loop over rings

  fprintf(log,"L %8.4f  R %8.4f\n",GridLeft[jcpu][dim],GridRight[jcpu][dim]);

  if (io_log)
  {
    for (i = 0; i < NumberOfParticles; i++)
    {
      fprintf(log,"%8.4f",outbuff[i]);
      if ( ((i+1) % 16) == 0 )
        fprintf(log,"\n");
    }
  }

  Slab_Rank = 2;
  Slab_Dims[0] = 3;
  Slab_Dims[1] = NumberOfParticles;

/*
  component_rank_attr = 3;
  component_size_attr = NumberOfParticles;
  field_rank_attr = 1;
  field_dims_attr = NumberOfParticles;
*/

// Data in memory is considered 1D, stride 1, with zero offset

  mem_stride = 1;                   // contiguous elements
  mem_count = NumberOfParticles;    // number of elements in field
  mem_offset = 0;                   // zero offset in buffer

// 1D memory model

  mem_dsp_id = H5Screate_simple(1, &mem_count, NULL);

  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);

// Data in the file is (1+Rank)D with Npart components per grid point.
// Offset[0] is the component Part of Npart components.  Data for each
// Part are contiguous in the file, so stride = 1.

  slab_stride[0] = 1;      // contiguous elements
  slab_count[0] = 1;       // one component per call
  slab_offset[0] = dim;    // component Part of Npart

  slab_stride[1] = 1;                   // contiguous elements
  slab_count[1] = NumberOfParticles;    // field dimensions
  slab_offset[1] = 0;                   // complete field, no offset

  file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);

  if ( dim == 0 )
  {
    file_id = H5Fcreate(PVel, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dset_id =  H5Dcreate(file_id, PVel, type_id, file_dsp_id, H5P_DEFAULT);

    attr_count = 1;
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
    attr_id = H5Acreate(dset_id, "NumberOfParticles", H5T_NATIVE_INT, attr_dsp_id, H5P_DEFAULT);
    h5_status = H5Awrite(attr_id, H5T_NATIVE_INT, &NumberOfParticles);
    h5_status = H5Aclose(attr_id);
    attr_id = H5Acreate(dset_id, "TotalParticleCount", H5T_NATIVE_INT, attr_dsp_id, H5P_DEFAULT);
    h5_status = H5Awrite(attr_id, H5T_NATIVE_INT, &TotalParticleCount);
    h5_status = H5Aclose(attr_id);
    h5_status = H5Sclose(attr_dsp_id);

  }
  else
  {
    file_id = H5Fopen(PVel, H5F_ACC_RDWR, H5P_DEFAULT);
    dset_id =  H5Dopen(file_id, PVel);
  }

  h5_status = H5Dwrite(dset_id, type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, outbuff);

  h5_status = H5Dclose(dset_id);
  h5_status = H5Sclose(mem_dsp_id);
  h5_status = H5Sclose(file_dsp_id);
  h5_status = H5Fclose(file_id);

  } // end of loop over dims

// CHOP  } 

  fclose(log);

  MPI_Finalize();
}

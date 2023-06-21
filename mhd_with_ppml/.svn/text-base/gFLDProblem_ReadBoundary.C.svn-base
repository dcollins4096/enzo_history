/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Implicit Problem Class read external 
/  boundary routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: This routine reads the radiation energy boundary data from 
/           the specified file(s).  If the preprocessor directive 
/           PARALLEL_HDF5 is on, each processor reads its local portion 
/           of the boundary data from the same HDF restart file.  This 
/           allows restarts using a different number of processors than 
/           the run that wrote the boundary conditions to file, but 
/           requires the same top-grid dimensions as the previous run.
/           
/           Otherwise, each processor that owns part of the boundary 
/           reads its own restart file, containing only the pieces of 
/           the boundary local to that processor.  This mode requires 
/           the same processor number/layout/dimensions as the run 
/           that wrote the files.
/
************************************************************************/
#ifdef RAD_HYDRO
#ifdef USE_MPI
#include <mpi.h>
#endif 

#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "hdf4.h"
 
#include "gFLDProblem_preincludes.h"
#include "gFLDProblem.h"


// HDF5 function prototypes
 
#include "extern_hdf5.h"

// function prototypes

int ReadListOfInts(FILE *fptr, int N, int nums[]);
int ReadListOfFloats(FILE *fptr, int N, float floats[]);
 


int gFLDProblem::ReadBoundary(FILE *fptr)
{
 
//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::ReadBoundary routine\n");

  int size, dim, j, face;
  float32 *buffer;
  float *bdryface;

  char hdfname[MAX_LINE_LENGTH];
 
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
  hid_t       float_type_id, file_type_id;
  hid_t       file_dsp_id, mem_dsp_id;
 
  hsize_t     mem_stride[2], file_stride[2];
  hsize_t     mem_count[2],  file_count[2];
 
  hssize_t    mem_offset[2], file_offset[2];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  hsize_t     MemDims[2];

  const char *dname_value = "RadiationEnergyBoundary";

#ifdef PARALLEL_HDF5
  hsize_t     OutDims[2], FOffset[2];
#endif
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
 
  // Set constants for HDF float/file sizes
  int ii = sizeof(float32);
  switch(ii) {
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

 
  // Check that this simulation has same rank, dimensions and boundary 
  // condition types as restart
  int itest[MAX_DIMENSION];
  if (fscanf(fptr, "BoundaryRank = %"ISYM"\n", itest) != 1) {
    fprintf(stderr, "Error reading BoundaryRank.\n");
    return FAIL;
  }
  if (itest[0] != rank) {
    fprintf(stderr, "Incompatible restart BoundaryRank = %"ISYM"\n",
	    itest[0]);
    return FAIL;
  }
 
  fscanf(fptr, "RadiationDimensions =");
  if (ReadListOfInts(fptr, rank, itest) == FAIL) {
    fprintf(stderr, "Error reading RadiationDimensions.\n");
    return FAIL;
  }
  if ((itest[0] != GlobDims[0]) || 
      (itest[1] != GlobDims[1]) || 
      (itest[2] != GlobDims[2])) {
    fprintf(stderr, "Incompatible restart RadiationDimensions = %"ISYM" %"ISYM" %"ISYM"\n",
	    itest[0], itest[1], itest[2]);
    return FAIL;
  }
 
  fscanf(fptr, "RadiationBoundaryTypes =");
  if (ReadListOfInts(fptr, rank, itest) == FAIL) {
    fprintf(stderr, "Error reading RadiationBoundaryTypes.\n");
    return FAIL;
  }
  if ((itest[0] != BdryType[0][0]) || 
      (itest[1] != BdryType[1][0]) || 
      (itest[2] != BdryType[2][0])) {
    fprintf(stderr, "Incompatible restart RadiationBoundaryTypes = %"ISYM" %"ISYM" %"ISYM"\n",
	    itest[0], itest[1], itest[2]);
    return FAIL;
  }
 

  // check if any local boundaries are non-periodic
  bool OwnDirichletBoundary = false;
  for (dim=0; dim<rank; dim++) {
    if ((BdryType[dim][0] != 0) && (OnBdry[dim][0] || OnBdry[dim][1]))
      OwnDirichletBoundary = true;
  }

  // return if I don't own any isolating boundaries
  if (!OwnDirichletBoundary)
    return SUCCESS;
  
  // read HDF file name
  if (fscanf(fptr, "RadiationBoundaryFileName = %s\n", hdfname) != 1) {
    fprintf(stderr, "Error reading RadiationBoundaryFileName.\n");
    return FAIL;
  }
 
  
  // Set up IO log file, output general information
  if (io_log) {
    char *logname = new char[strlen(hdfname)+6];
    strcpy(logname, hdfname);
    strcat(logname, ".log2");
    log_fptr = fopen(logname, "a");
    delete[] logname;

    fprintf(log_fptr, "Radiation ReadBoundary start, proc %"ISYM"\n",
	    MyProcessorNumber);
    fprintf(log_fptr, "  P%"ISYM", BoundaryRank %"ISYM"\n", 
	    MyProcessorNumber, rank);
    fprintf(log_fptr, "  P%"ISYM", RadiationBoundaryFileName = %s\n", 
	    MyProcessorNumber, hdfname);
    fprintf(log_fptr, "  P%"ISYM", LocalDimension = %"ISYM" %"ISYM" %"ISYM"\n", 
	    MyProcessorNumber, LocDims[0], LocDims[1], LocDims[2]);
    fprintf(log_fptr, 
	    "  P%"ISYM", RadiationBoundaryTypes = %"ISYM" %"ISYM" %"ISYM"\n",
	    MyProcessorNumber, BdryType[0][0], BdryType[1][0], BdryType[2][0]); 
    fprintf(log_fptr, "  P%"ISYM", OwnBoundary = [%i,%i]x[%i,%i]x[%i,%i]\n",
	    MyProcessorNumber, 
	    Eint32(OnBdry[0][0]), Eint32(OnBdry[0][1]), 
	    Eint32(OnBdry[1][0]), Eint32(OnBdry[1][1]), 
	    Eint32(OnBdry[2][0]), Eint32(OnBdry[2][1]));
  }
  

#ifdef PARALLEL_HDF5

  /****** Parallel I/O of boundary conditions from the same HDF5 file *****/
  
  // Open HDF file
  strcat(hdfname, ".hdf");
  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Fopen with Name = %s\n", 
		      MyProcessorNumber, hdfname);
  file_id = H5Fopen(hdfname, H5F_ACC_RDONLY, H5P_DEFAULT);
  assert( file_id != h5_error );
  
  
  // Loop over dims, reading in isolating potential boundary values
  for (dim = 0; dim < rank; dim++)
    if (BdryType[dim][0] != 0) {
      for (face = 0; face < 2; face++)
	if (OnBdry[dim][face]) {
	  
	  if (io_log) fprintf(log_fptr, 
			      "  P%"ISYM", dim %"ISYM" : face %"ISYM"\n", 
			      MyProcessorNumber, dim, face);
	  
	  // Calculate size and dims of flux plane
	  // (ensure correct permutations for each dim)
	  if (dim==0) {
	    MemDims[0] = LocDims[1];
	    MemDims[1] = LocDims[2];
	    OutDims[0] = GlobDims[1];
	    OutDims[1] = GlobDims[2];
	    FOffset[0] = EdgeIndices[1][0];
	    FOffset[1] = EdgeIndices[2][0];
	    size  = MemDims[0]*MemDims[1];
	  } else if (dim==1) {
	    MemDims[0] = LocDims[2];
	    MemDims[1] = LocDims[0];
	    OutDims[0] = GlobDims[2];
	    OutDims[1] = GlobDims[0];
	    FOffset[0] = EdgeIndices[2][0];
	    FOffset[1] = EdgeIndices[0][0];
	    size  = MemDims[0]*MemDims[1];
	  } else {
	    MemDims[0] = LocDims[0];
	    MemDims[1] = LocDims[1];
	    OutDims[0] = GlobDims[0];
	    OutDims[1] = GlobDims[1];
	    FOffset[0] = EdgeIndices[0][0];
	    FOffset[1] = EdgeIndices[1][0];
	    size  = MemDims[0]*MemDims[1];
	  }
	  
	  
	  // set name for dim/face pair
	  char *nfile = new char[2];
	  char *dname = new char[strlen(dname_value)+6];
	  nfile[0] = '\0';
	  dname[0] = '\0';
	  sprintf(nfile,"%"ISYM,dim);
	  strcat(strcat(strcat(dname,dname_value),"."),nfile);
	  sprintf(nfile,"%"ISYM,face);
	  strcat(strcat(dname,"."),nfile);
	  
	  
	  // Create simple HDF5 memory dataspace
	  mem_dsp_id = H5Screate_simple((Eint32) 2, MemDims, NULL);
	  assert( mem_dsp_id != h5_error );
	  
	  // Create simple HDF5 file dataspace
	  file_dsp_id = H5Screate_simple((Eint32) 2, OutDims, NULL);
	  assert( file_dsp_id != h5_error );
	  
	  // Open HDF5 dataset for Boundary values on this dim/face
	  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Dopen w/ Name = %s\n", 
			      MyProcessorNumber, dname);
	  dset_id =  H5Dopen(file_id, dname);
	  assert( dset_id != h5_error );
	  
	  
	  // Allocate temporary space
	  buffer = new float32[size];
	  
	  //   Define HDF5 memory hyperslab
	  mem_offset[0] = 0;
	  mem_offset[1] = 0;
	  mem_stride[0] = 1;
	  mem_stride[1] = 1;
	  mem_count[0]  = MemDims[0];
	  mem_count[1]  = MemDims[1];
	  h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, 
					  mem_offset, mem_stride, 
					  mem_count, NULL);
	  assert( h5_status != h5_error );
	  
	  //   Define HDF5 file hyperslab
	  file_offset[0] = FOffset[0];
	  file_offset[1] = FOffset[1];
	  file_stride[0] = 1;
	  file_stride[1] = 1;
	  file_count[0]  = OutDims[0];
	  file_count[1]  = OutDims[1];
	  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, 
					  file_offset, file_stride, 
					  file_count, NULL);
	  assert( h5_status != h5_error );
	  
	  // Read Boundary Values from dataset
	  h5_status = H5Dread(dset_id, float_type_id, mem_dsp_id, 
			      file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
	  assert( h5_status != h5_error );
	  
	  // copy buffer to float temporary array, store
	  bdryface = new float[size];
	  for (j=0; j<size; j++) 
	    bdryface[j] = float(buffer[j]);
	  this->SetupBoundary(dim, face, 0, bdryface);
	  
	  // Close HDF5 boundary value dataset
	  h5_status = H5Dclose(dset_id);
	  assert( h5_status != h5_error );
	  
	  // Close HDF5 memory dataspace
	  h5_status = H5Sclose(mem_dsp_id);
	  assert( h5_status != h5_error );
	  
	  // Close HDF5 file dataspace
	  h5_status = H5Sclose(file_dsp_id);
	  assert( h5_status != h5_error );
	  
	  // Close HDF5 transfer property
	  h5_status = H5Pclose(xfer_prop_list);
	  assert( h5_status != h5_error );
	  
	  // Free allocated memory
	  delete[] bdryface;
	  delete[] buffer;
	  delete[] nfile;
	  delete[] dname;
	  
	}  // end of loop over owned faces
    }  // end of loop over isolating dims
  
  // Close HDF file
  h5_status = H5Fclose(file_id);
  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Fclose: %"ISYM"\n", 
		      MyProcessorNumber, h5_status);
  assert( h5_status != h5_error );
  
  
#else
  /****** I/O of boundary conditions from separate HDF5 files *****/ 
  
  
  // Determine this proc's HDF output file
  char *procstr = new char[7];
  sprintf(procstr,"%"ISYM,MyProcessorNumber);
  strcat(strcat(strcat(hdfname, "."),procstr), ".hdf");
  delete[] procstr;
  
  // Open HDF5 potential boundary file
  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Fopen with Name = %s\n", 
		      MyProcessorNumber, hdfname);
  file_id = H5Fopen(hdfname, H5F_ACC_RDONLY, H5P_DEFAULT);
  assert( file_id != h5_error );
  
  
  // Loop over dims, reading in isolating potential boundary values 
  for (dim = 0; dim < rank; dim++)
    if (BdryType[dim][0] != 0) {
      for (face = 0; face < 2; face++)
	if (OnBdry[dim][face]) { 
	  
	  if (io_log) fprintf(log_fptr, 
			      "  P%"ISYM", dim %"ISYM" : face %"ISYM"\n", 
			      MyProcessorNumber, dim, face);
	  
	  // calculate size and dims of flux plane
	  // (ensure correct permutations for each dim)
	  if (dim==0) {
	    MemDims[0] = LocDims[1];
	    MemDims[1] = LocDims[2];
	    size  = MemDims[0]*MemDims[1];
	  } else if (dim==1) {
	    MemDims[0] = LocDims[2];
	    MemDims[1] = LocDims[0];
	    size  = MemDims[0]*MemDims[1];
	  } else {
	    MemDims[0] = LocDims[0];
	    MemDims[1] = LocDims[1];
	    size  = MemDims[0]*MemDims[1];
	  }
	  
	  // set name for dim/face pair
	  char *nfile = new char[2];
	  char *dname = new char[strlen(dname_value)+6];
	  nfile[0] = '\0';
	  dname[0] = '\0';
	  sprintf(nfile,"%"ISYM,dim);
	  strcat(strcat(strcat(dname,dname_value),"."),nfile);
	  sprintf(nfile,"%"ISYM,face);
	  strcat(strcat(dname,"."),nfile);
	  
	  
	  // Create simple HDF5 memory dataspace
	  mem_dsp_id = H5Screate_simple((Eint32) 2, MemDims, NULL);
	  assert( mem_dsp_id != h5_error );
	  
	  // Create simple HDF5 file dataspace
	  file_dsp_id = H5Screate_simple((Eint32) 2, MemDims, NULL);
	  assert( file_dsp_id != h5_error );
	  
	  // Open HDF5 dataset for Boundary values on this dim/face
	  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Dopen with Name = %s\n", 
			      MyProcessorNumber, dname);
	  dset_id =  H5Dopen(file_id, dname);
	  assert( dset_id != h5_error );
	  
	  
	  // Allocate temporary space
	  buffer = new float32[size];

	  //  Define HDF5 memory hyperslab
	  mem_offset[0] = 0;
	  mem_offset[1] = 0;
	  mem_stride[0] = 1;
	  mem_stride[1] = 1;
	  mem_count[0]  = MemDims[0];
	  mem_count[1]  = MemDims[1];
	  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, 
					   mem_offset, mem_stride, 
					   mem_count, NULL);
	  assert( h5_status != h5_error );
	  
	  //  Define HDF5 file hyperslab
	  file_offset[0] = 0;
	  file_offset[1] = 0;
	  file_stride[0] = 1;
	  file_stride[1] = 1;
	  file_count[0]  = MemDims[0];
	  file_count[1]  = MemDims[1];
	  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, 
					  file_offset, file_stride, 
					  file_count, NULL);
	  assert( h5_status != h5_error );
	  
	  // read Boundary Values from dataset
	  h5_status = H5Dread(dset_id, float_type_id, mem_dsp_id, 
			      file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
	  assert( h5_status != h5_error );
	  
	  // copy buffer to float temporary array, store
	  bdryface = new float[size];
	  for (j = 0; j < size; j++)
	    bdryface[j] = float(buffer[j]);
	  this->SetupBoundary(dim, face, 0, bdryface);

	  // Close HDF5 boundary value dataset
	  h5_status = H5Dclose(dset_id);
	  assert( h5_status != h5_error );
	  
	  //  Close HDF5 memory dataspace
	  h5_status = H5Sclose(mem_dsp_id);
	  assert( h5_status != h5_error );
	  
	  //  Close HDF5 file dataspace
	  h5_status = H5Sclose(file_dsp_id);
	  assert( h5_status != h5_error );
      
	  // Free allocated memory
	  delete[] bdryface;
	  delete[] buffer;
	  delete[] nfile;
	  delete[] dname;

	}  // end of loop over owned faces
    }   // end of loop over isolating dims


  // Close HDF5 file
  h5_status = H5Fclose(file_id);
  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Fclose: %"ISYM"\n", 
		      MyProcessorNumber, h5_status);
  assert( h5_status != h5_error );
  
#endif
  
  // close log file, wait for all procs to catch up
  if (io_log) fclose(log_fptr);
  MPI_Barrier(MPI_COMM_WORLD);
  
  return SUCCESS;
}
#endif

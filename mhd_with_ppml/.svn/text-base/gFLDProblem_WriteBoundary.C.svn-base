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
/  Gray Flux-Limited Diffusion Implicit Problem Class write external 
/  boundary routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: This routine writes the radiation energy boundary data 
/           to the specified file(s).  If the preprocessor directive 
/           PARALLEL_HDF5 is on, each processor writes its portion of 
/           the boundary data to the same HDF restart file, as given 
/           in the function argument list.  
/
/           Otherwise, each processor that owns part of the boundary 
/           writes its own restart file (with name given by the 
/           supplied restart file name, appended by ".<proc number>").  
/           These files contain only the pieces of the boundary local 
/           to the specific processor.
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
 
void WriteListOfInts(FILE *fptr, int N, int nums[]);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
int FileMover(char *file_to_move);


int gFLDProblem::WriteBoundary(FILE *fptr, char *hdfname)
{

//   if (debug)
//     fprintf(stdout,"Entering gFLDProblem::WriteBoundary routine\n");

  int dim, face, i, j, size;
  float32 *buffer;
  
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
  hid_t       float_type_id, file_type_id;
  hid_t       file_dsp_id, mem_dsp_id;

  hssize_t    mem_offset[2], file_offset[2];

  hsize_t     mem_stride[2], file_stride[2];
  hsize_t     mem_count[2],  file_count[2];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;

  hsize_t     MemDims[2];
 
  const char *dname_value = "RadiationEnergyBoundary";

#ifdef PARALLEL_HDF5
  hid_t       file_access_template;
  hid_t       xfer_prop_list;
  hsize_t     OutDims[2], FOffset[2];
#endif
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  // check if any local boundaries are non-periodic
  bool OwnDirichletBoundary = false;
  for (dim=0; dim<rank; dim++) {
    if ((BdryType[dim][0] != 0) && (OnBdry[dim][0] || OnBdry[dim][1]))
      OwnDirichletBoundary = true;
  }

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
 
  // Save general class data (root proc only)
  if (MyProcessorNumber == 0) {
    fprintf(fptr, "BoundaryRank = %"ISYM"\n", rank);
    
    fprintf(fptr, "RadiationDimensions = ");
    WriteListOfInts(fptr, rank, GlobDims);
    
    fprintf(fptr, "RadiationBoundaryTypes = %"ISYM" %"ISYM" %"ISYM"\n", 
	    BdryType[0][0], BdryType[1][0], BdryType[2][0]);
    
    // Write out information about the BoundaryValues
    fprintf(fptr, "RadiationBoundaryFileName = %s\n", hdfname);
  }
  
  // return if I don't own any isolating boundaries
  if (!OwnDirichletBoundary)
    return SUCCESS;
  
  // Set up IO log file, output general information
  if (io_log) {
    char *logname = new char[strlen(hdfname)+5];
    strcpy(logname, hdfname);
    strcat(logname, ".log");
    log_fptr = fopen(logname, "a");
    delete[] logname;
    
    fprintf(log_fptr, "Radiation WriteBoundary start, proc %"ISYM"\n",
	    MyProcessorNumber);
    fprintf(log_fptr, "  P%"ISYM", BoundaryRank %"ISYM"\n", 
	    MyProcessorNumber, rank);
    fprintf(log_fptr, "  P%"ISYM", RadiationBoundaryFileName = %s\n", 
	    MyProcessorNumber, hdfname);
    fprintf(log_fptr, "  P%"ISYM", LocalDimension = %"ISYM" %"ISYM" %"ISYM"\n",
	    MyProcessorNumber, LocDims[0], LocDims[1], LocDims[2]); 
    fprintf(log_fptr, 
	    "  P%"ISYM", RadiationBoundaryTypes = %"ISYM" %"ISYM" %"ISYM"\n",
	    MyProcessorNumber, BdryType[0][0], 
	    BdryType[1][0], BdryType[2][0]); 
    fprintf(log_fptr, "  P%"ISYM", OwnBoundary = [%i,%i]x[%i,%i]x[%i,%i]\n",
	    MyProcessorNumber, 
	    Eint32(OnBdry[0][0]), Eint32(OnBdry[0][1]), 
	    Eint32(OnBdry[1][0]), Eint32(OnBdry[1][1]), 
	    Eint32(OnBdry[2][0]), Eint32(OnBdry[2][1]));
  }


#ifdef PARALLEL_HDF5

  /****** Parallel I/O of boundary conditions to the same HDF5 file *****/
 
  // Set file access template for parallel HDF5
  file_access_template = H5Pcreate(H5P_FILE_ACCESS);
  assert( file_access_template != h5_error );

  // Store MPI communicator information to the file access prop. list
  h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, 
			       MPI_INFO_NULL);
  assert( h5_status != h5_error );

  // Create HDF file
  char *myhdf = new char[strlen(hdfname)+5];
  strcpy(myhdf, hdfname);
  strcat(myhdf, ".hdf");
  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Fcreate with Name = %s\n", 
		      MyProcessorNumber, myhdf);
  file_id = H5Fcreate(myhdf, H5F_ACC_TRUNC, H5P_DEFAULT, 
		      file_access_template);
  assert( file_id != h5_error );
  delete[] myhdf;

  // Terminate access to property list
  h5_status = H5Pclose(file_access_template);
  assert( h5_status != h5_error );


  // Loop over dims, writing out isolating potential boundary values
  for (dim=0; dim<rank; dim++)
    if (BdryType[dim][0] != 0) {
      for (face=0; face<2; face++)
	if (OnBdry[dim][face] && (EBdryVals[dim][face] != NULL)) {

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
	  delete[] nfile;


	  // Create simple HDF5 memory dataspace
	  mem_dsp_id = H5Screate_simple((Eint32) 2, MemDims, NULL);
	  assert( mem_dsp_id != h5_error );
 
	  // Create simple HDF5 file dataspace
	  file_dsp_id = H5Screate_simple((Eint32) 2, OutDims, NULL);
	  assert( file_dsp_id != h5_error );
 
	  // Create HDF5 dataset for Boundary values on this dim/face
	  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Dcreate w/ Name = %s\n", 
			      MyProcessorNumber, dname);
	  dset_id =  H5Dcreate(file_id, dname, file_type_id, 
			       file_dsp_id, H5P_DEFAULT);
	  assert( dset_id != h5_error );
	  delete[] dname;
 

	  // Allocate temporary space
	  buffer = new float32[size];

	  // fill temporary buffer
	  for (j = 0; j < int(MemDims[1]); j++)
	    for (i = 0; i < int(MemDims[0]); i++) 
	      buffer[j*MemDims[0] + i] = 
		float32(EBdryVals[dim][face][j*MemDims[0]+i]);
	  

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

	  // Create a new property for raw data transfer
	  xfer_prop_list = H5Pcreate(H5P_DATASET_XFER);
	  assert( xfer_prop_list != h5_error );

	  // Set data transfer mode to collective I/O access
	  h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
	  assert( h5_status != h5_error );
	  
	  // Write Boundary Values to dataset
	  h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, 
			       file_dsp_id,  xfer_prop_list, (VOIDP) buffer);
	  assert( h5_status != h5_error );
	  delete[] buffer;
	  
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

	}  // end of loop over owned faces
    }  // end of loop over isolating dims
  
  // Close HDF file
  h5_status = H5Fclose(file_id);
  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Fclose: %"ISYM"\n", 
		      MyProcessorNumber, h5_status);
  assert( h5_status != h5_error );
  

#else
  /****** I/O of boundary conditions to separate HDF5 files *****/ 

  // Determine this proc's HDF output file
  char *myhdf   = new char[strlen(hdfname)+11];
  char *procstr = new char[7];
  strcpy(myhdf, hdfname);
  sprintf(procstr,"%"ISYM,MyProcessorNumber);
  strcat(strcat(strcat(myhdf, "."),procstr), ".hdf");
  delete[] procstr;

  // Create HDF5 potential boundary file
  if (io_log) fprintf(log_fptr, "  P%"ISYM", H5Fcreate with Name = %s\n", 
		      MyProcessorNumber, myhdf);
  file_id = H5Fcreate(myhdf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  assert( file_id != h5_error );
  delete[] myhdf;

  // Loop over dims, writing out isolating potential boundary values
  for (dim = 0; dim < rank; dim++)
    if (BdryType[dim][0] != 0) {
      for (face = 0; face < 2; face++)
	if (OnBdry[dim][face] && (EBdryVals[dim][face] != NULL)) {

	  if (io_log) fprintf(log_fptr, 
			      "  P%"ISYM", dim %"ISYM" : face %"ISYM"\n", 
			      MyProcessorNumber, dim, face);
	  
	  // Calculate size and dims of flux plane
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
	  delete[] nfile;


	  // Create simple HDF5 memory dataspace
	  mem_dsp_id = H5Screate_simple((Eint32) 2, MemDims, NULL);
	  assert( mem_dsp_id != h5_error );
 
	  // Create simple HDF5 file dataspace
	  file_dsp_id = H5Screate_simple((Eint32) 2, MemDims, NULL);
	  assert( file_dsp_id != h5_error );
 
	  // Create HDF5 dataset for Boundary values on this dim/face
	  if (io_log) fprintf(log_fptr,
			      "  P%"ISYM", H5Dcreate w/ Name = %s\n", 
			      MyProcessorNumber, dname);
	  dset_id =  H5Dcreate(file_id, dname, file_type_id, 
			       file_dsp_id, H5P_DEFAULT);
	  assert( dset_id != h5_error );
	  delete[] dname;
 

	  // Allocate temporary space
	  buffer = new float32[size];

	  // fill temporary buffer
	  for (j = 0; j < int(MemDims[1]); j++)
	    for (i = 0; i < int(MemDims[0]); i++) 
	      buffer[j*MemDims[0] + i] = 
		float32(EBdryVals[dim][face][j*MemDims[0]+i]);
	  

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

	  // Write Boundary Values to dataset
	  h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, 
			       file_dsp_id, H5P_DEFAULT, (VOIDP) buffer);
	  assert( h5_status != h5_error );
	  delete[] buffer;
	  
	  // Close HDF5 boundary value dataset
	  h5_status = H5Dclose(dset_id);
	  assert( h5_status != h5_error );
      
	  // Close HDF5 memory dataspace
	  h5_status = H5Sclose(mem_dsp_id);
	  assert( h5_status != h5_error );
      
	  // Close HDF5 file dataspace
	  h5_status = H5Sclose(file_dsp_id);
	  assert( h5_status != h5_error );

	}  // end of loop over owned faces
    }  // end of loop over isolating dims
  

  // Close HDF file
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

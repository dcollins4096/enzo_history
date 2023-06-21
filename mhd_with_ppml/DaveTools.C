#include <assert.h>
#include "hdf5.h"
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "extern_hdf5.h"

#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

//returns the wall time.
void wall_time (char * string)
{

  // get current time in seconds (specifically the number of 
  // seconds since 00:00:00 UTC January 1 1970)              
  
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv,&tz);
  fprintf(stderr,"%d TIME %s %f \n", MyProcessorNumber, string, tv.tv_sec + 1e-6*tv.tv_usec);

}

//Writes an HDF5 cube.  Quick and dirty.  
// WriteCube(array, [nx,ny,nz], "ID string", dNum, gNum)
// prints out ID string to std err, along with file name.
// Filename is data111(dNum).grid(gNum)
// datasets are all called "monkey"
// Flow:
// 1.) create filename, print message
// 2.) define size of float 
// 3.) create file
// 3.5) invert array dimensions
// 4.) create dataspace, set
// 5.) close file,space,set.
void WriteCube_dcc(float * array, int Dims[], char* string, int dNum, int gNum){

  
  hid_t       file_id, dataset_id, dataspace_id, float_type_id, file_type_id;
  herr_t      status, h5_status, h5_error = -1;
  int FieldRankOut = 3;
  hsize_t     DimsInv[FieldRankOut];
  
  char filename[20];
  
  sprintf(filename, "data111%4.4d.grid%4.4d",dNum,gNum);
  fprintf(stderr,"GPFS WriteCube: %s %s [%d,%d,%d]\n", string, filename, Dims[0],Dims[1],Dims[2]);
  
#define floatdcc double  
  int ii = sizeof(floatdcc);
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
  
  
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  for(int moo=0;moo<FieldRankOut;moo++)
    DimsInv[moo]=Dims[FieldRankOut-moo-1];
  
  
  //Create Dataspace 
  dataspace_id=H5Screate_simple((Eint32) FieldRankOut, DimsInv, NULL);
  
  //create set
  //                       duh, name,      datatype,  shape of data, Something I dont get
  dataset_id = H5Dcreate(file_id, "monkey", file_type_id, dataspace_id, H5P_DEFAULT);
  
  //Write the Data Set
  //                (set, memory type, mem. space, file space, transfer info, actual data)
  fprintf(stderr,"Writing set\n");
  status = H5Dwrite(dataset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		    (VOIDP) array);
  
  
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  
  
}


// If 'flag' matches any element int MidWayDumpList, return true.
int MidWayDumpCheck( int flag ){
  int DumpHere = FALSE;
  for( int chk=0; chk< N_MidWayDumps; chk++){
    if( MidWayDumpList[ chk ] == flag )
      DumpHere = TRUE;
  }
  return DumpHere;
}

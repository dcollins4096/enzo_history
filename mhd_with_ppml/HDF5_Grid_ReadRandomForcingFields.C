/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/
/  PURPOSE: At restart reads initial velocity fields from the first output
/           data file(s) (e.g. data_0000.grid000?) and stores them as
/           RandomForcingField[].
/
************************************************************************/
 
#define USE_HDF5
 
#ifdef USE_HDF5
 
//  Input a grid from file pointer fpt
 
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
 
#include "hdf4.h"
 
//#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
 
// HDF5 function prototypes
 
#include "extern_hdf5.h"
 
// function prototypes
#ifdef PPML
int FindField(int field, int farray[], int numfields);
#include "PPML.h"
#endif //PPML 
int ReadListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
 
 
 
int grid::ReadRandomForcingFieldsHDF5(FILE *fptr)
{
 
  int i, j, k, dim, field, size, active_size;
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
 
  int ActiveDim[MAX_DIMENSION];
 
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_dsp_id;
  hid_t       num_type;
 
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[MAX_DIMENSION];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int         num_size;
 
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
      break;
 
    case 8:
      float_type_id = HDF5_R8;
      break;
 
    default:
      float_type_id = HDF5_R4;
 
  }
 
  int jj = sizeof(FLOAT);
 
  switch(jj)
  {
 
    case 4:
      FLOAT_type_id = HDF5_R4;
      break;
 
    case 8:
      FLOAT_type_id = HDF5_R8;
      break;
 
    case 16:
      FLOAT_type_id = HDF5_R16;
      break;
 
    default:
      printf("INCORRECT FLOAT DEFINITION\n");
 
  }
 
 
  /* Store the file position indicator. */
 
  fpos_t position;
  if (fgetpos(fptr, &position) != 0)
    WARNING_MESSAGE;
 
  /* Assume general grid class data are known and
     read velocity fields only (as RandomForcingFields).
     Do not modify the existing current BaryonFields and grid class data. */
 
  if (NumberOfBaryonFields <= GridRank + 1) {
    fprintf(stderr, "Error: No Baryon Fields => Nothing to Force. Right?.\n");
    ERROR_MESSAGE;
  }
 
  /* Read the filename where the current baryon fields are; assume that
     initial fields were sitting in a file with the same name but different
     number; change the current number into '0000' and, thus, prepare the
     required name. */
 
  if (fsetpos(fptr, &BaryonFileNamePosition) != 0)
    ERROR_MESSAGE;
  if (fscanf(fptr, "BaryonFileName = %s\n", name) != 1) {
    fprintf(stderr, "Error reading BaryonFileName.\n");
    ERROR_MESSAGE;
  }
 
  /* this wont work if file prefix contains anoter ".grid" and/or
     if restart is invoked for cycle # > 9999 */
 
  char * numberEnd = NULL;
  if ( (numberEnd = strstr(name, ".grid")) == NULL )
    ERROR_MESSAGE;
  *(numberEnd - 1) = '0';
  *(numberEnd - 2) = '0';
  *(numberEnd - 3) = '0';
  *(numberEnd - 4) = '0';
 
  /* check the name. */
 
  if (debug)
    printf("ReadRandomForcingFields from: %s\n", name);
  if ( strstr(name, "0000.grid") == NULL )
    ERROR_MESSAGE;
 
  /* read fields */
 
  if (MyProcessorNumber == ProcessorNumber) {
    char *logname = new char[strlen(name)+9];
    strcpy(logname, name);
    strcat(logname, ".hdf.log");
    if (io_log) log_fptr = fopen(logname, "a");
    fprintf(stderr,"H5Fopen with Name %s\n",name);
    file_id = H5Fopen(name,  H5F_ACC_RDONLY, H5P_DEFAULT);
    if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
    assert( file_id != h5_error );
  }
 
  if (MyProcessorNumber == ProcessorNumber) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                         Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "HDF5GRRFF: Error in IdentifyPhysicalQuantities.\n");
      return FAIL;
    }
    int vel = Vel1Num;
    printf("RandomForcing: Fields %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" \n",
           DensNum, TENum, Vel1Num, Vel2Num, Vel3Num);
 
 
    /* fill in ActiveDim for dims up to 3d */
 
    for (dim = 0; dim < 3; dim++)
      ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
 
    /* check dimensions of HDF file against this grid
       (note: we don't bother to check the coordinate arrays)  */
 
    size = 1;
    active_size = 1;
 
    for (dim = 0; dim < GridRank; dim++) {
      size *= GridDimension[dim];
      active_size *= ActiveDim[dim];
    }
 
//  CAUTION - are the coordinates reversed?
 
    for (dim = 0; dim < GridRank; dim++) {
      OutDims[GridRank-dim-1] = ActiveDim[dim];
      if (io_log) fprintf(log_fptr, "Outdims %"ISYM"\n", (int) OutDims[GridRank-dim-1]);
    }
 
    /* allocate temporary space */
 
    float32 *temp = new float32[active_size];
 
    /* skip fields written before the first velocity field;
       dirty trick, sorry; this is done only once at a restart. */
 
    for (i = 0; i < vel; i++) {
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error );
 
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", DataLabel[i]);
 
      dset_id =  H5Dopen(file_id, DataLabel[i]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
      //      JBPERF_COUNT_READ(dset_id, float_type_id, H5S_ALL);
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
    }
 
    /* loop over fields, reading each one */
 
    for (dim = 0; dim < GridRank; dim++) {
 
      /* get data into temporary array */
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error );
 
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", DataLabel[vel+dim]);
 
      dset_id =  H5Dopen(file_id, DataLabel[vel+dim]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
      //JBPERF_COUNT_READ(dset_id, float_type_id, H5S_ALL);
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
      /* copy velocity from active region into whole grid forcing field;
         put zeroes into ghost zones. */
 
      if (RandomForcingField[dim] == NULL)
        RandomForcingField[dim] = new float[size];
      for (i = 0; i < size; i++)
        RandomForcingField[dim][i] = 0;
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    RandomForcingField[dim][i + j*GridDimension[0] +
				      k*GridDimension[0]*GridDimension[1]] =
	      float(temp[(i-GridStartIndex[0])                         +
	                 (j-GridStartIndex[1])*ActiveDim[0]            +
	                 (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ]);
 
    } // end: loop over fields

#ifdef PPML
    //Need to read the rest of the driving fields.
    if( RandomForcingNumberOfFields == 0 ){
      fprintf(stderr," Severe Error: RandomForcingNumberOfFields == 0 (or undefined.)\n");
      fprintf(stderr," Should be set to 3 for PPM runs, many for PPM-L runs.\n");
      return FAIL;
    }
    if( HydroMethod == PPM_Local && RandomForcingNumberOfFields > GridRank ){
      PPML_InterfacePointerBundle Face( this );
      int RandomCounter = GridRank;
      int LabelIndex[6] = {FindField( Face_X_L_VX, FieldType, NumberOfBaryonFields ),
			   FindField( Face_X_R_VX, FieldType, NumberOfBaryonFields ),
			   FindField( Face_Y_L_VX, FieldType, NumberOfBaryonFields ),
			   FindField( Face_Y_R_VX, FieldType, NumberOfBaryonFields ),
			   FindField( Face_Z_L_VX, FieldType, NumberOfBaryonFields ),
			   FindField( Face_Z_R_VX, FieldType, NumberOfBaryonFields )};
      for(int face_part=0; face_part < PPML_NFaces; face_part++)
	for (dim = 0; dim < GridRank; dim++) {
	  
	  file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
	  assert( file_dsp_id != h5_error );
	  dset_id =  H5Dopen(file_id, DataLabel[ LabelIndex[face_part] + dim ]);
	  fprintf(stderr,"Data Label: %s\n",DataLabel[ LabelIndex[face_part] + dim ]);
	  assert( dset_id != h5_error );
	  h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
	  assert( h5_status != h5_error );
	  
	  h5_status = H5Sclose(file_dsp_id);
	  h5_status = H5Dclose(dset_id);
	  
	  if (RandomForcingField[RandomCounter] == NULL)
	    RandomForcingField[RandomCounter] = new float[size];
	  for (i = 0; i < size; i++)
	    RandomForcingField[RandomCounter][i] = 0;
	  
	  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
		RandomForcingField[RandomCounter][i + j*GridDimension[0] +
						    k*GridDimension[0]*GridDimension[1]] =
		  float(temp[(i-GridStartIndex[0])                         +
			     (j-GridStartIndex[1])*ActiveDim[0]            +
			     (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ]);
	  RandomCounter++;
	} // end: loop over fields
      }//hydro method == PPML
    //<dbg>
    
#endif //PPML
    delete [] temp;
 
  }  // end: if (MyProcessorNumber == ProcessorNumber)
 
  /* Close file. */
 
  if ( (MyProcessorNumber == ProcessorNumber) &&
       (NumberOfBaryonFields > 0) )
  {
     h5_status = H5Fclose(file_id);
       if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
       assert( h5_status != h5_error );
  }
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    if (io_log) fclose(log_fptr);
  }
 
  /* Set the file position indicator. */
 
  if (fsetpos(fptr, &position) != 0)
    WARNING_MESSAGE;
 
  return SUCCESS;
 
}
 
#endif
 

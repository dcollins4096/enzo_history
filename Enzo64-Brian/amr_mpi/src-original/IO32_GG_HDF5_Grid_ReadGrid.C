/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Alexei Kritsuk, Jan 2004   a trick for RandomForcing //AK
/
/  PURPOSE:
/
************************************************************************/
 
//  Input a grid from file pointer fpt
 
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
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
 
// HDF5 function prototypes
 
#include "extern_hdf5.h"
 
// function prototypes
 
int ReadListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
 
// extern int ParticleTypeInFile; // declared and set in ReadParameterFile
static int GridReadDataGridCounter = 0;
 
 
int grid::ReadGrid(FILE *fptr, int GridID)
{
 
  int i, j, k, dim, field, size, active_size;
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
  char logname[MAX_LINE_LENGTH];
  char procfilename[MAX_LINE_LENGTH];
 
  char id[10];
  char pid[5];
  char gpid[5];
 
  int ActiveDim[MAX_DIMENSION];
 
  FILE *log_fptr;
 
  hid_t       file_id, group_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_type_id, FILE_type_id;
  hid_t       file_dsp_id;
  hid_t       num_type;
 
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[MAX_DIMENSION];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int         num_size;
 
  char *ParticlePositionLabel[] =
     {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleVelocityLabel[] =
     {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
                                    "metallicity_fraction", "alpha_fraction"};
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
  /* make sure quantities defined at least for 3d */
 
  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }
 
  /* Read general grid class data */
 
  if (fscanf(fptr, "GridRank = %"ISYM"\n", &GridRank) != 1) {
    fprintf(stderr, "Error reading GridRank.\n");
    return FAIL;
  }
 
  if (fscanf(fptr, "GridDimension = ") != 0) {
    fprintf(stderr, "Error reading GridDimension(0).\n");
    return FAIL;
  }
 
  if (ReadListOfInts(fptr, GridRank, GridDimension) == FAIL) {
    fprintf(stderr, "Error reading GridDimension(1).\n");
    return FAIL;
  }
 
  fscanf(fptr, "GridStartIndex = ");
 
  if (ReadListOfInts(fptr, GridRank, GridStartIndex) == FAIL) {
    fprintf(stderr, "Error reading GridStartIndex.\n");
    return FAIL;
  }
 
  fscanf(fptr, "GridEndIndex = ");
 
  if (ReadListOfInts(fptr, GridRank, GridEndIndex) == FAIL) {
    fprintf(stderr, "Error reading GridEndIndex.\n");
    return FAIL;
  }
 
  fscanf(fptr, "GridLeftEdge = ");
 
  if (ReadListOfFloats(fptr, GridRank, GridLeftEdge) == FAIL) {
    fprintf(stderr, "Error reading GridLeftEdge.\n");
    return FAIL;
  }
 
  fscanf(fptr, "GridRightEdge = ");
 
  if (ReadListOfFloats(fptr, GridRank, GridRightEdge) == FAIL) {
    fprintf(stderr, "Error reading GridRightEdge.\n");
    return FAIL;
  }
 
  if (fscanf(fptr, "Time = %"PSYM"\n", &Time) != 1) {
    fprintf(stderr, "Error reading Time.\n");
    return FAIL;
  }
 
  if (fscanf(fptr, "SubgridsAreStatic = %"ISYM"\n", &SubgridsAreStatic) != 1) {
    fprintf(stderr, "Error reading SubgridsAreStatic.\n");
    return FAIL;
  }
 
  /* Compute Flux quantities */
 
  this->PrepareGridDerivedQuantities();
 
  /* Read baryon field quantities. */
 
  if (fscanf(fptr, "NumberOfBaryonFields = %"ISYM"\n",
	     &NumberOfBaryonFields) != 1) {
    fprintf(stderr, "Error reading NumberOfBaryonFields.\n");
    return FAIL;
  }
 
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
      FILE_type_id = HDF5_FILE_R4;
      break;
 
    case 8:
      FLOAT_type_id = HDF5_R8;
      FILE_type_id = HDF5_FILE_R8;
      break;
 
    case 16:
      FLOAT_type_id = HDF5_R16;
      FILE_type_id = H5Tcopy(HDF5_FILE_B8);
                     H5Tset_size(FILE_type_id,16);
      break;
 
    default:
      printf("INCORRECT FLOAT DEFINITION\n");
 
  }
 
  sprintf(id, "%8.8"ISYM, GridID);
 
  sprintf(pid, "%4.4"ISYM, MyProcessorNumber);
 
  sprintf(gpid, "%4.4"ISYM, ProcessorNumber);
 
  strcpy(name, "/Grid");
  strcat(name, id);
 
  if (NumberOfBaryonFields > 0) {
 
    fscanf(fptr, "FieldType = ");
 
    if (ReadListOfInts(fptr, NumberOfBaryonFields, FieldType) == FAIL) {
      fprintf(stderr, "Error reading FieldType.\n");
      return FAIL;
    }
 
    fgetpos(fptr, &BaryonFileNamePosition); //AK
 
    if (fscanf(fptr, "BaryonFileName = %s\n", procfilename) != 1) {
      fprintf(stderr, "Error reading BaryonFileName.\n");
      return FAIL;
    }
 
    fscanf(fptr, "CourantSafetyNumber    = %"FSYM"\n", &CourantSafetyNumber);
    fscanf(fptr, "PPMFlatteningParameter = %"ISYM"\n", &PPMFlatteningParameter);
    fscanf(fptr, "PPMDiffusionParameter  = %"ISYM"\n", &PPMDiffusionParameter);
    fscanf(fptr, "PPMSteepeningParameter = %"ISYM"\n", &PPMSteepeningParameter);
 
    if (MyProcessorNumber == ProcessorNumber)
    {
      strcpy(logname, procfilename);
      strcat(logname, ".in_log");
      if (io_log) log_fptr = fopen(logname, "a");
    }
 
    if (MyProcessorNumber == ProcessorNumber)
    {
      if (io_log) fprintf(log_fptr, "H5Fopen with Name %s\n", procfilename);
 
      file_id = H5Fopen(procfilename,  H5F_ACC_RDONLY, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
        assert( file_id != h5_error );
 
      if (io_log) fprintf(log_fptr, "H5Gopen with Name %s\n", name);
 
      group_id = H5Gopen(file_id, name);
        assert( group_id != h5_error );
    }
 
    if (MyProcessorNumber == ProcessorNumber) {
 
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
 
    /* loop over fields, reading each one */
 
    for (field = 0; field < NumberOfBaryonFields; field++) {
 
      /* get data into temporary array */
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        assert( file_dsp_id != h5_error );
 
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", DataLabel[field]);
 
      dset_id =  H5Dopen(group_id, DataLabel[field]);
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        assert( dset_id != h5_error );
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      /* copy active region into whole grid */
 
      BaryonField[field] = new float[size];
 
      for (i = 0; i < size; i++)
	BaryonField[field][i] = 0;
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    BaryonField[field][i + j*GridDimension[0] +
		                   k*GridDimension[0]*GridDimension[1]] =
	      float(temp[(i-GridStartIndex[0])                         +
	                 (j-GridStartIndex[1])*ActiveDim[0]            +
	                 (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ]);
 
    } // end: loop over fields
 
    delete [] temp;
 
  }  // end: if (MyProcessorNumber == ProcessorNumber)
  }  // end read baryon fields
 
 
 
 
 
  /* 3) Read particle info */
 
  if (fscanf(fptr, "NumberOfParticles = %"ISYM"\n", &NumberOfParticles) != 1) {
    fprintf(stderr, "error reading NumberOfParticles.\n");
    return FAIL;
  }
 
 
  if (NumberOfParticles > 0) {
 
    /* Read particle file name. */
 
    if (fscanf(fptr, "ParticleFileName = %s\n", name) != 1) {
      fprintf(stderr, "Error reading ParticleFileName.\n");
      return FAIL;
    }
 
    if (MyProcessorNumber == ProcessorNumber) {
 
    /* Open file if not already done (note: particle name must = grid name). */
 
    if (NumberOfBaryonFields == 0) {
 
      if (MyProcessorNumber == ProcessorNumber)
      {
        strcpy(logname, procfilename);
        strcat(logname, ".in_log");
        if (io_log) log_fptr = fopen(logname, "a");
      }
 
      if (io_log) fprintf(log_fptr, "H5Fopen with Name %s\n", procfilename);
 
      file_id = H5Fopen(procfilename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
        assert( file_id != h5_error );
 
      if (io_log) fprintf(log_fptr, "H5Gopen with Name %s\n", name);
 
      group_id = H5Gopen(file_id, name);
        assert( group_id != h5_error );
 
    } // end: if (NumberOfBaryonFields == 0)
 
    /* Allocate room for particles. */
 
    this->AllocateNewParticles(NumberOfParticles);
 
    TempIntArray[0] = NumberOfParticles;
 
    /* Create a temporary buffer (32 bit or twice the size for 64). */
 
    float32 *temp = NULL;
 
    jj = sizeof(FLOAT);
 
    switch(jj)
    {
 
      case 4:
        temp = new float32[NumberOfParticles];
        break;
 
      case 8:
        temp = new float32[NumberOfParticles*2];
        break;
 
      case 16:
        temp = new float32[NumberOfParticles*4];
        break;
 
      default:
        printf("INCORRECT FLOAT DEFINITION\n");
 
    }
 
    if (temp == NULL)
      temp = new float32[NumberOfParticles];
 
    /* Read ParticlePosition (use temporary buffer). */
 
    for (dim = 0; dim < GridRank; dim++) {
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        assert( file_dsp_id != h5_error );
 
      if (io_log) fprintf(log_fptr,"H5Dopen with Name = %s\n", ParticlePositionLabel[dim]);
 
      dset_id =  H5Dopen(group_id, ParticlePositionLabel[dim]);
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        assert( dset_id != h5_error );
 
      num_type = H5Dget_type(dset_id);
      num_size = H5Tget_size(num_type);
 
      if (sizeof(FLOAT) == 16)
      {
 
//                                 NOTE: for 128bits this must be FILE_type_id and NOT FLOAT_type_id!
      h5_status = H5Dread(dset_id, FILE_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticlePosition[dim]);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      }
      else
      {
 
      h5_status = H5Dread(dset_id, FLOAT_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticlePosition[dim]);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      }
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
    }
 
 
    /* Read ParticleVelocity. */
 
    for (dim = 0; dim < GridRank; dim++) {
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        assert( file_dsp_id != h5_error );
 
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", ParticleVelocityLabel[dim]);
 
      dset_id =  H5Dopen(group_id, ParticleVelocityLabel[dim]);
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        assert( dset_id != h5_error );
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      for (i = 0; i < NumberOfParticles; i++)
	ParticleVelocity[dim][i] = float(temp[i]);
    }
 
 
    /* Read ParticleMass into temporary buffer and Copy to ParticleMass. */
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error );
 
    if (io_log) fprintf(log_fptr,"H5Dopen with Name = particle_mass\n");
 
    dset_id =  H5Dopen(group_id, "particle_mass");
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );
 
    h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    for (i = 0; i < NumberOfParticles; i++)
      ParticleMass[i] = float(temp[i]);
 
    /* Read ParticleNumber into temporary buffer and Copy to ParticleNumber. */
 
    int *tempint = new int[NumberOfParticles];
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      assert( file_dsp_id != h5_error);
 
    if (io_log) fprintf(log_fptr,"H5Dopen  with Name = particle_index\n");
 
    dset_id =  H5Dopen(group_id, "particle_index");
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );
 
    h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempint);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    for (i = 0; i < NumberOfParticles; i++)
      ParticleNumber[i] = tempint[i];
 
 
// Read ParticleType if present
 
    if (ParticleTypeInFile == TRUE) {
 
      /* Read ParticleType into temporary buffer and Copy to ParticleType. */
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        assert( file_dsp_id != h5_error);
 
      if (io_log) fprintf(log_fptr,"H5Dopen  with Name = particle_type\n");
 
      dset_id =  H5Dopen(group_id, "particle_type");
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        assert( dset_id != h5_error );
 
      h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempint);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      for (i = 0; i < NumberOfParticles; i++)
        ParticleType[i] = tempint[i];
 
      for (i = 0; i < NumberOfParticles; i++)
        if (ParticleType[i] < PARTICLE_TYPE_GAS ||
            ParticleType[i] > NUM_PARTICLE_TYPES-1) {
          fprintf(stderr, "file: %s: particle %"ISYM" has unknown type %"ISYM"\n",
                  name, i, ParticleType[i]);
          return FAIL;
        }
 
    } else {
 
      /* Otherwise create the type. */
 
      for (i = 0; i < NumberOfParticles; i++)
        ParticleType[i] = ReturnParticleType(i);
 
    }
 
 
    /* Read ParticleAttributes. */
 
    for (j = 0; j < NumberOfParticleAttributes; j++) {
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        assert( file_dsp_id != h5_error );
 
      if (io_log) fprintf(log_fptr,"H5Dopen with Name = %s\n",ParticleAttributeLabel[j]);
 
      dset_id =  H5Dopen(group_id, ParticleAttributeLabel[j]);
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        assert( dset_id != h5_error );
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        assert( h5_status != h5_error );
 
      for (i = 0; i < NumberOfParticles; i++)
	ParticleAttribute[j][i] = float(temp[i]);
 
    }
 
    delete [] temp;
    delete [] tempint;
 
  } // end: if (MyProcessorNumber == ProcessorNumber)
  } // end: if (NumberOfParticles > 0)
 
  /* Close file. */
 
  if ( (MyProcessorNumber == ProcessorNumber) &&
       (NumberOfParticles > 0 || NumberOfBaryonFields > 0) ){
 
     h5_status = H5Gclose(group_id);
       if (io_log) fprintf(log_fptr, "H5Gclose: %"ISYM"\n", h5_status);
       assert( h5_status != h5_error );
 
     h5_status = H5Fclose(file_id);
       if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
       assert( h5_status != h5_error );
  }
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    if (io_log) fclose(log_fptr);
  }
 
  /* 4) Read gravity info */
 
  if (SelfGravity)
    if (fscanf(fptr, "GravityBoundaryType = %"ISYM"\n",&GravityBoundaryType) != 1) {
      fprintf(stderr, "Error reading GravityBoundaryType.\n");
      return FAIL;
    }
 
  return SUCCESS;
 
}

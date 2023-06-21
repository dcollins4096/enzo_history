/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//  Input a grid from file pointer fpt
//

#include <string.h>
#include <stdio.h>
#include "mfhdf.h"
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
#ifdef USE_FLEXIO
#include "IO.hh"
#include "IEEEIO.hh"
#endif /* USE_FLEXIO */

/* function prototypes */

int ReadListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
int ReadField(float32 *pointer, int32 Size[], int size, int sd_id, 
	      int32 &sds_index, char *name, int DataFileType);

/* Raw file pointer and FlexIO reader object. */

static FILE *Rawfptr;
#ifdef USE_FLEXIO
static IObase *FlexIOreader;
#endif


int grid::ReadGrid(FILE *fptr)
{
 
  /* declarations */

  int i, j, k, dim, field, DataFileType = 0;
  int32 TempIntArray[MAX_DIMENSION], ActiveDim[MAX_DIMENSION];
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
  int32 sd_id, sds_id, sds_index = 0, start[] = {0, 0, 0},
        num_type, attributes, TempInt;

  /* make sure quantities defined at least for 3d */

  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }

  /* Read general grid class data */

  if (fscanf(fptr, "GridRank = %d\n", &GridRank) != 1) {
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

  if (fscanf(fptr, "SubgridsAreStatic = %d\n", &SubgridsAreStatic) != 1) {
    fprintf(stderr, "Error reading SubgridsAreStatic.\n");
    return FAIL;
  }

  /* Compute Flux quantities */

  this->PrepareGridDerivedQuantities();

  /* Read baryon field quantities. */

  if (fscanf(fptr, "NumberOfBaryonFields = %d\n", 
	     &NumberOfBaryonFields) != 1) {
    fprintf(stderr, "Error reading NumberOfBaryonFields.\n");
    return FAIL;
  }

  if (NumberOfBaryonFields > 0) {

    fscanf(fptr, "FieldType = ");
    if (ReadListOfInts(fptr, NumberOfBaryonFields, FieldType) == FAIL) {
      fprintf(stderr, "Error reading FieldType.\n");
      return FAIL;
    }

    if (fscanf(fptr, "BaryonFileName = %s\n", name) != 1) {
      fprintf(stderr, "Error reading BaryonFileName.\n");
      return FAIL;
    }

    fscanf(fptr, "CourantSafetyNumber    = %"FSYM"\n", &CourantSafetyNumber);
    fscanf(fptr, "PPMFlatteningParameter = %d\n", &PPMFlatteningParameter);
    fscanf(fptr, "PPMDiffusionParameter  = %d\n", &PPMDiffusionParameter);
    fscanf(fptr, "PPMSteepeningParameter = %d\n", &PPMSteepeningParameter);

    if (MyProcessorNumber == ProcessorNumber) {

    /* Set the data file type: 1 - HDF, 2 - Raw, 3 - FlexIO and open file. */

    DataFileType = 1;
    if (strstr(name, ".raw") != NULL) DataFileType = 2;
    if (strstr(name, ".flexio") != NULL) DataFileType = 3;

    /* 1) HDF Type */

    if (DataFileType == 1) {

      if ((sd_id = SDstart(name, DFACC_RDONLY)) == HDF_FAIL) {
	fprintf(stderr, "Error opening file %s.\n", name);
	return FAIL;
      }

      sds_id = SDselect(sd_id, sds_index++);
      while (SDiscoordvar(sds_id)) {
	SDendaccess(sds_id);
	sds_id = SDselect(sd_id, sds_index++);
      }
      if (SDgetinfo(sds_id, dummy, &TempInt, TempIntArray, &num_type, 
		    &attributes) == HDF_FAIL) {
	fprintf(stderr, "Error reading dims from %s.\n", name);
	return FAIL;
      }
      SDendaccess(sds_id);
      sds_index--;

      /* check rank against this grid */

      if (TempInt != GridRank) {
	fprintf(stderr, "HDF rank (%d) does not match GridRank.\n", TempInt);
	return FAIL;
      }
    }

    /* 2) Raw Type */

    if (DataFileType == 2)
      if ((Rawfptr = fopen(name, "rb")) == NULL) {
	fprintf(stderr, "Error opening raw file %s.\n", name);
	return FAIL;
      }

    /* 3) FlexIO */

#ifdef USE_FLEXIO
    if (DataFileType == 3) {
      FlexIOreader = new IEEEIO(name, IObase::Read);
      if (!FlexIOreader->isValid()) {
	fprintf(stderr, "Error opening FlexIO file %s.\n", name);
	return FAIL;
      }
    }
#endif /* USE_FLEXIO */

    /* fill in ActiveDim for dims up to 3d */

    for (dim = 0; dim < 3; dim++)
      ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;

    /* check dimensions of HDF file against this grid
       (note: we don't bother to check the coordinate arrays)  */

    int size = 1, active_size = 1;
    for (dim = 0; dim < GridRank; dim++) {
      if (TempIntArray[GridRank-dim-1] != ActiveDim[dim] && 
	  DataFileType == 1) {
	fprintf(stderr, "HDF file dimensions do not match GridDimensions.\n");
	return FAIL;
      }
      size *= GridDimension[dim];
      active_size *= ActiveDim[dim];
    }

    /* allocate temporary space */

    float32 *temp = new float32[active_size];

    /* loop over fields, reading each one */

    for (field = 0; field < NumberOfBaryonFields; field++) {

      /* get data into temporary array */

      if (ReadField(temp, TempIntArray, active_size, sd_id, sds_index, name, 
		    DataFileType) == FAIL) {
	fprintf(stderr, "Error reading field %d.\n", field);
	return FAIL;
      }

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

  if (fscanf(fptr, "NumberOfParticles = %d\n", &NumberOfParticles) != 1) {
    fprintf(stderr, "error reading NumberOfParticles.\n");
    return FAIL;
  }

//  NumberOfParticles = 0;

  if (NumberOfParticles > 0) {

    /* Read particle file name. */

    if (fscanf(fptr, "ParticleFileName = %s\n", name) != 1) {
      fprintf(stderr, "Error reading ParticleFileName.\n");
      return FAIL;
    }

    if (MyProcessorNumber == ProcessorNumber) {

    /* Set the data file type: 1 - HDF, 2 - Raw, 3 - FlexIO and open file. */

    DataFileType = 1;
    if (strstr(name, ".raw") != NULL) DataFileType = 2;
    if (strstr(name, ".flexio") != NULL) DataFileType = 3;

    /* Open file if not already done (note: particle name must = grid name). */

    if (NumberOfBaryonFields == 0) {

      if (DataFileType == 1)
	if ((sd_id = SDstart(name, DFACC_RDONLY)) == HDF_FAIL) {
	  fprintf(stderr, "Error opening file %s.\n", name);
	  return FAIL;
	}

      if (DataFileType == 2)
	if ((Rawfptr = fopen(name, "rb")) == NULL) {
	  fprintf(stderr, "Error opening raw file %s.\n", name);
	  return FAIL;
	}

#ifdef USE_FLEXIO
      if (DataFileType == 3){
	FlexIOreader = new IEEEIO(name, IObase::Read);
	if (!FlexIOreader->isValid()) {
	  fprintf(stderr, "Error opening FlexIO file %s.\n", name);
	  return FAIL;
	}
      }
#endif /* USE_FLEXIO */

    } // end: if (NumberOfBaryonFields == 0)

    /* Check dims if this is an HDF file. */

    if (DataFileType == 1) {

      /* Read dims.
	 If Rank != 1, we may have just read some other field SDS.  If so,
	 then try again. */
      
      TempInt = 0;
      while (TempInt != 1) {
	sds_id = SDselect(sd_id, sds_index++);
	while (SDiscoordvar(sds_id)) {
	  SDendaccess(sds_id);
	  sds_id = SDselect(sd_id, sds_index++);
	}
	if (SDgetinfo(sds_id, dummy, &TempInt, TempIntArray, &num_type, 
		      &attributes) == HDF_FAIL) {
	  fprintf(stderr, "Error reading dims from %s.\n", name);
	  return FAIL;
	}
	SDendaccess(sds_id);
      }
      sds_index--; 

      /* Check dims. */

      if (TempInt != 1 || TempIntArray[0] != NumberOfParticles) {
	fprintf(stderr, "HDF particle dims do not match NumberOfParticles.\n");
	fprintf(stderr, "  (HDF dim[0] = %d, NumberOfParticles = %d)\n",
		int(TempIntArray[0]), NumberOfParticles);
	return FAIL;  
      }

    } // end: if (RawDataFile)
    
    /* Allocate rooms for particles. */

    this->AllocateNewParticles(NumberOfParticles);

    /* Set the Data type we are using here. */

    int32 HDFDataType = (sizeof(FLOAT) == 4) ? DFNT_FLOAT32 : DFNT_FLOAT64;
    if (sizeof(FLOAT) == 16)
      HDFDataType = DFNT_FLOAT128;
    TempIntArray[0] = int32(NumberOfParticles);

    /* Create a temporary buffer (32 bit or twice the size for 64). */

    float32 *temp = NULL;
    if (num_type != HDFDataType) {
      if (num_type == DFNT_FLOAT64)
	temp = new float32[NumberOfParticles*2];
      if (num_type == DFNT_FLOAT128)
	temp = new float32[NumberOfParticles*4];
    }
    if (temp == NULL)
      temp = new float32[NumberOfParticles];

    /* Read ParticlePosition (use temporary buffer). */ 
      
    for (dim = 0; dim < GridRank; dim++) {

      if (DataFileType == 1 && num_type == HDFDataType) {

	/* same data type: just read. */

	sds_id = SDselect(sd_id, sds_index++);
	while (SDiscoordvar(sds_id)) {
	  SDendaccess(sds_id);
	  sds_id = SDselect(sd_id, sds_index++);
	}
	if (SDreaddata(sds_id, start, (int32 *) NULL, TempIntArray, 
		       (VOIDP) ParticlePosition[dim]) == HDF_FAIL) {
	  fprintf(stderr, "Error reading data from %s.\n", name);
	  return FAIL;
	}
	SDendaccess(sds_id);
	
      } else {

	/* convert data: Read into temporary buffer and copy. */

	if (ReadField(temp, TempIntArray, NumberOfParticles, sd_id, sds_index, 
		      name, DataFileType) == FAIL) {
	  fprintf(stderr, "Error reading ParticlePosition %d\n", dim);
	  return FAIL;
	}

	float64 *temp64 = (float64 *) temp;
	long_double *temp128 = (long_double *) temp;

	if (num_type == DFNT_FLOAT32)
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticlePosition[dim][i] = FLOAT(temp[i]);
	if (num_type == DFNT_FLOAT64)
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticlePosition[dim][i] = FLOAT(temp64[i]);
	if (num_type == DFNT_FLOAT128)
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticlePosition[dim][i] = FLOAT(temp128[i]);
      }
    }

    /* Read ParticleVelocity. */

    for (dim = 0; dim < GridRank; dim++) {
      if (ReadField(temp, TempIntArray, NumberOfParticles, sd_id, sds_index, 
		    name, DataFileType) == FAIL) {
	fprintf(stderr, "Error reading ParticleVelocity %d\n", dim);
	return FAIL;
      }
      for (i = 0; i < NumberOfParticles; i++)
	ParticleVelocity[dim][i] = float(temp[i]);
    }

    /* Read ParticleMass into temporary buffer and Copy to ParticleMass. */

    if (ReadField(temp, TempIntArray, NumberOfParticles, sd_id, sds_index, 
		  name, DataFileType) == FAIL)
      return FAIL;
    for (i = 0; i < NumberOfParticles; i++)
      ParticleMass[i] = float(temp[i]);

    /* Read ParticleNumber into temporary buffer and Copy to ParticleNumber. */

    int32 *tempint = new int32[NumberOfParticles];

    if (DataFileType == 1) {
      sds_id = SDselect(sd_id, sds_index++);
      while (SDiscoordvar(sds_id)) {
	SDendaccess(sds_id);
	sds_id = SDselect(sd_id, sds_index++);
      }
      if (SDreaddata(sds_id, start, NULL, TempIntArray, (VOIDP) tempint) ==
	  HDF_FAIL) {
	fprintf(stderr, "Error reading data from %s.\n", name);
	return FAIL;
      }
      SDendaccess(sds_id);
    }

    if (DataFileType == 2)
      if (fread((VOIDP) tempint, sizeof(int32), NumberOfParticles,
		Rawfptr) != NumberOfParticles) {
	perror("fread error");
	return FAIL;
      }

#ifdef USE_FLEXIO
    if (DataFileType == 3) {
      IObase::DataType InputDataType;
      int InputRank, InputDims[3];
      FlexIOreader->readInfo(InputDataType, InputRank, InputDims);
      if (InputDims[0] != NumberOfParticles || InputRank != 1 ||
	  InputDataType != IObase::Int32) {
	fprintf(stderr, "FlexIO dim or type mismatch: particle index\n");
	return FAIL;
      }
      FlexIOreader->read((VOIDP) temp);
    }
#endif /* USE_FLEXIO */

    for (i = 0; i < NumberOfParticles; i++)
      ParticleNumber[i] = int(tempint[i]);

    /* Read ParticleAttributes. */

    for (j = 0; j < NumberOfParticleAttributes; j++) {
      if (ReadField(temp, TempIntArray, NumberOfParticles, sd_id, sds_index, 
		    name, DataFileType) == FAIL) {
	fprintf(stderr, "Error reading ParticleAttribute %d\n", j);
	return FAIL;
      }
      for (i = 0; i < NumberOfParticles; i++)
	ParticleAttribute[j][i] = float(temp[i]);
    }

    delete [] temp;
    delete [] tempint;

  } // end: if (MyProcessorNumber == ProcessorNumber)
  } // end: if (NumberOfParticles > 0)

  /* Close file. */

  if (MyProcessorNumber == ProcessorNumber) {
    if (DataFileType == 1)
      SDend(sd_id);
    if (DataFileType == 2)
      fclose(Rawfptr);
#ifdef USE_FLEXIO
    if (DataFileType == 3)
      delete FlexIOreader;
#endif /* USE_FLEXIO */
  }

  /* 4) Read gravity info */

  if (SelfGravity)
    if (fscanf(fptr, "GravityBoundaryType = %d\n",&GravityBoundaryType) != 1) {
      fprintf(stderr, "Error reading GravityBoundaryType.\n");
      return FAIL;
    }

  return SUCCESS;
}


int ReadField(float32 *temp, int32 Size[], int size, int sd_id, 
	      int32 &sds_index, char *name, int DataFileType)
{

  int32 sds_id, start[] = {0, 0, 0};

  /* 1) Read from HDF file, skipping past any coordinate axis. */

  if (DataFileType == 1) {
    sds_id = SDselect(sd_id, sds_index++);
    while (SDiscoordvar(sds_id)) {
      SDendaccess(sds_id);
      sds_id = SDselect(sd_id, sds_index++);
    }
    if (SDreaddata(sds_id, start, (int32 *) NULL, Size, (VOIDP) temp)
	== HDF_FAIL) {
      fprintf(stderr, "Error reading data from %s.\n", name);
      return FAIL;
    }
    SDendaccess(sds_id);
  }

  /* 2) Read from binary file. */

  if (DataFileType == 2)
    if (fread((VOIDP) temp, sizeof(float32), size, Rawfptr) != size) {
      perror("fread error");
      return FAIL;
    }

  /* 3) Read From FlexIO file. */

#ifdef USE_FLEXIO
  if (DataFileType == 3) {
    IObase::DataType InputDataType;
    int InputRank, InputDims[3] = {0,0,0};
    FlexIOreader->readInfo(InputDataType, InputRank, InputDims);
    //    for (int dim = 0; dim < InputRank; dim++)
    //      if (InputDims[dim] != Size[InputRank-1-dim]) {
    //	fprintf(stderr, "FlexIO dims mismatch: %d %d (dim = %d)\n", 
    //		InputDims[dim], Size[InputRank-1-dim], dim);
    //	return FAIL;
    //      }
    FlexIOreader->read((VOIDP) temp);
  }
#endif /* USE_FLEXIO */

  return SUCCESS;
}

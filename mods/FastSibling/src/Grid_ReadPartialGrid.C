/***********************************************************************
/
/  GRID CLASS (READS THE GRID DATA FROM A FILE CONTAINING THE WHOLE TOPGRID)
/
/  written by: Greg Bryan
/  date:       
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

extern int ParticleTypeInFile; // declared and set in ReadParameterFile

/* function prototypes */

int ReadListOfInts(FILE *fptr, int N, int nums[]);
int ReadPartialField(float32 *temp, int32 Size[], int32 start[], int Rank, 
		     int size, int sd_id, int32 &sds_index, char *name, 
		     int DataFileType);


int grid::ReadPartialGrid(FILE *fptr)
{
 
  /* declarations */
  
  int i, j, k, dim, field, DataFileType = 0, One = 1;
  int32 TempIntArray[MAX_DIMENSION], ActiveDim[MAX_DIMENSION];
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
  int32 sd_id, sds_id, sds_index = 0, ReadStartIndex[] = {0, 0, 0},
    num_type, attributes, TempInt;

  /* Read 10 lines to skip over the information we're not interested
     in (which has already been read in) */

  for (i = 0; i < 10; i++)
    fgets(dummy, MAX_LINE_LENGTH, fptr);

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

    fscanf(fptr, "CourantSafetyNumber    = %f\n", &CourantSafetyNumber);
    fscanf(fptr, "PPMFlatteningParameter = %d\n", &PPMFlatteningParameter);
    fscanf(fptr, "PPMDiffusionParameter  = %d\n", &PPMDiffusionParameter);
    fscanf(fptr, "PPMSteepeningParameter = %d\n", &PPMSteepeningParameter);

    if (MyProcessorNumber == ProcessorNumber) {

      //      fprintf(stderr, "P(%d): reading %d %d %\n", MyProcessorNumber, GridDimension[0],
      //	      GridDimension[1], GridDimension[2]);

    /* Set the data file type: 1 - HDF, 2 - Raw, 3 - FlexIO, 4 - FORTIO,
       and open file. */

    DataFileType = 1;
    if (strstr(name, ".raw") != NULL) DataFileType = 2;
    if (strstr(name, ".flexio") != NULL) DataFileType = 3;
    if (strstr(name, ".fort") != NULL) DataFileType = 4;

    if (DataFileType != 1) {
      fprintf(stderr, "ReadPartialGrid only supports HDF type.\n");
      return FAIL;
    }

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

    /* fill in ActiveDim for dims up to 3d */

    int size = 1, active_size = 1;
    for (dim = 0; dim < 3; dim++) {
      ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
      size *= GridDimension[dim];
      active_size *= ActiveDim[dim];
    }

    /* allocate temporary space */

    float32 *temp = new float32[active_size];

    /* Calculate the region to read (note that we reverse the dimension
       order for HDF). */

    for (dim = 0; dim < GridRank; dim++) {
      ReadStartIndex[GridRank-dim-1] = 
	nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/ CellWidth[dim][0]);
      TempIntArray[GridRank-dim-1] = ActiveDim[dim];
    }
    /* fprintf(stderr, "P(%d) Reading %d %d %d  %d %d %d\n", MyProcessorNumber,
	    ReadStartIndex[0], ReadStartIndex[1], ReadStartIndex[2],
	    ActiveDim[0], ActiveDim[1], ActiveDim[2]); */

    /* loop over fields, reading each one */

    for (field = 0; field < NumberOfBaryonFields; field++) {

      /* get data into temporary array */

      if (ReadPartialField(temp, TempIntArray, ReadStartIndex, GridRank, 
		    active_size, sd_id, sds_index, name, 
		    DataFileType) == FAIL) {
	fprintf(stderr, "Error reading field %d.\n", field);
	return FAIL;
      }

      /* copy active region into whole grid unless not read grid data */

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

    /* If necessary, read extra 3d fields (temperature and dm density). */

    if (DataFileType == 4) {
      if (ComovingCoordinates)
	ReadPartialField(temp, TempIntArray, ReadStartIndex, GridRank, 
		  active_size, sd_id, sds_index, name, DataFileType);
      if (SelfGravity && NumberOfParticles > 0)
	ReadPartialField(temp, TempIntArray, ReadStartIndex, GridRank, 
		  active_size, sd_id, sds_index, name, DataFileType);
    }

    delete [] temp;

  }  // end: if (MyProcessorNumber == ProcessorNumber)
  }  // end read baryon fields

  /* 3) Read particle info */

  if (fscanf(fptr, "NumberOfParticles = %d\n", &NumberOfParticles) != 1) {
    fprintf(stderr, "error reading NumberOfParticles.\n");
    return FAIL;
  }

  if (NumberOfParticles > 0) {

    /* Read particle file name. */

    if (fscanf(fptr, "ParticleFileName = %s\n", name) != 1) {
      fprintf(stderr, "Error reading ParticleFileName.\n");
      return FAIL;
    }

    /* Determine number of particles we are actually going to read.  Note that
       these particles may not live on this grid but this is ok because we
       do a rebuild immediately after which moves everything around. */

    int Count = int(NumberOfParticles/NumberOfProcessors);
    int ParticleStartIndex = Count*MyProcessorNumber;
    if (MyProcessorNumber == NumberOfProcessors-1)
      NumberOfParticles = NumberOfParticles - ParticleStartIndex;
    else
      NumberOfParticles = Count;

    if (MyProcessorNumber == ProcessorNumber) {

    /* Set the data file type: 1 - HDF, 2 - Raw, 3 - FlexIO and open file. */

    DataFileType = 1;
    if (strstr(name, ".raw") != NULL) DataFileType = 2;
    if (strstr(name, ".flexio") != NULL) DataFileType = 3;
    if (strstr(name, ".fort") != NULL) DataFileType = 4;

    /* Open file if not already done (note: particle name must = grid name). */

    if (NumberOfBaryonFields == 0) {

      if (DataFileType == 1)
	if ((sd_id = SDstart(name, DFACC_RDONLY)) == HDF_FAIL) {
	  fprintf(stderr, "Error opening file %s.\n", name);
	  return FAIL;
	}

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

      if (TempInt != 1) {
	fprintf(stderr, "HDF particle dims do not match NumberOfParticles.\n");
	fprintf(stderr, "  (HDF dim[0] = %d, NumberOfParticles = %d)\n",
		int(TempIntArray[0]), NumberOfParticles);
	return FAIL;  
      }

    } // end: if (HDFDataFile)

    ReadStartIndex[0] = ParticleStartIndex;
    TempIntArray[0] = int32(NumberOfParticles);
    /* fprintf(stderr, "P(%d) particle start %d  count %d\n", 
       MyProcessorNumber, ReadStartIndex[0], NumberOfParticles); */
    
    /* Allocate rooms for particles. */

    this->AllocateNewParticles(NumberOfParticles);

    /* Set the Data type we are using here. */

    int32 HDFDataType = (sizeof(FLOAT) == 4) ? DFNT_FLOAT32 : DFNT_FLOAT64;
    if (sizeof(FLOAT) == 16)
      HDFDataType = DFNT_FLOAT128;

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
	if (SDreaddata(sds_id, ReadStartIndex, (int32 *) NULL, TempIntArray, 
		       (VOIDP) ParticlePosition[dim]) == HDF_FAIL) {
	  fprintf(stderr, "Error reading data from %s.\n", name);
	  return FAIL;
	}
	SDendaccess(sds_id);
	
      } else {

	/* convert data: Read into temporary buffer and copy. */

	if (ReadPartialField(temp, TempIntArray, ReadStartIndex, 1, 
		      NumberOfParticles, sd_id, sds_index, 
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
      if (ReadPartialField(temp, TempIntArray, ReadStartIndex, 1, 
		    NumberOfParticles, sd_id, sds_index, 
		    name, DataFileType) == FAIL) {
	fprintf(stderr, "Error reading ParticleVelocity %d\n", dim);
	return FAIL;
      }
      for (i = 0; i < NumberOfParticles; i++)
	ParticleVelocity[dim][i] = float(temp[i]);
    }

    /* Read ParticleMass into temporary buffer and Copy to ParticleMass. */

    if (ReadPartialField(temp, TempIntArray, ReadStartIndex, 1, 
		  NumberOfParticles, sd_id, sds_index, 
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
      if (SDreaddata(sds_id, ReadStartIndex, NULL, TempIntArray, 
		     (VOIDP) tempint) == HDF_FAIL) {
	fprintf(stderr, "Error reading data from %s.\n", name);
	return FAIL;
      }
      SDendaccess(sds_id);
    }

    for (i = 0; i < NumberOfParticles; i++)
      ParticleNumber[i] = int(tempint[i]);

    /* Read ParticleType into temporary buffer and Copy to ParticleNumber. */

    if (ParticleTypeInFile == TRUE) {

      if (DataFileType == 1) {
	sds_id = SDselect(sd_id, sds_index++);
	while (SDiscoordvar(sds_id)) {
	  SDendaccess(sds_id);
	  sds_id = SDselect(sd_id, sds_index++);
	}
	if (SDreaddata(sds_id, ReadStartIndex, NULL, TempIntArray, 
		       (VOIDP) tempint) == HDF_FAIL) {
	  fprintf(stderr, "Error reading data from %s.\n", name);
	  return FAIL;
	}
	SDendaccess(sds_id);
      }

      for (i = 0; i < NumberOfParticles; i++)
	ParticleNumber[i] = int(tempint[i]);

    } else {
      
      /* Otherwise create the type. */

      for (i = 0; i < NumberOfParticles; i++)
	ParticleType[i] = ReturnParticleType(i);

    }


    /* Read ParticleAttributes. */

#define NO_RESTART_WITH_ATTRIBUTES
    for (j = 0; j < NumberOfParticleAttributes; j++) {
#ifdef RESTART_WITH_ATTRIBUTES
      for (i=0; i < NumberOfParticles; i++)
	ParticleAttribute[j][i] = 0;
#else
      if (ReadPartialField(temp, TempIntArray, ReadStartIndex, 1, 
		    NumberOfParticles, sd_id, sds_index, 
		    name, DataFileType) == FAIL) {
	fprintf(stderr, "Error reading ParticleAttribute %d\n", j);
	return FAIL;
      }
      for (i = 0; i < NumberOfParticles; i++)
	ParticleAttribute[j][i] = float(temp[i]);
#endif
    }

    delete [] temp;
    delete [] tempint;

  } // end: if (MyProcessorNumber == ProcessorNumber)
  } // end: if (NumberOfParticles > 0)

  /* Close file. */

  if (MyProcessorNumber == ProcessorNumber) {
    if (DataFileType == 1)
      SDend(sd_id);
  }

  return SUCCESS;
}


int ReadPartialField(float32 *temp, int32 Size[], int32 start[], int Rank, 
		     int size, int sd_id, int32 &sds_index, char *name, 
		     int DataFileType)
{

  int32 sds_id;

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

  return SUCCESS;
}

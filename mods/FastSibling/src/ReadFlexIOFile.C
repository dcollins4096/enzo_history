/***********************************************************************
/
/  READ A FlexIO FILE
/
/  written by: Greg Bryan
/  date:       August, 1999
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine reads the parameter file in the argument and sets parameters
//   based on it.

#include <stdio.h>
#include <df.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
#ifdef USE_FLEXIO
#include "IO.hh"
#include "IEEEIO.hh"
#endif

/* FlexIO reader object */

#ifdef USE_FLEXIO
static IObase *FlexIOreader = NULL;
static char *LastUsedFileName = NULL;
#endif

/* function prototypes */

int ReadFlexIOFile(char *name, int Rank, int Dim[], int StartIndex[],
		   int EndIndex[], int BufferOffset[], float *buffer)
{

#ifdef USE_FLEXIO

  /* declarations */

  int dim, i, j, k, TempInt;
  int TempIntArray[MAX_DIMENSION];

  /* Error check name */

  if (name == NULL)
    return FAIL;

  /* If the name is different, then close the file. */

  if (LastUsedFileName != NULL)
    if (strcmp(name, LastUsedFileName) != 0) {
      delete FlexIOreader;
      delete [] LastUsedFileName;
      LastUsedFileName = NULL;
      FlexIOreader = NULL;
    }

  /* Open file if needed. */

  if (FlexIOreader == NULL) {
    FlexIOreader = new IEEEIO(name, IObase::Read);
    if (!FlexIOreader->isValid()) {
      fprintf(stderr, "Error opening FlexIO file %s.\n", name);
      return FAIL;
    }
  }

  /* Set last name if not already set. */

  if (LastUsedFileName == NULL) {
    LastUsedFileName = new char[strlen(name)];
    strcpy(LastUsedFileName, name);
  }

  /* Read dimensions. */

  IObase::DataType InputDataType;
  FlexIOreader->readInfo(InputDataType, TempInt, TempIntArray);

  /* Error check. */

  if (Rank < 1 || Rank > 3) {
    fprintf(stderr, "Rank %d not supported.\n", Rank);
    return FAIL;
  }

  if (TempInt != Rank) {
    fprintf(stderr, "Rank mismatch in %s.\n", name);
    return FAIL;
  }

  /* Check dimensions. */

  if (ParallelRootGridIO == FALSE)
    for (dim = 0; dim < Rank; dim++)
      if (TempIntArray[dim] != (EndIndex[dim]-StartIndex[dim]+1)) {
	fprintf(stderr, "Dimension mismatch in %s.\n", name);
	return FAIL;
      }

  /* Compute size offield. */

  int size = 1;
  for (dim = 0; dim < Rank; dim++)
    size *= int(TempIntArray[dim]);

  /* Allocate space for temp buffer. */

  float32 *tempbuffer = new float32[size];

  /* Read FlexIO file. */

  FlexIOreader->read((VOIDP) tempbuffer);

  /* clear buffer (primarily to prevent errors in unused area). */

  size = 1;
  for (dim = 0; dim < Rank; dim++)
    size *= Dim[dim];
  for (i = 0; i < size; i++)
    buffer[i] = 0;
    
  /* Copy field into real array. */

  if (Rank == 1)
    for (i = StartIndex[0]; i <= EndIndex[0]; i++)
      buffer[i] = tempbuffer[i-StartIndex[0]+BufferOffset[0]];

  if (Rank == 2)
    for (j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	buffer[j*Dim[0] + i] =
	  tempbuffer[(j-StartIndex[1]+BufferOffset[1])*TempIntArray[0] + 
		     (i-StartIndex[0]+BufferOffset[0])];

  if (Rank == 3)
    for (k = StartIndex[2]; k <= EndIndex[2]; k++)
      for (j = StartIndex[1]; j <= EndIndex[1]; j++)
	for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	  buffer[k*Dim[0]*Dim[1] + j*Dim[0] + i] =
	    tempbuffer[(k-StartIndex[2]+BufferOffset[2])*
		                          TempIntArray[0]*TempIntArray[1] +
	               (j-StartIndex[1]+BufferOffset[1])*TempIntArray[0] +
	               (i-StartIndex[0]+BufferOffset[0])];

  /* clean up */

  delete tempbuffer;

#endif /* USE_FLEXIO */

  return SUCCESS;
}

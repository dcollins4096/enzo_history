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
/***********************************************************************
/
/  READ AN HDF FILE
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine reads the parameter file in the argument and sets parameters
//   based on it.

#ifdef USE_HDF4

#include <stdio.h>
#include <df.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* function prototypes */

int ReadHDF4File(char *name, int Rank, int Dim[], int StartIndex[],
		int EndIndex[], int BufferOffset[], float *buffer,
		float32 **tempbuffer)
{
  /* declarations */

  int dim, i, j, k, TempInt;
  int32 TempIntArray[MAX_DIMENSION];

  /* Error check name */

  if (name == NULL)
    return FAIL;

  /* Read dimensions. */

  if (DFSDgetdims(name, &TempInt, TempIntArray, MAX_DIMENSION) == HDF_FAIL) {
    fprintf(stderr, "Error reading dims from %s.\n", name);
    return FAIL;
  }

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

  /* Compute size of HDF field. */

  int size = 1;
  for (dim = 0; dim < Rank; dim++)
    size *= (EndIndex[dim]-StartIndex[dim]+1);

  /* Allocate space for temp buffer. */

  (*tempbuffer) = new float32[size];

  /* Read HDF file. */

  if (ParallelRootGridIO == TRUE) {

    /* If doing parallel IO, then just read this grid. */

    int32 start[3], slab_size[3];
    for (dim = 0; dim < Rank; dim++) {
      start[dim] = BufferOffset[dim]+1;  // one-based
      slab_size[dim] = EndIndex[dim]-StartIndex[dim]+1;
      TempIntArray[dim] = slab_size[dim];
    }
    /*    fprintf(stderr, "P(%d): start = %d %d %d   slab_size = %d %d %d\n",
	    MyProcessorNumber, start[0], start[1], start[2], slab_size[0],
	    slab_size[1], slab_size[2]); */
    if (DFSDreadslab(name, start, slab_size, start, (VOIDP) (*tempbuffer), 
		     TempIntArray) == HDF_FAIL) {
      fprintf(stderr, "Error reading slab data from %s.\n", name);
      HEprint(stderr, 0);
      return FAIL;
    }
    /*    fprintf(stderr, "P(%d): TempIntArray = %d %d %d\n", 
	    MyProcessorNumber, TempIntArray[0], TempIntArray[1], 
	    TempIntArray[2]); */

  } else

    /* otherwise, read whole file. */

    if (DFSDgetdata(name, TempInt, TempIntArray, (VOIDP) (*tempbuffer)) 
	== HDF_FAIL) {
      fprintf(stderr, "Error reading data from %s.\n", name);
      return FAIL;
    }

  /* If buffer is not defined, then just return w/o clearing (*tempbuffer). */

  if (buffer == NULL)
    return SUCCESS;

  /* clear buffer (primarily to prevent errors in unused area). */

  size = 1;
  for (dim = 0; dim < Rank; dim++)
    size *= Dim[dim];
  for (i = 0; i < size; i++)
    buffer[i] = 0;
    
  /* Copy field into real array. */

  if (Rank == 1)
    for (i = StartIndex[0]; i <= EndIndex[0]; i++)
      buffer[i] = (*tempbuffer)[i-StartIndex[0]];

  if (Rank == 2)
    for (j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	buffer[j*Dim[0] + i] =
	  (*tempbuffer)[(j-StartIndex[1])*TempIntArray[0] + 
		     (i-StartIndex[0])];

  if (Rank == 3)
    for (k = StartIndex[2]; k <= EndIndex[2]; k++)
      for (j = StartIndex[1]; j <= EndIndex[1]; j++)
	for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	  buffer[k*Dim[0]*Dim[1] + j*Dim[0] + i] =
	    (*tempbuffer)[(k-StartIndex[2])*TempIntArray[0]*TempIntArray[1] +
	               (j-StartIndex[1])*TempIntArray[0] +
	               (i-StartIndex[0])];

  /* clean up */

  delete [] (*tempbuffer);

  return SUCCESS;
}

#endif /* USE_HDF4 */

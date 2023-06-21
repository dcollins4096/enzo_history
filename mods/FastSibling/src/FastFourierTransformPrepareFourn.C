/***********************************************************************
/
/  COMPUTE A FAST FOURIER TRANSFORM (FORWARD OR INVERSE)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/      buffer - field to be FFTed
/      Rank   - rank of FFT
/      DimensionReal[] - declared dimensions of buffer
/      Dimension[]     - active dimensions of buffer
/      direction       - +1 forward, -1 inverse
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"

/* Defines */

extern "C" void FORTRAN_NAME(prefourn)(float *data, float *work, 
				       int *dim1, int *dim2, int *dim3,
				       int *rdim1, int *rdim2, int *rdim3,
				       int *rank, int *isign);
extern "C" void FORTRAN_NAME(fourn)(float *data, int *nn, int *ndim, 
				    int *isign);

/* Start routine. */

int FastFourierTransformPrepareFourn(float *buffer, int Rank, 
                                     int DimensionReal[], 
                                     int Dimension[], int direction, int type)
{

  int dim, Dim[MAX_DIMENSION], DimReal[MAX_DIMENSION];

  /* Error check. */

  if (Rank < 1 || Rank > 3) {
    fprintf(stderr, "Does not support Rank = %d\n", Rank);
    return FAIL;
  }

  /* Copy passed dims to Real dims to make sure they are at least 3d. */

  for (dim = 0; dim < Rank; dim++) {
    Dim[dim] = Dimension[dim];
    DimReal[dim] = DimensionReal[dim];
  }
  for (dim = Rank; dim < MAX_DIMENSION; dim++) {
    Dim[dim] = 1;
    DimReal[dim] = 1;
  }

  /* Error check for power of two. */

  for (dim = 0; dim < Rank; dim++)
    if (nint(pow(2.0, log10(float(Dim[dim]))/log10(2.0))) != Dim[dim]) {
      fprintf(stderr, "Dim[%d] = %d is not 2^n.\n", dim, Dim[dim]);
      return FAIL;
    }

  /* Compute size of required temporary buffer. */

  int size = 1;
  for (dim = 0; dim < Rank; dim++)
    size *= Dimension[dim];

  /* ---------------------------------------------------------- */
  /* a) real-to-complex or reverse.  This uses a complex-to-complex
        FFT so we must convert real part to complex; done in prefourn */

  if (type == REAL_TO_COMPLEX) {

    /* Allocate temporary buffer. */

    float *tempbuffer = new float[2*size];

    /* Prepare and carry out transform. */

    direction *= -1;
    FORTRAN_NAME(prefourn)(buffer, tempbuffer, &Dim[0], &Dim[1], &Dim[2],
			   &DimReal[0], &DimReal[1], &DimReal[2], &Rank, 
			   &direction);
    direction *= -1;

    /* clean up */

    delete [] tempbuffer;

  }

  /* ---------------------------------------------------------- */
  /* b) complex-to-complex.  Just pass the data to fourn */

  if (type == COMPLEX_TO_COMPLEX) {

    /* Error check. */

    for (dim = 0; dim < Rank; dim++)
      if (Dim[dim] != DimReal[dim]) {
	fprintf(stderr, "FFT does not support Dim != RealDim.\n");
	return FAIL;
      }

    direction *= -1;
    FORTRAN_NAME(fourn)(buffer, Dim, &Rank, &direction);
    direction *= -1;

    /* Scale on inverse. */

    if (direction == FFT_INVERSE) {
      float factor = 1.0/size;
      for (int i = 0; i < size*2; i++)
        buffer[i] *= factor;
    }

  }

  return SUCCESS;

}

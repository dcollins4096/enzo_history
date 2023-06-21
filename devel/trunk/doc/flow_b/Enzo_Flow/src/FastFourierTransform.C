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
#include "macros_and_parameters.h"

#ifdef GOT_FFT
# undef GOT_FFT
#endif

/* Defines */

int FastFourierTransformSGIMATH(float *buffer, int Rank, int DimensionReal[], 
				int Dimension[], int direction, int type);
int FastFourierTransformSGIMATHComplex(float *buffer, int Rank, 
				       int DimensionReal[], 
				       int Dimension[], int direction);
int FastFourierTransformVECLIB(float *buffer, int Rank, int DimensionReal[], 
			       int Dimension[], int direction);
int FastFourierTransformPrepareFourn(float *buffer, int Rank, 
				     int DimensionReal[], 
                                     int Dimension[], int direction, int type);

/* Start routine. */

int FastFourierTransform(float *buffer, int Rank, int DimensionReal[], 
			 int Dimension[], int direction, int type)
{

#if defined(IRIS4) && defined(SGI_MATH)

  /* --------------------------------------------------------------------- */
  /* Use SGI's real-to-complex routines. */

  if (FastFourierTransformSGIMATH(buffer, Rank, DimensionReal, 
				  Dimension, direction, type) == FAIL) {
    fprintf(stderr, "Error in FastFourierTransformSGIMATH.\n");
    return FAIL;
  }

#define GOT_FFT
#endif /* IRIS4 && SGI_MATH */


#if defined(IRIS4) && defined(SGI_MATH_COMPLEX)

  /* --------------------------------------------------------------------- */
  /* Use SGI's complex-to-complex routines. */

  if (FastFourierTransformSGIMATHComplex(buffer, Rank, DimensionReal, 
				  Dimension, direction) == FAIL) {
    fprintf(stderr, "Error in FastFourierTransformSGIMATHComplex.\n");
    return FAIL;
  }

#define GOT_FFT
#endif /* IRIS4 && SGI_MATH_COMPLEX */


#ifndef GOT_FFT

  /* --------------------------------------------------------------------- */
  /* Catchall: if there is no available FFT, use FOURN. */

  if (FastFourierTransformPrepareFourn(buffer, Rank, DimensionReal,
				       Dimension, direction, type) == FAIL) {
    fprintf(stderr, "Error in FastFourierTransformPrepareFourn.\n");
    return FAIL;
  }
  
#endif /* GOT_FFT */  

  return SUCCESS;
}

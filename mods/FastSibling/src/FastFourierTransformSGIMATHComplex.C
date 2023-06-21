/***********************************************************************
/
/  COMPUTE A FAST FOURIER TRANSFORM USING SGIMATH's COMPLEX-TO-XOMPLEX
/
/  written by: Greg Bryan
/  date:       June, 1995
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
#include "global_data.h"
#include "typedefs.h"

/* Defines */

int SGIMATHComplexTransform(float *buffer, int Rank, int Dimension[],
			    int direction);

#if defined(IRIS4) && defined(SGI_MATH_COMPLEX)
extern "C" float *cfft1di(int dim1, float *Work);
extern "C" int    cfft1d (int dir, int dim1, float *data, int skip, 
		 	  float *Work);
extern "C" void   cscal1d(int dim1, float scale, float *data, int skip);

extern "C" float *cfft2di(int dim1, int dim2, float *Work);
extern "C" int    cfft2d (int dir, int dim1, int dim2, float *data, 
			  int ld1, float *Work);
extern "C" void   cscal2d(int dim1, int dim2, float scale, float *data, 
			  int ld1);

extern "C" float *cfft3di(int dim1, int dim2, int dim3, float *Work);
extern "C" int    cfft3d (int dir, int dim1, int dim2, int dim3, float *data, 
			  int ld1, int ld2, float *Work);
extern "C" void   cscal3d(int dim1, int dim2, int dim3, float scale, 
			  float *data, int ld1, int ld2);
#endif /* IRIS4 && SGI_MATH_COMPLEX */

/* Start routine. */

int FastFourierTransformSGIMATHComplex(float *buffer, int Rank, 
				       int DimensionReal[], 
				       int Dimension[], int direction)
{

  int i, i1, j, j1, k, k1, dim;
  int Dim[MAX_DIMENSION], DimReal[MAX_DIMENSION];

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

  /* Compute size of required temporary buffer. */

  int size = 1;
  for (dim = 0; dim < Rank; dim++)
    size *= Dimension[dim];

  /* Allocate temporary buffer. */

  float *tempbuffer = new float[2*size];

  /* Forward transform: */

  if (direction == FFT_FORWARD) {

    /* copy real data into complex temporary. */

    for (k = 0; k < Dim[2]; k++)
      for (j = 0; j < Dim[1]; j++)
        for (i = 0; i < Dim[0]; i++) {
          *(tempbuffer + k*Dim[0]*Dim[1]*2 + j*Dim[0]*2 + i*2) =
            *(buffer + k*DimReal[0]*DimReal[1] + j*DimReal[0] + i);
          *(tempbuffer + k*Dim[0]*Dim[1]*2 + j*Dim[0]*2 + i*2 + 1) = 0.0;
        }

    /* Complex-to-complex transform. */

    if (SGIMATHComplexTransform(tempbuffer, Rank, Dimension, direction)
	== FAIL) {
      fprintf(stderr, "Error in SGIMATHComplexTransform (forward).\n");
      return FAIL;
    }

    /* Copy unique part of complex temporary back to buffer. */

    for (k = 0; k < Dim[2]; k++)
      for (j = 0; j < Dim[1]; j++)
        for (i = 0; i < Dim[0]+2; i++)
          *(buffer + k*DimReal[0]*DimReal[1] + j*DimReal[0] + i) =
            *(tempbuffer + k*Dim[0]*Dim[1]*2 + j*Dim[0]*2 + i);

  } // end of forward transform

  /* Inverse transform: */

  if (direction == FFT_INVERSE) {

    /* Recreate the other half of the complex data. */

    for (k = 0; k < Dim[2]; k++)
      for (j = 0; j < Dim[1]; j++) {

        /* First copy the 1/2 that is already there. */

        for (i = 0; i < Dim[0]+2; i++)
          *(tempbuffer + k*Dim[0]*Dim[1]*2 + j*Dim[0]*2 + i) =
            *(buffer + k*DimReal[0]*DimReal[1] + j*DimReal[0] + i);

        /* Use conjugate relations to make up other 1/2. */

        j1 = (Dim[1]-j) % Dim[1];
        k1 = (Dim[2]-k) % Dim[2];
        for (i = Dim[0]/2 + 1; i < Dim[0]; i++) {
          i1 =  Dim[0]-i;
          *(tempbuffer + k*Dim[0]*Dim[1]*2 + j*Dim[0]*2 + 2*i  ) =
            *(buffer + k1*DimReal[0]*DimReal[1] + j1*DimReal[0] + i1*2);
          *(tempbuffer + k*Dim[0]*Dim[1]*2 + j*Dim[0]*2 + 2*i+1) = -
            *(buffer + k1*DimReal[0]*DimReal[1] + j1*DimReal[0] + i1*2+1);
        }
      }

    /* Complex-to-complex transform. */

    if (SGIMATHComplexTransform(tempbuffer, Rank, Dimension, direction)
	== FAIL) {
      fprintf(stderr, "Error in SGIMATHComplexTransform (inverse).\n");
      return FAIL;
    }

    /* Copy real part of tempbuffer to buffer & scale. */

    for (k = 0; k < Dim[2]; k++)
      for (j = 0; j < Dim[1]; j++)
        for (i = 0; i < Dim[0]; i++)
          *(buffer + k*DimReal[0]*DimReal[1] + j*DimReal[0] + i) =
            *(tempbuffer + k*Dim[0]*Dim[1]*2 + j*Dim[0]*2 + i*2)
              /float(Dim[0]*Dim[1]*Dim[2]);

  } // end of inverse transform

  /* clean up */

  delete tempbuffer;

  return SUCCESS;

}



int SGIMATHComplexTransform(float *buffer, int Rank, int Dimension[],
			    int direction)
{

#if defined(IRIS4) && defined(SGI_MATH_COMPLEX)

  /* 1D transform. */

  if (Rank == 1) {
    float *Work = new float[2*(Dimension[0]+15)];
    cfft1di(Dimension[0], Work);
    cfft1d(direction*-1, Dimension[0], buffer, 1, Work);
    if (direction == FFT_INVERSE)          // scale by 1/n for inverse FFT
      cscal1d(Dimension[0], 1.0/float(Dimension[0]), buffer, 1);
    delete Work;
  } // end of 1D FFT

  /* 2D transform. */

  if (Rank == 2) {
    float *Work = new float[2*(Dimension[0]+15 + 2*Dimension[1]+15)];
    cfft2di(Dimension[0], Dimension[1], Work);
    cfft2d(direction*-1, Dimension[0], Dimension[1], buffer, 
	   Dimension[0], Work);
    if (direction == FFT_INVERSE)          // scale by 1/n for inverse FFT
      cscal2d(Dimension[0], Dimension[1], 1.0/float(Dimension[0]+Dimension[1]),
	      buffer, Dimension[0]);
    delete Work;
  } // end of 2D FFT

  /* 3D transform. */

  if (Rank == 3) {
    float *Work = new float[2*(  Dimension[0]+15 + 2*Dimension[1]+15 +
			       2*Dimension[2]+15)];
    cfft3di(Dimension[0], Dimension[1], Dimension[2], Work);
    cfft3d(direction*-1, Dimension[0], Dimension[1], Dimension[2], buffer, 
	    Dimension[0], Dimension[1], Work);
    if (direction == FFT_INVERSE)          // scale by 1/n for inverse FFT
      cscal3d(Dimension[0], Dimension[1], Dimension[2], 
	      1.0/float(Dimension[0]+Dimension[1]+Dimension[2]),
	      buffer, Dimension[0], Dimension[1]);
    delete Work;
  } // end of 3D FFT

#endif /* IRIS4 && SGI_MATH_COMPLEX */

  return SUCCESS;
}

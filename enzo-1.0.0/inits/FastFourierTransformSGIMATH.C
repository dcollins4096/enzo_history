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
/  COMPUTE A FAST FOURIER TRANSFORM USING THE SGI_MATH LIBRARY 
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

/* Defines */

#if defined(IRIS4) && defined(SGI_MATH)
extern "C" float *sfft1di(int dim1, float *Work);
extern "C" int    sfft1du(int dir, int dim1, float *data, int skip, 
		 	  float *Work);
extern "C" void   sscal1d(int dim1, float scale, float *data, int skip);

extern "C" float *sfft2dui(int dim1, int dim2, float *Work);
extern "C" int    sfft2du(int dir, int dim1, int dim2, float *data, 
			  int ld1, float *Work);
extern "C" void   sscal2d(int dim1, int dim2, float scale, float *data, 
			  int ld1);

extern "C" float *sfft3dui(int dim1, int dim2, int dim3, float *Work);
extern "C" int    sfft3du(int dir, int dim1, int dim2, int dim3, float *data, 
			  int ld1, int ld2, float *Work);
extern "C" void   sscal3d(int dim1, int dim2, int dim3, float scale, 
			  float *data, int ld1, int ld2);
#define GOT_FFT
#endif /* IRIS4 && SGI_MATH */

/* Start routine. */

int FastFourierTransformSGIMATH(float *buffer, int Rank, int DimensionReal[], 
				int Dimension[], int direction)
{

#if defined(IRIS4) && defined(SGI_MATH)

  /* Use SGI's routines.  Note: the direction flag is defined backwards to the
     usual sense (duh). */

  /* 1D transform. */

  if (Rank == 1) {
    float *Work = new float[2*(DimensionReal[0]+15)];
    sfft1di(Dimension[0], Work);
    sfft1du(direction*-1, Dimension[0], buffer, 1, Work);
    if (direction == FFT_INVERSE)          // scale by 1/n for inverse FFT
      sscal1d(Dimension[0], 1.0/float(Dimension[0]), buffer, 1);
    delete Work;
  } // end of 1D FFT

  /* 2D transform. */

  if (Rank == 2) {
    float *Work = new float[2*(DimensionReal[0]+15 + 2*DimensionReal[1]+15)];
    sfft2dui(Dimension[0], Dimension[1], Work);
    sfft2du(direction*-1, Dimension[0], Dimension[1], buffer, 
	    DimensionReal[0], Work);
    if (direction == FFT_INVERSE)          // scale by 1/n for inverse FFT
      sscal2d(Dimension[0], Dimension[1], 1.0/float(Dimension[0]*Dimension[1]),
	      buffer, DimensionReal[0]);
    delete Work;
  } // end of 2D FFT

  /* 3D transform. */

  if (Rank == 3) {
    float *Work = new float[2*(  DimensionReal[0]+15 + 2*DimensionReal[1]+15 +
			       2*DimensionReal[2]+15)];
    sfft3dui(Dimension[0], Dimension[1], Dimension[2], Work);
    sfft3du(direction*-1, Dimension[0], Dimension[1], Dimension[2], buffer, 
	    DimensionReal[0], DimensionReal[1], Work);
    if (direction == FFT_INVERSE)          // scale by 1/n for inverse FFT
      sscal3d(Dimension[0], Dimension[1], Dimension[2], 
	      1.0/float(Dimension[0]*Dimension[1]*Dimension[2]),
	      buffer, DimensionReal[0], DimensionReal[1]);
    delete Work;
  } // end of 3D FFT

  return SUCCESS;

#else /* IRIS4 && SGI_MATH */

  /* This is an error. */

  fprintf(stderr, "What are we doing here!?!\n");
  return FAIL;

#endif /* IRIS4 && SGI_MATH */

}

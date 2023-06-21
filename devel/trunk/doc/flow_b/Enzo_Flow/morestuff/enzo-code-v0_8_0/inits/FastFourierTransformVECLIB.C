/***********************************************************************
/
/  COMPUTE A FAST FOURIER TRANSFORM USING THE CONVEX VECLIB LIBRARY 
/
/  written by: Greg Bryan
/  date:       June, 1995
/  modified1:
/
/  PURPOSE:  This routine is an interface to the FORTRAN-callable CONVEX
/            Real-to-complex 1/2/3 dimension FFT routines.
/            NOTE: For some reason, the 2/3 D routines are quite different
/                  from the 1D routine.  The 1D routines requires a seperate
/                  initialization call and only supports 2^n lengths, while
/                  the 2 and 3 routines don't need a seperate call.
/            For this reason, we'll just fake calls to the 3D routine even
/            if it's 1 or 2 D.
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

void VECLIBInterpretError3D(int ErrorNumber, int l1, int l2, int l3,
			     int ldx, int mdx);

#if (defined(CONVEX) || defined(SPP)) && defined(VECLIB)

extern "C" void FORTRAN_NAME(src3ft)(float *data, int *l1, int *l2, int *l3, 
				     int *ldx, int *mdx, int *iopt, int *ier);

#define GOT_FFT
#endif /* (CONVEX || SPP) && VECLIB */

/* Start routine. */

int FastFourierTransformVECLIB(float *buffer, int Rank, int DimensionReal[], 
				int Dimension[], int direction)
{

#if (defined(CONVEX) || defined(SPP)) && defined(VECLIB)

  /* declarations. */

  int iopt, ier, l1, ldx, l2 = 1, l3 = 1, mdx = 1;

  /* Error check. */

  if (Rank < 1 || Rank > 3) {
    fprintf(stderr, "Rank = %d not supported.\n", Rank);
    return FAIL;
  }

  /* Setup for 1/2/3D transform. */

  l1  = Dimension[0];
  ldx = DimensionReal[0];

  /* Setup for 2/3D transform. */

  if (Rank >= 2) {
    l2 = Dimension[1];
    mdx = DimensionReal[1];
  }

  /* Setup for 3D transform. */

  if (Rank >= 3)
    l3 = Dimension[2];

  /* Set iopt according to direction of transform. */

  if (direction == FFT_FORWARD) iopt = 0;
  if (direction == FFT_INVERSE) iopt = -1;

  /* 1/2/3D transform. */

  if (Rank >= 1 && Rank <= 3) {

    /* Perform transform (and scale if necessary). */

    FORTRAN_NAME(src3ft)(buffer, &l1, &l2, &l3, &ldx, &mdx, &iopt, &ier);

    /* Error handling. */

    if (ier != 0) {
      VECLIBInterpretError3D(ier, l1, l2, l3, ldx, mdx);
      return FAIL;
    }
    
  } // end of 1/2/3D FFT

  return SUCCESS;

#else /* (CONVEX || SPP) && VECLIB */

  /* This is an error. */

  fprintf(stderr, "What are we doing here!?!\n");
  return FAIL;

#endif /* (CONVEX || SPP) && VECLIB */

}


void VECLIBInterpretError3D(int ErrorNumber, int l1, int l2, int l3,
			     int ldx, int mdx)
{
  fprintf(stderr, "VECLIB: l1 = %d, l2 = %d, l3 = %d, ldx = %d, mdx = %d.\n",
	  l1, l2, l3, ldx, mdx);

  switch (ErrorNumber) {

  case -1:
    fprintf(stderr, "VECLIB: l1 not of the required form.\n");
    break;

  case -2:
    fprintf(stderr, "VECLIB: l2 not of the required form.\n");
    break;

  case -3:
    fprintf(stderr, "VECLIB: l3 not of the required form.\n");
    break;

  case -4:
    fprintf(stderr, "VECLIB: ldx < l1+2.\n");
    break;

  case -5:
    fprintf(stderr, "VECLIB: mdx < l2.\n");
    break;

  case -6:
    fprintf(stderr, "VECLIB: Probable error in ldx or mdx\n");
    break;

  default:
    fprintf(stderr, "VECLIB: unknown error = %d.\n", ErrorNumber);
    break;

  }
}

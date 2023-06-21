/***********************************************************************
/
/  COMPUTE A FAST FOURIER TRANSFORM USING THE FFTW LIBRARY
/
/  written by: Elizabeth Tasker
/  date:       April, 2005
/  modified1:
/
/  PURPOSE:
/     This routine calls the C FFTW routines. 
/     1/2/3 dimensions are supported with real-to-complex and
/     complex-to-complex transforms.
/
/  INPUTS:
/     buffer - field to be FFTed
/     Rank   - rank of FFT
/     DimensionReal[] - declared dimensions of buffer
/     Dimension[]     - active dimensions of buffer
/     dir             - +1 forward, -1 inverse
/     type            - REAL_TO_COMPLEX or COMPLEX_TO_COMPLEX
/
/ 
/  COMMENTS:
/     FFTW requires row-ordered, not column-ordered input so dimensions 
/     reversed in routine calls.
/     Forwards/Backwards FT oppositely defined, so -1 is forward.
/
***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"

#if defined(FFTW)
#include "fftw3.h"
#endif /* FFTW */

int FastFourierTransformFFTW(float *buffer, int Rank, int DimensionReal[], 
			     int Dimension[], int dir, int type)
{

#if defined(FFTW)

  fftwf_plan plan;
  fftwf_complex *out_complex = (fftwf_complex*) buffer;
  fftwf_complex *in_complex = (fftwf_complex*) buffer;
  float *out_real = buffer;
  float *in_real = buffer;
  int i, size, rsize;

  /* 1D transform */

  if (Rank == 1){
    size = Dimension[0];
    rsize = DimensionReal[0];
    if (type == REAL_TO_COMPLEX) {
      if (dir == FFT_FORWARD){
	plan =  fftwf_plan_dft_r2c_1d(Dimension[0], in_real, out_complex, FFTW_ESTIMATE);
      }
      else {
	plan = fftwf_plan_dft_c2r_1d(Dimension[0], in_complex, out_real, FFTW_ESTIMATE);
      }
    }
    else {
      plan = fftwf_plan_dft_1d(Dimension[0], in_complex, out_complex, (-1)*dir, FFTW_ESTIMATE);
    }
  } // end of 1D FFT

  /* 2D transform */

  if (Rank == 2){
    size = Dimension[0]*Dimension[1];
    rsize = DimensionReal[0]*DimensionReal[1];
    if (type == REAL_TO_COMPLEX) {
      if (2*(Dimension[0]/2+1) != DimensionReal[0]) {
	fprintf(stderr, "FFTW array size incorrect.\n");
	return FAIL;
      }
      if (dir == FFT_FORWARD){
	plan = fftwf_plan_dft_r2c_2d(Dimension[1], Dimension[0], in_real, out_complex, FFTW_ESTIMATE);
      }
      else {
	plan = fftwf_plan_dft_c2r_2d(Dimension[1], Dimension[0], in_complex, out_real, FFTW_ESTIMATE);
      }
    }
    else {
      plan = fftwf_plan_dft_2d(Dimension[1], Dimension[0], in_complex, out_complex, (-1)*dir, FFTW_ESTIMATE);
    }
  } // end of 2D FFT

  if (Rank == 3) {
    size = Dimension[0]*Dimension[1]*Dimension[2];
    rsize = DimensionReal[0]*DimensionReal[1]*DimensionReal[2];
    if (type == REAL_TO_COMPLEX) {
      if (2*(Dimension[0]/2+1) != DimensionReal[0] ||
	  Dimension[1] != DimensionReal[1]) {
	fprintf(stderr, "FFTW(3d) array size incorrect.\n");
	return FAIL;
      }
      if (dir == FFT_FORWARD){
        plan = fftwf_plan_dft_r2c_3d(Dimension[2], Dimension[1], Dimension[0], in_real, out_complex, FFTW_ESTIMATE);
      }
      else {
        plan = fftwf_plan_dft_c2r_3d(Dimension[2], Dimension[1], Dimension[0], in_complex, out_real, FFTW_ESTIMATE);
      }
    }
    else {
      plan = fftwf_plan_dft_3d(Dimension[2], Dimension[1], Dimension[0], in_complex, out_complex, (-1)*dir, FFTW_ESTIMATE);
    }
  } // end of 3d FFT

  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  // Normalisation (use proper FFT size to calculate renormalization factor 
  //  but make sure to loop over the actual array size; should be the same
  //  for the complex FFT but may differ in the leading dimension in r2c).

  float factor = 1.0/size;
  if (dir == FFT_INVERSE) {
    if (type == REAL_TO_COMPLEX) {
      for (i = 0; i < rsize; i++)
	out_real[i] *= factor;
    }
    else {
      for (i = 0; i < rsize; i++) {
	out_complex[i][0] *= factor;
        out_complex[i][1] *= factor;

      }
    }
  }
  
#endif /* FFTW */

  return SUCCESS;
}

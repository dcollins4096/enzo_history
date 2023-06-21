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
/  GENERATES A RANDOM FIELD REALIZATION
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "global_data.h"
#include "PowerSpectrumParameters.h"
#include "CosmologyParameters.h"

/* function prototypes */

extern "C" void FORTRAN_NAME(make_field_kpreserving)
             (float *field, int *nx, int *ny, int *nz,
	      int *in, int *jn, int *kn, int *itype, int *iseed, float *box,
	      float *PSTable, float *PSMin, float *PSStep, int *kfcutoff);
extern "C" void FORTRAN_NAME(make_field)
             (float *field, int *nx, int *ny, int *nz,
	      int *nxmax, int *nymax, int *nzmax,
	      int *in, int *jn, int *kn, int *itype, int *iseed, float *box,
	      float *PSTable, float *PSMin, float *PSStep, int *kfcutoff);
#ifdef SHIFT_FOR_LARS
extern "C" void FORTRAN_NAME(shift)
             (float *field, int *nx, int *ny, int *nz,
	      int *in, int *jn, int *kn, float *temp);
#endif /* SHIFT_FOR_LARS */
int FastFourierTransform(float *buffer, int Rank, int DimensionReal[], 
			 int Dimension[], int direction, int type);


int GenerateField(int Rank, int Dims[3], int MaxDims[3], int WaveNumberCutoff,
		  float *Field, int FieldType, int NewCenter[3], 
		  int Refinement, int StartIndex[3], int Species)
{

  float *Temp, k1, k2, delk, box;
  int dim, i, j, k, ii, jj, kk, Extract = FALSE, size = 1, findex, tindex;
  int RealDims[3] = {1,1,1}, TempDims[3] = {1,1,1}, MoveBy[3] = {0,0,0};

  /* If Field doesn't correspond to the whole volume, we must allocate a 
     temporary field which does. */

  for (dim = 0; dim < Rank; dim++) {
    if (StartIndex[dim] != 0 || Dims[dim]*Refinement != MaxDims[dim])
      Extract = TRUE;
    TempDims[dim] = MaxDims[dim]/Refinement;
    RealDims[dim] = TempDims[dim] + ((dim == 0) ? 2 : 0);
    size *= RealDims[dim];
  }
  if ((Temp = new float[size]) == NULL) {
    fprintf(stderr, "GenerateField: malloc failure (%d).\n", size);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i<size; i++) {
    Temp[i] = 0.0;
  }

  /* Fill in Temp with random phase complex values. */

  if (debug) printf("filling k-space...\n");
  k1 = log(kmin);
  k2 = log(kmax);
  delk = (k2 - k1)/(NumberOfkPoints-1);
  box = ComovingBoxSize/HubbleConstantNow;
  if (Rank == 3)
    FORTRAN_NAME(make_field_kpreserving)(Temp, TempDims, TempDims+1,TempDims+2, 
			   RealDims, RealDims+1, RealDims+2, &FieldType,
			      &RandomSeed, &box,
			   PSLookUpTable[Species], &k1, &delk, 
			      &WaveNumberCutoff);
  else
    FORTRAN_NAME(make_field)(Temp, TempDims, TempDims+1, TempDims+2, 
			     MaxDims, MaxDims+1, MaxDims+2,
			   RealDims, RealDims+1, RealDims+2, &FieldType,
			      &RandomSeed, &box,
			   PSLookUpTable[Species], &k1, &delk, 
			      &WaveNumberCutoff);

  /* Perform an inverse FFT. */

  if (debug) printf("transforming...\n");
  if (FastFourierTransform(Temp, Rank, RealDims, TempDims, FFT_INVERSE,
			   REAL_TO_COMPLEX)  == FAIL) {
    fprintf(stderr, "GenerateField: FFT error.\n");
    exit(EXIT_FAILURE);
  }

  /* Shift by 1/2 a cell. */

#ifdef SHIFT_FOR_LARS
  if (debug) printf("shifting...\n");
  float *Temp2 = new float[size];
  FORTRAN_NAME(shift)(Temp, TempDims, TempDims+1, TempDims+2, 
		      RealDims, RealDims+1, RealDims+2, Temp2);
  delete Temp2;
#endif /* SHIFT_FOR_LARS */

  if (Extract) {

    /* Extract the requested volume (if necessary). */

    if (debug) printf("extracting...\n");
    for (dim = 0; dim < Rank; dim++)
      if (StartIndex[dim] % Refinement != 0) {
	fprintf(stderr, "Extract: StartIndex[%d] = %d must be divisible by refinement.\n", dim, StartIndex[dim]);
	exit(EXIT_FAILURE);
      }
    for (k = 0; k < Dims[2]; k++) {
      kk = k+StartIndex[2]/Refinement;
      if (kk >= MaxDims[2]/Refinement) kk -= MaxDims[2]/Refinement;
      for (j = 0; j < Dims[1]; j++) {
	jj = j+StartIndex[1]/Refinement;
	if (jj >= MaxDims[1]/Refinement) jj -= MaxDims[1]/Refinement;
	tindex = (kk*RealDims[1] + jj)*RealDims[0];
	findex = (k*Dims[1] + j)*Dims[0];
	for (i = 0; i < Dims[0]; i++) {
	  ii = i+StartIndex[0]/Refinement;
	  if (ii >= MaxDims[0]/Refinement) ii -= MaxDims[0]/Refinement;
	  Field[findex+i] = Temp[tindex+ii];
	}
      }
    }

  } else {

    /* Recenter the volume if requested (otherwise just compactify). */

    for (dim = 0; dim < Rank; dim++) {
      if (NewCenter[dim] == INT_UNDEFINED)
	MoveBy[dim] = 0;
      else
	MoveBy[dim] = NewCenter[dim] - (MaxDims[dim]/2 - 1);
      if (abs(MoveBy[dim]) % Refinement != 0) {
	fprintf(stderr, "Centering: move by %d %d %d must be divisible by refinement.\n", MoveBy[0], MoveBy[1], MoveBy[2]);
	exit(EXIT_FAILURE);
      }
      MoveBy[dim] /= Refinement;
    }

    if (debug) 
      printf("Centering: move by %d %d %d\n", MoveBy[0], MoveBy[1], MoveBy[2]);
    for (k = 0; k < Dims[2]; k++) {
      kk = (k + MoveBy[2] + Dims[2]) % Dims[2];
      for (j = 0; j < Dims[1]; j++) {
	jj = (j + MoveBy[1] + Dims[1]) % Dims[1];
	findex = (k*Dims[1] + j)*(Dims[0]+0);
	tindex = (kk*Dims[1] + jj)*(Dims[0]+2);
	ii = (0 + MoveBy[0] + Dims[0]) % Dims[0];
	for (i = 0; i < Dims[0]; i++)
	  //	  Field[tindex+((ii+i) % Dims[0])] = Temp[findex + i];
	  Field[findex + i] = Temp[tindex+((ii+i) % Dims[0])];
      }
    }

  } /* end: if (Extract) */

  /* deallocate temporary. */

  delete Temp;

  return SUCCESS;
}

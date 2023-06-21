/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Self-Gravity Problem Class setup routine
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: This routine is used to set isolating boundary values for the 
/           self-gravity potential field solve on any dimension {0,1,2}, 
/           and for either low or high face {0,1}.
/
/           For setting the face values, either scalar or array-valued 
/           arguments may be given to set the values on that face, 
/           depending on the value of BdryConst (zero implies 
/           array-valued, nonzero implies constant).
/
/           If a scalar is desired, it will be taken as the value referenced
/           by the pointer BdryValue (pass the constant by reference).
/
/           If array input is chosen, the BdryValue array will 
/           be copied in column-major (Fortran-style) ordering:
/
/                   dim  |  fast index  |  slow index
/                 -------------------------------------
/                    0   |      x1      |      x2
/                    1   |      x2      |      x0
/                    2   |      x0      |      x1
/                 -------------------------------------
/
/           These boundary values for the potential field should be 
/           given over the full extent of the gravitational potential 
/           boundary, i.e. they must extend one cell further out at 
/           both the left and right boundaries than the hydrodynamics 
/           domain (including hydro ghost zones).
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
#ifdef ISO_GRAV
#include "SelfGravityProblem_preincludes.h"
#include "SelfGravityProblem.h"



int SelfGravityProblem::SetupBoundary(int Dim, int Face, 
				      int BdryConst, float *BdryData) 
{
  // local variables
  int size, facesize, index, i;

  // Error check
  if ((Dim < 0) || (Dim >= rank)) {
    fprintf(stderr, "SetupBoundary: Dim %"ISYM" out of bounds.\n", Dim);
    return FAIL;
  }
  if ((Face != 0) && (Face != 1)) {
    fprintf(stderr, "SetupBoundary: Face %"ISYM" != {0,1}.\n", Face);
    return FAIL;
  }

  // compute size of local mesh and relevant faces
  size = 1;
  for (i=0; i<rank; i++)  size *= LocDims[i];
  facesize = size/LocDims[Dim];

  // delete previous boundary values arrays
  if (BdryVals[Dim][Face] != NULL)
    delete[] BdryVals[Dim][Face];
  
  // initialize and set boundary conditions
  BdryVals[Dim][Face] = new float[facesize];
  
  // set constant face conditions
  if (BdryConst) {
    for (index=0; index<facesize; index++) 
      BdryVals[Dim][Face][index] = *BdryData;
  }
  
  // set spatially-varying face conditions
  else {
    for (index=0; index<facesize; index++) 
      BdryVals[Dim][Face][index] = BdryData[index];
  }

  return SUCCESS;
}
#endif

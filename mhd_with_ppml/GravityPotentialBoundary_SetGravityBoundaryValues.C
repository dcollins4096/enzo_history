/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Daniel R. Reynolds                                         *
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRAVITY POTENTIAL BOUNDARY CLASS (SET BOUNDARY FACE VALUES)
/
/  written by: Daniel R. Reynolds
/  date:       September, 2005
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

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "GravityPotentialBoundary.h"
#include "ExternalBoundary.h"
#include "Grid.h"


#ifdef ISO_GRAV
int GravityPotentialBoundary::SetGravityBoundaryValues(int Dim, int Face, 
						       int BdryConst, 
						       float *BdryValues)
{
  int size, facesize, index, i;

  /* Error check */
  if (Dim > BoundaryRank) {
    fprintf(stderr, "SetGravityBoundaryValues: Dim %"ISYM" > Rank %"ISYM".\n", 
	    Dim, BoundaryRank);
    return FAIL;
  }
  if ((Face != 0) && (Face != 1)) {
    fprintf(stderr, "SetGravityBoundaryValues: Face %"ISYM" != {0,1}.\n", Face);
    return FAIL;
  }

  /* compute size of local mesh and relevant faces */
  size = 1;
  for (i=0; i<BoundaryRank; i++) size *= LocalDimension[i];
  facesize = size/LocalDimension[Dim];
//   printf("GravityPotentialBoundary::SetGravityBoundaryValues dim %"ISYM", face %"ISYM", local mesh is %"ISYM", face size is %"ISYM"\n", Dim, Face, size, facesize);
  
  /* initialize and set boundary conditions */
  BoundaryValues[Dim][Face] = new float[facesize];
  
  /* set constant face conditions */
  if (BdryConst) {
    for (index=0; index<facesize; index++) 
      BoundaryValues[Dim][Face][index] = *BdryValues;
  }
  
  /* set spatially-varying face conditions */
  else {
    for (index=0; index<facesize; index++) 
      BoundaryValues[Dim][Face][index] = BdryValues[index];
  }
  
  return SUCCESS;

}
#endif

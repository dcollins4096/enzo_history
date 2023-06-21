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
/  EXTERNAL BOUNDARY CLASS (SET A GRID'S BOUNDARY)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1: Hao Xu
/   date:       July,2006
/          Setting boundary of Electronic fields, Using it to update 
/          magnetic fields, so preserve the divergence-free   
/
/  PURPOSE:
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
#include "ExternalBoundary.h"
#include "Grid.h"

// This is used to set the corners (which are not really used) of the
//   grid to something reasonable in the case of periodic B.C.'s

#define USE_PERIODIC

#ifdef EMF_BOUNDARY

int ExternalBoundary::SetExternalBoundary_EMF(int FieldRank, int GridDims[],
					  int GridOffset[],
					  int StartIndex[], int EndIndex[],
					  float *Field, int direction)
{

  /* declarations */

  int i, j, k, dim, Sign, bindex, ver = FALSE, EMFEndIndex[3], EMFGridDims[3];
  float *index;

  switch(direction){
  case 0: 
     EMFEndIndex[0] = EndIndex[0];
     EMFEndIndex[1] = EndIndex[1]+1;
     EMFEndIndex[2] = EndIndex[2]+1;
     EMFGridDims[0] = GridDims[0];
     EMFGridDims[1] = GridDims[1]+1;
     EMFGridDims[2] = GridDims[2]+1;
   break;
  case 1:
     EMFEndIndex[0] = EndIndex[0]+1;
     EMFEndIndex[1] = EndIndex[1];
     EMFEndIndex[2] = EndIndex[2]+1;
     EMFGridDims[0] = GridDims[0]+1;
     EMFGridDims[1] = GridDims[1];
     EMFGridDims[2] = GridDims[2]+1;
   break;
   case 2:
     EMFEndIndex[0] = EndIndex[0]+1;
     EMFEndIndex[1] = EndIndex[1]+1;
     EMFEndIndex[2] = EndIndex[2];
     EMFGridDims[0] = GridDims[0]+1;
     EMFGridDims[1] = GridDims[1]+1;
     EMFGridDims[2] = GridDims[2];
   break;
  default:
   break;
  }

  /* error check: grid ranks */

  if (FieldRank != BoundaryRank) {
    fprintf(stderr, "FieldRank(%d) != BoundaryRank(%d).\n",
            FieldRank, BoundaryRank);
    return FAIL;
  }

int field;
   
   field=direction;

  /* error check: make sure the boundary type array exists */

  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] != 1) {
      if (MagneticBoundaryType[field][dim][0] == NULL) {
        fprintf(stderr, "BoundaryType not yet declared.\n");
        return FAIL;
      }
    }

  

  /* set Boundary conditions */

  Sign = 1;
  
  if (BoundaryDimension[0] > 1 && GridOffset[0] == 0) {

    /* set x inner (left) face */

    for (i = 0; i < StartIndex[0]; i++)
      for (j = 0; j < EMFGridDims[1]; j++)
	for (k = 0; k < EMFGridDims[2]; k++) {
	  index = Field + i + j*EMFGridDims[0] + k*EMFGridDims[1]*EMFGridDims[0];
	  bindex = j+GridOffset[1] + (k+GridOffset[2])*BoundaryDimension[1];
	  switch (MagneticBoundaryType[field][0][0]) {
	  case reflecting:
	    *index = Sign*(*(index + (2*StartIndex[0] - 1 - 2*i))); 
	    break;
	  case outflow:
	    *index =       *(index + (  StartIndex[0]     -   i)) ; 
	    break;
	  case inflow:
	    *index = 0.0; //EMFBoundaryValue[field][0][0][bindex];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index + (EMFEndIndex[0] - StartIndex[0] + 1)); 
#endif /* USE_PERIODIC */
	    break;
	  case BoundaryUndefined:
	  default:
	    fprintf(stderr, "BoundaryType not recognized (x-left).\n");
	    return FAIL;
	  }
	}
  }
  
  if (BoundaryDimension[0] > 1 && GridOffset[0]+GridDims[0] == BoundaryDimension[0]) {

    /* set x outer (right) face */

    for (i = 0; i < EMFGridDims[0]-EMFEndIndex[0]-1; i++)
      for (j = 0; j < EMFGridDims[1]; j++)
	for (k = 0; k < EMFGridDims[2]; k++) {
	  index = Field + i + EMFEndIndex[0]+1 + 
	    j*EMFGridDims[0] + k*EMFGridDims[1]*EMFGridDims[0];
	  bindex = j+GridOffset[1] + (k+GridOffset[2])*BoundaryDimension[1];
	  switch (MagneticBoundaryType[field][0][1]) {
	  case reflecting:
	    *index = Sign*(*(index - (2*i + 1))); 
	    break;
	  case outflow:
	    *index =       *(index + (-1 - i)) ; 
	    break;
	  case inflow:
	    *index = 0.0;//EMFBoundaryValue[field][0][1][bindex];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index - (EMFEndIndex[0] - StartIndex[0] + 1)); 
#endif /* USE_PERIODIC */
	    break;
	  case BoundaryUndefined:
	  default:
	    fprintf(stderr, "BoundaryType not recognized (x-right).\n");
	    return FAIL;
	  }
	}							 


  }
  
  /* set y inner (left) face */

  Sign = 1;
  
  if (BoundaryDimension[1] > 1 && GridOffset[1] == 0) {

    for (j = 0; j < StartIndex[1]; j++)
      for (i = 0; i < EMFGridDims[0]; i++)
	for (k = 0; k < EMFGridDims[2]; k++) {
	  index = Field + i + j*EMFGridDims[0] + k*EMFGridDims[1]*EMFGridDims[0];
	  bindex = i+GridOffset[0] + (k+GridOffset[2])*BoundaryDimension[0];
	  switch (MagneticBoundaryType[field][1][0]) {
	  case reflecting:
	    *index = Sign*(*(index + (2*StartIndex[1] - 1 - 2*j)*EMFGridDims[0])); 
	    break;
	  case outflow:
	    *index =       *(index + (  StartIndex[1]     - j)*EMFGridDims[0]) ;
	    break;
	  case inflow:
	    *index = 0.0; //EMFBoundaryValue[field][1][0];
	     break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index + (EMFEndIndex[1] - StartIndex[1] + 1)*EMFGridDims[0]);
#endif /* USE_PERIODIC */
	     break;
	  case BoundaryUndefined:
	  default:
	    fprintf(stderr, "BoundaryType not recognized (y-left).\n");
	    return FAIL;
	  }
	}


  }

  if (BoundaryDimension[1] > 1 && GridOffset[1]+GridDims[1] == BoundaryDimension[1]) {

    /* set y outer (right) face */

    for (j = 0; j < EMFGridDims[1]-EMFEndIndex[1]-1; j++)
      for (i = 0; i < EMFGridDims[0]; i++)
	for (k = 0; k < EMFGridDims[2]; k++) {
	  index = Field + i + (j + EMFEndIndex[1]+1)*EMFGridDims[0] + 
	    k*EMFGridDims[1]*EMFGridDims[0];
	  bindex = i+GridOffset[0] + (k+GridOffset[2])*BoundaryDimension[0];
	  switch (MagneticBoundaryType[field][1][1]) {
	  case reflecting:
	    *index = Sign*(*(index - (2*j + 1)*EMFGridDims[0])); 
	    break;
	  case outflow:
	    *index =       *(index + (-1 - j)*EMFGridDims[0]) ; 
	    break;
	  case inflow:
	    *index = 0.0;//EMFBoundaryValue[field][1][1]; 
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index - (EMFEndIndex[1] - StartIndex[1] + 1)*EMFGridDims[0]);
#endif /* USE_PERIODIC */
	    break;
	  case BoundaryUndefined:
	  default:
	    fprintf(stderr, "BoundaryType not recognized (y-right).\n");
	    return FAIL;
	  }
	}							 


  }
  
  /* set z inner (left) face */

  Sign = 1;
  
  if (BoundaryDimension[2] > 1 && GridOffset[2] == 0) {

    for (k = 0; k < StartIndex[2]; k++)
      for (i = 0; i < EMFGridDims[0]; i++)
	for (j = 0; j < EMFGridDims[1]; j++) {
	  index = Field + i + j*EMFGridDims[0] + k*EMFGridDims[1]*EMFGridDims[0];
	  bindex = i+GridOffset[0] + (j+GridOffset[1])*BoundaryDimension[0];
	  switch (MagneticBoundaryType[field][2][0]) {
	  case reflecting:
	    *index = Sign*(*(index + (2*StartIndex[2]-1 - 2*k)*EMFGridDims[0]*EMFGridDims[1]));
	    break;
	  case outflow:
	    *index =       *(index + (  StartIndex[2]   - k)*EMFGridDims[0]*EMFGridDims[1]) ;
	    break;
	  case inflow:
	    *index = 0.0;//EMFBoundaryValue[field][2][0];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index + (EMFEndIndex[2]-StartIndex[2]+1)*EMFGridDims[0]*EMFGridDims[1]);
#endif /* USE_PERIODIC */
	    break;
	  case BoundaryUndefined:
	  default:
	    fprintf(stderr, "BoundaryType not recognized (z-left).\n");
	    return FAIL;
	  }
	}

  }

  if (BoundaryDimension[2] > 1 && GridOffset[2]+GridDims[2] == BoundaryDimension[2]) {

    /* set z outer (right) face */

    for (k = 0; k < EMFGridDims[2]-EMFEndIndex[2]-1; k++)
      for (i = 0; i < EMFGridDims[0]; i++)
	for (j = 0; j < EMFGridDims[1]; j++) {
	  index = Field + i + j*EMFGridDims[0] + 
	    (k + EMFEndIndex[2]+1)*EMFGridDims[1]*EMFGridDims[0];
	  bindex = i+GridOffset[0] + (j+GridOffset[1])*BoundaryDimension[0];
	  switch (MagneticBoundaryType[field][2][1]) {
	  case reflecting:
	    *index = Sign*(*(index - (2*k + 1)*EMFGridDims[0]*EMFGridDims[1]));
	    break;
	  case outflow:
	    *index =       *(index + (-1 - k)*EMFGridDims[0]*EMFGridDims[1]) ;
	    break;
	  case inflow:
	    *index = 0.0;//EMFBoundaryValue[field][2][1];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index - (EMFEndIndex[2]-StartIndex[2]+1)*EMFGridDims[0]*EMFGridDims[1]);
#endif /* USE_PERIODIC */
	    break;
	  case BoundaryUndefined:
	  default:
	    fprintf(stderr, "BoundaryType not recognized (z-right).\n");
	    return FAIL;
	  }
	}							 
  }
  

  return SUCCESS;
}
#endif //EMF_BOUNDARY


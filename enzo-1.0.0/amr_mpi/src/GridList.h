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
/  GRID LIST STRUCTURE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"

struct GridList
{
  GridList *NextGrid;
  int GridRank;
  int GridDimension[MAX_DIMENSION];
  FLOAT GridLeftEdge[MAX_DIMENSION];
  FLOAT GridRightEdge[MAX_DIMENSION];

  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];
  int NumberFlagged;
};

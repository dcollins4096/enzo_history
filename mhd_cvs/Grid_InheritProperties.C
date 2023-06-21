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
/  GRID CLASS (INHERIT BASIC GRID PROPERTIES FROM THE PARENT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//
//  Assign basic values to a grid (allocate fields)
//

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void grid::InheritProperties(grid *ParentGrid)
{

  /*  Set rank and current time */

  GridRank = ParentGrid->GridRank;
  Time     = ParentGrid->Time;

  /*  Baryons only: set up field quantities and allocate fields
       (we assume here the grid is uniform in each dimension) */

  NumberOfBaryonFields = ParentGrid->NumberOfBaryonFields;

  for (int field = 0; field < NumberOfBaryonFields; field++)
    FieldType[field]      = ParentGrid->FieldType[field];

  CourantSafetyNumber    = ParentGrid->CourantSafetyNumber;
  PPMFlatteningParameter = ParentGrid->PPMFlatteningParameter;
  PPMDiffusionParameter  = ParentGrid->PPMDiffusionParameter;
  PPMSteepeningParameter = ParentGrid->PPMSteepeningParameter;

  /* For gravity, this must be a subgrid, so set the boundary accordingly. */

  GravityBoundaryType = SubGridIsolated;

}

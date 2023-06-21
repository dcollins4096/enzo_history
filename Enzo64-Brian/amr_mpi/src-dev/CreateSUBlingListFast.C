
/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Rick Wagner                                                *
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  CREATE SUBLINGLIST
/
/  written by: David Collins.
/  date:       September, 2005
/                
/
/  PURPOSE: Create a list, for each grid, of all grids on the next finer level
/           that share a face ONLY.  
/           Note that this configuration does NOT recognize grids that are strict subgrids
/           AT ALL, so if you're running with fewer than 2 processors in each direction
/           AND periodic boundary conditions, there will be grids at the domain edge that may
/           cause conservation errors (since they'll both share a face AND be subgrids.)
/           In that case, the user will need to do some work to a.) get those grids included
/           in this list AND b.) avoid correcting those cells multiple times.  Contact me
/           if you need help.
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"


int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
                                 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);


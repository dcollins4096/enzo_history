
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

#ifdef DC_OPT_SIBSUB
int CreateSUBlingList(TopGridData *MetaData, 
			  HierarchyEntry *Grids[],
			  int NumberOfGrids,
			  LevelHierarchyEntry ***SUBlingList)
{


  //For historical reasons, this is a two pass ordeal.  First the SiblingList is made matching ThisLevel with the
  //finer level.  Then that list is refined to only match faces.

  //First create chaining mesh for all subgrids 
  int grid, grid2;
  ChainingMeshStructure SubgridMesh;
  HierarchyEntry * NextSubgrid;
  SiblingGridList * PotentialSUBlingList = new SiblingGridList[NumberOfGrids];


  FastSiblingLocatorInitialize(&SubgridMesh, MetaData->TopGridRank,
                               MetaData->TopGridDims);

  //generate chaining mesh for subgrids. 
  for( grid=0;grid<NumberOfGrids;grid++){
    NextSubgrid = Grids[grid]->NextGridNextLevel;
    while( NextSubgrid != NULL ){
      NextSubgrid->GridData->FastSiblingLocatorAddGrid(&SubgridMesh);
      NextSubgrid = NextSubgrid->NextGridThisLevel;
    }
  }

  //Now match "siblings" on the finer level to this level.  Nieces and Nephews, really.
    
  for( grid=0; grid<NumberOfGrids;grid++){
    if( Grids[grid]->GridData->FastSiblingLocatorFindSiblings(
					   &SubgridMesh, &PotentialSUBlingList[grid],
					   MetaData->LeftFaceBoundaryCondition,
					   MetaData->RightFaceBoundaryCondition) == FAIL) {
      fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
      return FAIL;
    }
  }


  FastSiblingLocatorFinalize(&SubgridMesh);

  //Now re-examine that list to check for only shared faces.  
  
  for (grid = 0; grid < NumberOfGrids; grid++){
    LevelHierarchyEntry *LastEntry = NULL;
    
    
    for (grid2 = 0; grid2 < PotentialSUBlingList[grid].NumberOfSiblings; grid2++){

      if( Grids[grid]->GridData->CheckForSharedFace( PotentialSUBlingList[grid].GridList[grid2],
						     MetaData->LeftFaceBoundaryCondition,
						     MetaData->RightFaceBoundaryCondition
						     ) == TRUE){
	  
	if( LastEntry != NULL ){
	  LastEntry->NextGridThisLevel = new LevelHierarchyEntry;
	  LastEntry = LastEntry->NextGridThisLevel;
	}else{
	  (*SUBlingList)[grid] = new LevelHierarchyEntry;
	  LastEntry = (*SUBlingList)[grid];
	}
	LastEntry->GridHierarchyEntry = NULL; //Don't have that info.
	LastEntry->GridData = PotentialSUBlingList[grid].GridList[grid2];
	LastEntry->NextGridThisLevel = NULL;
	
      }// if grids overlap 
      
      
    }//grid2
  } // loop over other grids this level 
  
  return SUCCESS;
  
}
#endif 

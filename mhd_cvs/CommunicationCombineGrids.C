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
/  COMMUNICATION ROUTINE: COMBINE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
//#include "pout.h"
/* function prototypes */


int CommunicationCombineGrids(HierarchyEntry *OldHierarchy,
			      HierarchyEntry **NewHierarchyPointer,
			      FLOAT WriteTime)
{

  /* If there is only one proc, then just point the new one at the old one. */

  //Pout( " CommbineGrids (PRGIO, TRUE)", ParallelRootGridIO, TRUE);

  if (NumberOfProcessors == 1 || ParallelRootGridIO == TRUE) {
    *NewHierarchyPointer = OldHierarchy;
    return SUCCESS;
  }

  /* Otherwise generate a new hierarchy entry and proceed. */

  HierarchyEntry *NewHierarchy = new HierarchyEntry;
  *NewHierarchyPointer = NewHierarchy;

  int Rank, dim, Dims[MAX_DIMENSION], NewDims[MAX_DIMENSION], 
      SendOffset[MAX_DIMENSION], TempDims[MAX_DIMENSION], 
      StartIndex[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION], CellSize[MAX_DIMENSION];

  /* Compute dims, etc. for new grid assuming it fills entire domain. */

  OldHierarchy->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (dim < Rank)
      Dims[dim] -= 2*DEFAULT_GHOST_ZONES;
    CellSize[dim] = (Right[dim] - Left[dim])/FLOAT(Dims[dim]);
    NewDims[dim] = nint((DomainRightEdge[dim] - DomainLeftEdge[dim])/
			CellSize[dim])
                   + ((dim < Rank) ? 2*DEFAULT_GHOST_ZONES : 0);
  }
  if (debug)
    printf("CombineGrids: NewDims = %d %d %d\n", 
	   NewDims[0], NewDims[1], NewDims[2]);

  /* Generate a new grid. */

  //Pout( "  Make new grid");

  NewHierarchy->GridData = new grid;
  NewHierarchy->GridData->InheritProperties(OldHierarchy->GridData);
  NewHierarchy->GridData->SetGravityParameters(
		       OldHierarchy->GridData->ReturnGravityBoundaryType());
  NewHierarchy->GridData->PrepareGrid(Rank, NewDims, DomainLeftEdge,
				      DomainRightEdge, 0);

  /* Loop over old grids and copy info. */

  HierarchyEntry *Temp = OldHierarchy;
  grid *NewGrid = NewHierarchy->GridData;
  //Pout( "  Gathering Grid Info");
  while (Temp != NULL) {
    /* Compute region to send. */

    grid *OldGrid = Temp->GridData;
    OldGrid->ReturnGridInfo(&Rank, TempDims, Left, Right);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      SendOffset[dim] = (dim < Rank)? DEFAULT_GHOST_ZONES : 0;
      TempDims[dim] -= 2*SendOffset[dim];
      StartIndex[dim] = nint((Left[dim] - DomainLeftEdge[dim])/CellSize[dim])
	              + SendOffset[dim];
    }

    //Sometimes E isn't created by the time this code is run. 

    if(MHD_Used == TRUE  && MyProcessorNumber == OldGrid->ReturnProcessorNumber() ){    
      for(int field=0;field<3;field++)
	if( OldGrid->ElectricField[field] == NULL ){
	  OldGrid->ElectricField[field] = new float[OldGrid->ElectricSize[field]];
	  for(int i=0;i<OldGrid->ElectricSize[field];i++)
	    OldGrid->ElectricField[field][i]=0.0;
	}
	  
    }

    /* Copy grid region. */
    
    int OldProc = OldGrid->ReturnProcessorNumber(),
        NewProc = NewGrid->ReturnProcessorNumber();
    //Pout( "  Receiving Region (NewGrid, OldGrid)", NewProc, OldProc);
//    printf("(%d): %d --> %d\n", MyProcessorNumber, OldProc, NewProc);
    if (MyProcessorNumber == NewProc || MyProcessorNumber == OldProc)
      if (NewGrid->CommunicationReceiveRegion(OldGrid, OldProc, ALL_FIELDS, 
			      ((WriteTime < 0) ? NEW_ONLY : NEW_AND_OLD),
			      StartIndex, TempDims, FALSE) == FAIL) {
	fprintf(stderr, "Error in grid->CommunicationReceiveRegion.\n");
	return FAIL;
      }

    /* Copy particles. */

//    if (MyProcessorNumber == NewProc || MyProcessorNumber == OldProc)
      if (OldGrid->CommunicationSendParticles(NewGrid, NewProc, 0, 
			    OldGrid->ReturnNumberOfParticles(), -1) == FAIL) {
	fprintf(stderr, "Error in grid->CommunicationSendParticles.\n");
	return FAIL;
      }

    /* Next Grid */

    Temp = Temp->NextGridThisLevel;
  }
//  printf("(%d): done\n", MyProcessorNumber);

  /* Create a new first level of hierarchy entries that are all below the
     new one.  Below that, just point back into the old hierarchy. */
  //Pout( "  Link list");

  NewHierarchy->ParentGrid = OldHierarchy->ParentGrid;
  NewHierarchy->NextGridNextLevel = NULL;
  NewHierarchy->NextGridThisLevel = NULL;
  Temp = OldHierarchy;
  HierarchyEntry *Previous = NULL;
  while (Temp != NULL) {

    HierarchyEntry *Temp2 = Temp->NextGridNextLevel;

    while (Temp2 != NULL) {

      HierarchyEntry *NewEntry = new HierarchyEntry;
      NewEntry->NextGridThisLevel = Previous;
      NewEntry->NextGridNextLevel = Temp2->NextGridNextLevel;
      NewEntry->ParentGrid        = NewHierarchy;
      NewEntry->GridData = Temp2->GridData;
      NewHierarchy->NextGridNextLevel = NewEntry;
      Previous = NewEntry;

      Temp2 = Temp2->NextGridThisLevel;
    }

    Temp = Temp->NextGridThisLevel;
  }
  //Pout( "  End of Combine Grids");

  return SUCCESS;
}


#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

int grid::MHD_SendOldFineGrids(LevelHierarchyEntry * OldFineLevel, grid *ParentGrid){

  int StartIndex[MAX_DIMENSION], EndIndex[MAX_DIMENSION], Offset[MAX_DIMENSION],TempDim[3];
  int Refinement[MAX_DIMENSION];
  int i,dim;

  ParentGrid->ComputeRefinementFactors(this, Refinement);

  for(dim=0;dim<GridRank;dim++){

    StartIndex[dim] = nint((CellLeftEdge[dim][0] - 
			    ParentGrid->CellLeftEdge[dim][0])/
			   CellWidth[dim][0]);
    EndIndex[dim] = nint((CellLeftEdge[dim][GridDimension[dim]-1] +
			  CellWidth[dim][GridDimension[dim]-1]
			  - ParentGrid->CellLeftEdge[dim][0])/
			 CellWidth[dim][0])
      - 1;
    
    /* Record the offset between the real StartIndex and the one adjusted
       to conform with the parent. */
    
    Offset[dim] = StartIndex[dim] -
      int(StartIndex[dim]/Refinement[dim])*Refinement[dim];
    
    /* Adjust Start/end index to conform with parental cell edges. */
    
    StartIndex[dim] = int(StartIndex[dim]/Refinement[dim])*Refinement[dim];
    EndIndex[dim] = (int(EndIndex[dim]/Refinement[dim])+1)*Refinement[dim]-1;
    
    TempDim[dim]            = EndIndex[dim] - StartIndex[dim] + 1;

  }

  CommunicationDirection = COMMUNICATION_SEND;
  //ComputeInterpolation Derivatives
  if( MHD_CID(OldFineLevel, Offset, TempDim, Refinement) == FAIL ){
    fprintf(stderr, "Shit!  MHD_CID failed.\n");
    return FAIL;
  }




  

  return SUCCESS;
}

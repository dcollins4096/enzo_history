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
/  GRID CLASS (COPY OVERLAPPING ZONES FROM GRID IN ARGUMENT TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

// This routine copies zones which overlap from the grid in the argument
//   to the current grid.  We use only the active region of the OtherGrid,
//   but copy into the entire region (including boundaries) of this grid.
//
// The input argument EdgeOffset is the amount the corner of this grid is
//   considered to have moved for grid comparison and copying purposes.
//   See Grid_CheckForOverlappingZones for more details.

#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void Pout( char * string, int i1 = -12345, int i2 = -12345,
           int i3 = -12345, int i4 = -12345, int i5 = -12345 );

int grid::CopyZonesFromGrid(grid *OtherGrid, FLOAT EdgeOffset[MAX_DIMENSION])
{
  

  /* Return if this doesn't involve us. */

  
  if (ProcessorNumber != MyProcessorNumber &&
      OtherGrid->ProcessorNumber != MyProcessorNumber){
    Pout( "             exit because I don't have any information");
    return SUCCESS;
  }
  
  if (NumberOfBaryonFields == 0){
    Pout( "             Exit because NumberOfBaryonFields == 0.  WTF?");
    return SUCCESS;
  }
  this->DebugCheck("CopyZonesFromGrid (before)");
  
  /* declarations */
  
  int dim;
  
  /* Compute the left and right edges of this grid (including ghost zones). */
  
  FLOAT GridLeft[MAX_DIMENSION], GridRight[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    GridLeft[dim]  = CellLeftEdge[dim][0] + EdgeOffset[dim];
    GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] + 
      CellWidth[dim][GridDimension[dim]-1]    +
      EdgeOffset[dim];
  }
  
  /* Do a quick check to see if there is any overlap. */
  
  for (dim = 0; dim < GridRank; dim++)
    if (GridLeft[dim]  >= OtherGrid->GridRightEdge[dim] ||
        GridRight[dim] <= OtherGrid->GridLeftEdge[dim]   )
      return SUCCESS;
  
  /* There is some overlap, so copy overlapping region */
  
  FLOAT Left, Right;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  int StartOther[MAX_DIMENSION], Dim[MAX_DIMENSION];
  int OtherDim[MAX_DIMENSION];
  
  /* compute start and stop indicies of overlapping region for both this
     grid and the Other grid. */
  
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Start[dim]      = 0;
    End[dim]        = 0;
    StartOther[dim] = 0;
    OtherDim[dim]   = 1;
  }
  
  for (dim = 0; dim < GridRank; dim++) 
    if (GridDimension[dim] > 1) {
      
      /* Compute left and right positions in problem space.
	 note: include buffer zones of this grid but not the other grid. */
      
      //I need to ensure corners are done correctly.
      //This is only really usefull for periodic grids where one face ISNT 
      //periodic.  
      
      if( this == OtherGrid && EdgeOffset[dim] == 0 && 0==1 ){
	Left  = max(GridLeft[dim], OtherGrid->CellLeftEdge[dim][0]);
	Right = min(GridRight[dim], OtherGrid->CellLeftEdge[dim][GridDimension[dim]-1]+
		    CellWidth[dim][GridDimension[dim]-1]);
      }else{
	Left  = max(GridLeft[dim], OtherGrid->GridLeftEdge[dim]);
	Right = min(GridRight[dim], OtherGrid->GridRightEdge[dim]);
      }
      
      
      
      /* Convert this to index positions in this grid */
      
      Start[dim] = nint((Left  - GridLeft[dim]) / CellWidth[dim][0]);
      End[dim]   = nint((Right - GridLeft[dim]) / CellWidth[dim][0]) - 1;
      
      if (End[dim] - Start[dim] < 0)
	return SUCCESS;
      
      /* Compute index positions in the other grid */
      
      StartOther[dim] = nint((Left - OtherGrid->CellLeftEdge[dim][0])/
			     CellWidth[dim][0]);
      
      /* Copy dimensions into temporary space */
      
      OtherDim[dim] = OtherGrid->GridDimension[dim];
    }
  
  /* Calculate dimensions */
  
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Dim[dim] = End[dim] - Start[dim] + 1;
  
  /* Copy data from other processor if needed (modify OtherDim and
     StartOther to reflect the fact that we are only coping part of
     the grid. */
  if(MHD_Verbose == TRUE ) Pout("                 SendRegion");
  
  // && OffProcessorHasRegion != TRUE
  if (ProcessorNumber != OtherGrid->ProcessorNumber) {
    OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber, 
				       ALL_FIELDS, NEW_ONLY, StartOther, Dim);
    for (dim = 0; dim < GridRank; dim++) {
      OtherDim[dim] = Dim[dim];
      StartOther[dim] = 0;
    }
  }
  
  
  
  /* Return if this is not our concern. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  
  /* Copy zones */
  if(MHD_Verbose == TRUE ) Pout("                  CopyZones (baryon)");
  if(MHD_Verbose == TRUE ) Pout("                     Start:",Start[0],Start[1],Start[2]);
  if(MHD_Verbose == TRUE ) Pout("                     StartOther:",StartOther[0],StartOther[1],StartOther[2]);
  if(MHD_Verbose == TRUE ) Pout("                     Dim:", Dim[0],Dim[1],Dim[2]);
  if(MHD_Verbose == TRUE ) Pout("                     EdgeOffset:", EdgeOffset[0],EdgeOffset[1],EdgeOffset[2]);
  int thisindex, otherindex, field, i, j, k;

  //For shifted periodic boundary conditions.
  //If this is the right outer boundary, on the Y or Z faces, 
  //shift by the appropriate ammount.
  int Shift[3] = {0,0,0};
  if( MHDBlastNormal[0] != 0 ){
    for(dim=1;dim<GridRank;dim++){
      if( Start[dim] > GridEndIndex[dim] &&
	  fabs( GridRightEdge[dim] - DomainRightEdge[dim] ) < CellWidth[dim][0]/100 ){
	Shift[0] = (DomainRightEdge[dim] - DomainLeftEdge[dim])*
	  (MHDBlastNormal[dim]/(MHDBlastNormal[0]*CellWidth[0][0]) );
      }
      if( End[dim] < GridStartIndex[dim] &&
	  fabs( GridLeftEdge[dim] - DomainLeftEdge[dim] ) < CellWidth[dim][0]/100 ){
	Shift[0] = -(DomainRightEdge[dim] - DomainLeftEdge[dim])*
	  (MHDBlastNormal[dim]/(MHDBlastNormal[0]*CellWidth[0][0]) );
      }//shift check
    }//dim
  }//if shift has been done.

  for (field = 0; field < NumberOfBaryonFields; field++)
    for (k = 0; k < Dim[2]; k++)
      for (j = 0; j < Dim[1]; j++) {
	thisindex = (0 + Start[0]) + (j + Start[1])*GridDimension[0] + 
	  (k + Start[2])*GridDimension[0]*GridDimension[1];
	otherindex = (0 + StartOther[0] + Shift[0]) + (j + StartOther[1] + Shift[1])*OtherDim[0] + 
	  (k + StartOther[2] + Shift[2])*OtherDim[0]*OtherDim[1];
	for(i = 0; i < Dim[0]; i++, thisindex++, otherindex++)
	    BaryonField[field][thisindex] =
	      OtherGrid->BaryonField[field][otherindex];
      }
  
  /* Centered Magentic Field */
  if( MHD_Used == TRUE ) {   
    for (field = 0; field < 3; field++)
      for (k = 0; k < Dim[2]; k++)
	for (j = 0; j < Dim[1]; j++) {
	  thisindex = (0 + Start[0]) + (j + Start[1])*GridDimension[0] + 
	    (k + Start[2])*GridDimension[0]*GridDimension[1];
	  otherindex = (0 + StartOther[0] + Shift[0]) + (j + StartOther[1]+Shift[1])*OtherDim[0] + 
	    (k + StartOther[2]+Shift[2])*OtherDim[0]*OtherDim[1];
	  for(i = 0; i < Dim[0];i++, thisindex++, otherindex++)
	    CenteredB[field][thisindex] = 
	      OtherGrid->CenteredB[field][otherindex];
	}
    
    //
    // Face Centered Magnetic Field
    //
    
    if( MHD_SendFace == TRUE ){
      
      // set up some Magnetic Field properties
      // The DimToAdd business controlls whether or not to add the first active
      // face centered field.  (i.e., Bzf(z=0) on the z boundary.)
      // Only the right edge, face centered field (i.e. Bzf(z=1) on the z face) is coppied
      // for boundary calls.  This ensures correct periodicity. dcc.
      
      int MHDDim[3][3], MHDOtherDim[3][3], MHDShift[3][3], DimToAdd=0;

      for(field=0;field<3;field++)
	for(dim=0;dim<3;dim++){
	  if(0==1){
	    MHDOtherDim[field][dim] = OtherDim[dim] + MHDAdd[field][dim];
	    
	    //These variables denote differences from the baryon fields above.
	    
	    MHDShift[field][dim] = 0;
	    MHDDim[field][dim] = Dim[dim];
	    
	    //If this routine is the fill of the Right Boundary zones,
	    //shift the Face On field by one, so as not to copy the active zone.
	    //(i.e., the x face of Bx has one extra zone, due to centering.)
	    
	    if( End[dim] == GridDimension[dim]-1 && 
		Start[dim]==GridEndIndex[dim]+1 && 
		field == dim )
	      MHDShift[field][dim] = 1;
	    
	    //If this is NOT a boundary set, then it's a bulk set, like used in
	    //Move Grid, or Rebuild Hierarchy.  In this case, the full size of the Magnetic
	    //Field needs to be accounted for.
	    //if( Start[dim] != 0 && End[dim] != DEFAULT_GHOST_ZONES - 1 )
	    
	    if(End[dim] != DEFAULT_GHOST_ZONES - 1 &&
	       Start[dim] != GridEndIndex[dim] +1 && 
	       field == dim)
	      MHDDim[field][dim]++;
	    
	    if(MHD_Verbose == TRUE ) 
	      Pout("                     (field, dim, MHDShift, MHDDim)", 
		   field, dim, MHDShift[field][dim], MHDDim[field][dim]);
	    
	  }else{
	    MHDShift[field][dim]=0;
	    DimToAdd = (End[dim] == DEFAULT_GHOST_ZONES-1 && field == dim ) ? 0 : 1;
	    //MHDDim[field][dim] = Dim[dim]+ DimToAdd;
	    MHDDim[field][dim] = Dim[dim] + ( (field == dim) ? 1 : 0 );
	    MHDOtherDim[field][dim] = OtherDim[dim] + ((field == dim) ?1:0);
	  }
	}
      
      //
      // Output stuff
      //
      if( 0==1){
	fprintf(stderr,"***********************************\n");
	fprintf(stderr,"Dim %d %d %d\n", Dim[0], Dim[1], Dim[2]);
	fprintf(stderr,"Start %d %d %d\n", Start[0], Start[1], Start[2]);
	fprintf(stderr,"OtherStart %d %d %d\n", StartOther[0], StartOther[1], StartOther[2]);
	
	for( field=0; field<3; field++){
	  
	  fprintf(stderr, "MagneticDims[%d] %d %d %d\n", field, 
		  MagneticDims[field][0],MagneticDims[field][1],MagneticDims[field][2]);
	  fprintf(stderr, "MHDDimToCopy[%d] %d %d %d\n", field, 
		  MHDDim[field][0],MHDDim[field][1],MHDDim[field][2]);
	  fprintf(stderr, "MHDOtherDim[%d] %d %d %d\n", field, 
		  MHDOtherDim[field][0],MHDOtherDim[field][1],MHDOtherDim[field][2]);
	  
	}
	
	
	fprintf(stderr,"OtherStart %d %d %d\n", StartOther[0], StartOther[1], StartOther[2]);
	fprintf(stderr,"***********************************\n");
      }
      //
      // do the copy.
      //
      //FILE *stu = fopen("file.stu", "a");    
      


      //fprintf(stderr, "kludge: Copy Zones From Gridsometimes setting Magnetic Field to Zero.\n");
      int othersize[3]={1,1,1};
      for (field =0; field<3; field++){
	
	if( MagneticField[field] == NULL ){
	  fprintf(stderr,"Severe Error: Grid_CopyZonesFromGrid.  MagneticField[%d] == NULL..\n", field);
	  return FAIL;
	}
	
	othersize[field] = MHDOtherDim[field][0]*MHDOtherDim[field][1]*MHDOtherDim[field][2];
	for( k=0; k<MHDDim[field][2]; k++)
	  for( j=0; j<MHDDim[field][1]; j++)
	    for( i=0; i<MHDDim[field][0]; i++){
	      thisindex = ( (i + Start[0]+MHDShift[field][0])
			    +(j+ Start[1]+MHDShift[field][1])*MagneticDims[field][0]
			    +(k+ Start[2]+MHDShift[field][2])*MagneticDims[field][1]*MagneticDims[field][0] );
	      
	      otherindex= ( (i + StartOther[0]+MHDShift[field][0]+Shift[0] )
			    +(j+ StartOther[1]+MHDShift[field][1]+Shift[1] )*(MHDOtherDim[field][0])
			    +(k+ StartOther[2]+MHDShift[field][2]+Shift[2] )*(MHDOtherDim[field][0]*MHDOtherDim[field][1]));
	      
	      //MagneticField[field][thisindex] = ((dccCounter < 0 )? 0 : OtherGrid->MagneticField[field][otherindex]);
	      MagneticField[field][thisindex] = OtherGrid->MagneticField[field][otherindex];
	      if( thisindex >= MagneticSize[field] ){
		fprintf(stderr, "Severe Error: Grid_CopyZonesFromGrid overstepped array bounds.\n");
		Pout(" ---------- fuck: this index >= Magnetic Size\n");
		return FAIL;
	      }
	      if(otherindex >= othersize[field]){
		fprintf(stderr, "Severe Error: Grid_CopyZonesFromGrid overstepped array bounds (other field).\n");
		return FAIL;
	      }
	    }//i
      }//field

    }//MHD_SendFace
    /* Clean up if we have transfered data. */
  }//MHD_Used       


    if (MyProcessorNumber != OtherGrid->ProcessorNumber)
      OtherGrid->DeleteAllFields();
    
  this->DebugCheck("CopyZonesFromGrid (after)");


  
  return SUCCESS;
  
}

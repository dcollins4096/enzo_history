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
#include "pout.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

extern "C" void FORTRAN_NAME(mhd_interpolate)
                               (float *parentx, float *parenty, float *parentz,
				int parentdim[], int refine[],
				float * childx, float * childy, float *childz,
				int childdim[], int childstart[], int refinedim[],
				float * otherbx, float * otherby, float * otherbz,
				int otherdim[], int otherstart[],
				float* DyBx, float* DzBx, float* DyzBx,
				float* DxBy, float* DzBy, float* DxzBy,
				float* DxBz, float* DyBz, float* DxyBz,
				int* DBxFlag,int* DByFlag,int* DBzFlag,
				FLOAT *cellsizex, FLOAT *cellsizey, FLOAT *cellsizez,
                                int *method, int * step, int * counter);


//In order to reuse CheckForOverlap, I need to get around it's fixed 
//signature by putting some extra parameters in this structure.
struct CID_Parameters{
  int Offset[MAX_DIMENSION];
  int TempDim[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION];
};
CID_Parameters CID_Params;

//MHD_ Calculate Interpolation Derivatives
int grid::MHD_CIDWorker(grid* OldFineGrid, FLOAT EdgeOffset[MAX_DIMENSION]){

  //This is a flag for the interpolator.  Step = 1 tells mhd_interpolate to only take derivatives
  //don't actually interpolate.
  int Step = 1, ver = FALSE, ver2=FALSE,Overlap, dim, i, OnlyOneFace, field;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  int StartOther[MAX_DIMENSION], Dim[MAX_DIMENSION];
  int OtherDim[MAX_DIMENSION];
  
  FLOAT GridLeft[MAX_DIMENSION], GridRight[MAX_DIMENSION], Some;
  FLOAT OverlapLeft, OverlapRight;

  char Label[30] = "MHD_CID";

  //In order to reuse CheckForOverlap, I need to get around it's fixed 
  //signature by putting some extra parameters in this structure.
  int Offset[MAX_DIMENSION];
  int TempDim[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION];

  for(dim=0;dim<GridRank;dim++){
    Offset[dim] =CID_Params.Offset[dim];
    TempDim[dim]=CID_Params.TempDim[dim];
    Refinement[dim]=CID_Params.Refinement[dim];
  }

  Pout("Pig: CID : Good processor combo.  Checking");
  OnlyOneFace = -1;
  Overlap = TRUE;
  
  //Overlap will be changed during this loop, so it's important to check at the beginning.
  for(dim=0;dim<3;dim++){
    if(Overlap == TRUE){
      
      
      //Calculate left and right edge of the entire New Fine Grid Expanded.
      //The expansion generalizes the interpolation to any refinement level, 
      //eliminating the need to re-allign ghost zones that might not line up.
      //
      //EdgeOffset is the ammount this grid has been 'periodically shifted.'
      GridLeft[dim]=CellLeftEdge[dim][0]-CellWidth[dim][0]*Offset[dim] + EdgeOffset[dim];
      GridRight[dim]=CellLeftEdge[dim][GridDimension[dim]-1]+
	CellWidth[dim][GridDimension[dim]-1]*(1+Offset[dim]) + EdgeOffset[dim];
      
      
      //"Some" is to make the ensuing floating point logic more rigorous.
      //Just comparing two floats for equality is a dumb idea, but there's no native
      //global integer position for the grid boundaries.
      
      Some = 0.1*CellWidth[dim][0];
      
      if( GridLeft[dim] > OldFineGrid->GridRightEdge[dim]+Some ) {
	Overlap = FALSE;
	if(ver==TRUE) fprintf(stderr,"MHD_CID: left too right (%d %f %f %f)\n",
			      dim, GridLeft[dim],OldFineGrid->GridRightEdge[dim], Some );
      }
      
      if( GridRight[dim] < OldFineGrid->GridLeftEdge[dim]-Some ) {
	Overlap = FALSE;
	if(ver==TRUE) fprintf(stderr,"MHD_CID: right too left (%d %f %f %f)\n",
			      dim, GridRight[dim],OldFineGrid->GridLeftEdge[dim], Some );
	
      }
      
      if( fabs(GridLeft[dim]-OldFineGrid->GridRightEdge[dim]) < Some ){
	if( OnlyOneFace == -1 )
	  OnlyOneFace = dim;
	else{
	  if(ver==TRUE) fprintf(stderr, "PFG: too many faces: %d, %d\n",OnlyOneFace, dim);
	  Overlap = FALSE;
	}
      }
      
      if( fabs(GridRight[dim]-OldFineGrid->GridLeftEdge[dim]) < Some ){
	if( OnlyOneFace == -1 )
	  OnlyOneFace = 10+dim;
	else{
	  if(ver==TRUE) fprintf(stderr, "PFG: too many faces: %d, %d\n",OnlyOneFace, dim);
	  Overlap = FALSE;
	}
      }
      
      Pout("Pig: CID: There is an overlap, dim ", dim);
      // Clearly the code has gotten this far, so some prolonging must happen.
      
      //
      // Compute the start and stop indicies of the overlapping region.
      //
      
      //initialize
      Start[dim]      = 0;
      End[dim]        = 0;
      StartOther[dim] = 0;
      OtherDim[dim]   = 1;
      
      //Determine position in Spatial coordinates
      OverlapLeft  = max(GridLeft[dim], OldFineGrid->GridLeftEdge[dim]);
      OverlapRight = min(GridRight[dim], OldFineGrid->GridRightEdge[dim]);
      //convert to GridCoordinates
      Start[dim] = nint((OverlapLeft  - GridLeft[dim]) / CellWidth[dim][0]);
      End[dim]   = nint((OverlapRight - GridLeft[dim]) / CellWidth[dim][0]) - 1;
      
      if (End[dim] - Start[dim] < 0){
	if(ver==TRUE) fprintf(stderr, "PFG: End < Start\n");
	
	Overlap = FALSE;

      }
      
      if(ver==TRUE) 
	fprintf(stderr, "MHD_CID:     dim %d OverlapLeft %f OverlapRight %f Ss %d %d\n", 
		dim, 16*OverlapLeft, 16*OverlapRight, Start[dim], End[dim]);
      
      StartOther[dim] = nint((OverlapLeft - OldFineGrid->CellLeftEdge[dim][0])/
			     CellWidth[dim][0]);
      
      //This is more than just simplification: if the two grids aren't on the same processor,
      //only the relevant information is coppied.  OtherDim must reflect that.
      OtherDim[dim] = OldFineGrid->GridDimension[dim];
      
      
    }//Overlap == TRUE
  }// dim<3
  
  //If there's no overlap, get out of this damned loop.
  if( Overlap == FALSE ){
    if(ver==TRUE) fprintf(stderr,"MHD_CID: continuing.  No overlap.\n");
    Pout("PIG: CID: No Overlap.  Next Grid");
    return SUCCESS;
    //OldFineLevelIterator = OldFineLevelIterator->NextGridThisLevel;
    //continue;
    
  }

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Dim[dim] = End[dim] - Start[dim] + 1;//dcf
  
  //If the controlling processor doesn't own the Old Subgrid, allocate space for it.
  if( MyProcessorNumber != OldFineGrid->ProcessorNumber) {
    Pout("Pig: CID: Allocating field");
    for(field=0;field<3;field++){
      if( OldFineGrid->MagneticField[field] != NULL )
	fprintf(stderr, "WTF? MHD_CID, there's a magnetic field where there shouldn't be.\n");
      
      int Size=1;
      for(dim=0;dim<3;dim++)
	Size*=Dim[dim] + ((field==dim)? 1:0);
      
      OldFineGrid->MagneticField[field] = new float[Size];
      for(i=0;i<Size;i++)
	OldFineGrid->MagneticField[field][i]=0.0;
    }
  }
  //Communicate, if I need to.
  if (ProcessorNumber != OldFineGrid->ProcessorNumber) {
    Pout("Pig: CID: Communication. (MyProcessor, Processor, OFG Processor)",
	 MyProcessorNumber, ProcessorNumber, OldFineGrid->ProcessorNumber);
    OldFineGrid->CommunicationSendRegion(OldFineGrid, ProcessorNumber, 
					 ALL_FIELDS, NEW_ONLY, StartOther, Dim);
    for (dim = 0; dim < GridRank; dim++) {
      OtherDim[dim] = Dim[dim];
      StartOther[dim] = 0;
    }
  }
  
  //The processor that doesn't have ThisGrid has nothing else to do.
  //We've gotten what we want from it, now we'll leave it to die alone in the desert.
  //Bwa ha ha ha!
  
  if (ProcessorNumber != MyProcessorNumber){
    return SUCCESS;
    //OldFineLevelIterator = OldFineLevelIterator->NextGridThisLevel;
    //continue;
  }
  
  Pout("Pig: CID: Doing Interpolation.  Should NOT be here if not in IFV call.");
  
  if(MHD_Verbose == TRUE ) Pout("pfg                  CopyZones (baryon)");
  if(MHD_Verbose == TRUE ) Pout("pfg                     Start:",Start[0],Start[1],Start[2]);
  if(MHD_Verbose == TRUE ) Pout("pfg                     StartOther:"
				,StartOther[0],StartOther[1],StartOther[2]);
  if(MHD_Verbose == TRUE ) Pout("pfg                     Dim:", Dim[0],Dim[1],Dim[2]);
  
  
  // yo
  
  //doom is a dummy field used for debugging. Feel free to remove it.
  //fprintf(stderr,"kludge:  no prolongation\n");
  int doom[6] = {1,1,1,1,1,1};
  //int doom[6] = {1,0,0,0,0,0};
  
  //Oh yeah, in my notation "Left" and "Right" refer to the "other" grid    
  int ProlongStart[3], ProlongStartOther[3], ProlongDim[3] = {1,1,1}, Method = -12;    
  int Prolong = FALSE;
  
  for( int Looper=1;Looper<=6;Looper++){
    Prolong = FALSE;
    switch(Looper){
    case 1:
      
      //This is the important flag.
      
      if( OnlyOneFace == -1 || OnlyOneFace == 0  )
	if( Start[0] > 0 && doom[0] == 1){
	  
	  for(i=0;i<3;i++){
	    ProlongStart[i] = Start[i];
	    ProlongDim[i] = End[i] - Start[i]+1;
	    ProlongStartOther[i] = StartOther[i];
	  }
	  
	  ProlongDim[0] = 2;
	  ProlongStart[0] -= 2;
	  Method = -1;
	  
	  Prolong = TRUE;
	}//x left
      
      break;
      
    case 2:
      if( OnlyOneFace == -1 || OnlyOneFace == 10 )
	if( End[0] < GridDimension[0] - 1 && doom[1] == 1){
	  
	  if( ver2==TRUE)fprintf(stderr, "Prolong x r %s\n", Label);
	  
	  for(i=0;i<3;i++){
	    ProlongStart[i] = Start[i];
	    ProlongStartOther[i] = StartOther[i];
	    ProlongDim[i] = End[i] - Start[i]+1;
	  }
	  
	  ProlongDim[0] = 2;
	  ProlongStart[0] = End[0]+1;
	  ProlongStartOther[0] = StartOther[0]+Dim[0];
	  Method = -2;
	  Prolong = TRUE;
	}//x Right
      
      break;
      
    case 3:
      if( OnlyOneFace == -1 || OnlyOneFace == 1 )
	if( Start[1] > 0 && doom[2] == 1){
	  if( ver2==TRUE)fprintf(stderr, "Prolong y l %s\n", Label);
	  
	  for(i=0;i<3;i++){
	    ProlongStart[i] = Start[i];
	    ProlongDim[i] = End[i] - Start[i]+1;
	    ProlongStartOther[i] = StartOther[i];
	  }
	  
	  ProlongDim[1] = 2;
	  ProlongStart[1] -= 2;
	  Method = -3;
	  
	  Prolong = TRUE;
	}//Y left
      
      break;
      
    case 4:
      if(OnlyOneFace == -1 || OnlyOneFace == 11 )
	if( End[1] < GridDimension[1]-1 && doom[3] == 1 ){
	  if( ver2==TRUE)fprintf(stderr, "Prolong y r %s\n", Label);
	  
	  for(i=0;i<3;i++){
	    ProlongStart[i] = Start[i];
	    ProlongStartOther[i] = StartOther[i];
	    ProlongDim[i] = End[i] - Start[i]+1;
	  }
	  
	  ProlongDim[1] = 2;
	  ProlongStart[1] = End[1]+1;
	  ProlongStartOther[1] = StartOther[1]+Dim[1];
	  Method = -4;
	  Prolong = TRUE;
	}//Y right
      
      break;
      
    case 5:
      if(OnlyOneFace == -1 || OnlyOneFace == 2)
	if( Start[2] > 0 && doom[4] == 1 ){
	  if( ver2==TRUE)fprintf(stderr, "Prolong z l %s\n", Label);
	  
	  for(i=0;i<3;i++){
	    ProlongStart[i] = Start[i];
	    ProlongDim[i] = End[i] - Start[i]+1;
	    ProlongStartOther[i] = StartOther[i];
	  }
	  
	  ProlongDim[2] = 2;
	  ProlongStart[2] -= 2;
	  Method = -5;
	  
	  Prolong = TRUE;
	}//Z left
      
      break;
      
    case 6:
      if(OnlyOneFace == -1 || OnlyOneFace == 12 )
	if( End[2] < GridDimension[2] -1 && doom[5] == 1 ){
	  if( ver2==TRUE)fprintf(stderr, "Prolong z r %s\n", Label);
	  
	  for(i=0;i<3;i++){
	    ProlongStart[i] = Start[i];
	    ProlongStartOther[i] = StartOther[i];
	    ProlongDim[i] = End[i] - Start[i]+1;
	  }
	  
	  
	  //Yes, the value for prolongstartother is correct.  It's a magneticfield.
	  ProlongDim[2] = 2;
	  ProlongStart[2] = End[2]+1;
	  ProlongStartOther[2] = StartOther[2]+Dim[2];
	  Method = -6;
	  Prolong = TRUE;
	}//Z right
      
      break;
      
    default:
      break;
      
    }//switch
    
    if(Prolong == TRUE){
      
      Pout("                     ProlongStart ",
	   ProlongStart[0],ProlongStart[1],ProlongStart[2]);
      Pout("                     ProlongStartOther ",
	   ProlongStartOther[0],ProlongStartOther[1],ProlongStartOther[2]);
      Pout("                     ProlongDim ",
	   ProlongDim[0],ProlongDim[1],ProlongDim[2]);
      
      //fprintf(stderr,"Heycf2 here\n");
      
      Step = 1;

      FORTRAN_NAME(mhd_interpolate)(this->MHDParentTemp[0], 
				    this->MHDParentTemp[1], 
				    this->MHDParentTemp[2], 
				    this->MHDParentTempPermanent,
				    this->MHDRefinementFactors,
				    MagneticField[0],
				    MagneticField[1],
				    MagneticField[2],
				    TempDim, 
				    ProlongStart,
				    ProlongDim,
				    OldFineGrid->MagneticField[0],
				    OldFineGrid->MagneticField[1],
				    OldFineGrid->MagneticField[2],
				    //OldFineGrid->GridDimension,
				    OtherDim,
				    ProlongStartOther,
				    DyBx, DzBx, DyzBx,
				    DxBy, DzBy, DxzBy,
				    DxBz, DyBz, DxyBz,
				    DBxFlag,DByFlag,DBzFlag,
				    &this->ParentDx,
				    &this->ParentDy,
				    &this->ParentDz,
				    &Method,
				    &Step,
				    &dccCounter12);
      
      char oot[30];
      sprintf(oot,"CID, loop %d",Looper);
      if( MHD_DCheck(TempDim, oot) == FAIL ){
	fprintf(stderr, "dcheck2\n");
	Pout("crap: dcheck failure, looper",Looper);
	return FAIL;
      }
    }	  
  }//Loop	
  
  // 10/4: the copy to the smaller grid will go here.
  //Remove this delete.
  // set "OffProcessorHasRegion = TRUE "
  if( MyProcessorNumber != OldFineGrid->ProcessorNumber) {
    for(field=0;field<NumberOfBaryonFields;field++){
      delete OldFineGrid->BaryonField[field];
      OldFineGrid->BaryonField[field] = NULL;
    }
    
    for(field=0;field<3;field++){
      delete OldFineGrid->MagneticField[field];
      OldFineGrid->MagneticField[field] = NULL;
    }
  }

  
  if(ver==TRUE) fprintf(stderr, "MHD_CID: OldFineLevel \n");
  return SUCCESS;
}
int grid::MHD_CID(LevelHierarchyEntry * OldFineLevel, int Offset[], 
		  int TempDim[], int Refinement[])
{
  
  if( MHD_Used != TRUE )
    return SUCCESS;

  grid * OldFineGrid;
  LevelHierarchyEntry *OldFineLevelIterator = OldFineLevel;

  int ver = FALSE;
  
  int dim;
  //In order to reuse CheckForOverlap, I need to get around it's fixed 
  //signature by putting some extra parameters in this structure.
  for(dim=0;dim<GridRank;dim++){
    CID_Params.Offset[dim] = Offset[dim];
    CID_Params.TempDim[dim]= TempDim[dim];
    CID_Params.Refinement[dim] = Refinement[dim];
  }
  //if(ver==TRUE) fprintf(stderr, "MHD_CID: Enter\n");
  if(ver==TRUE) fprintf(stderr, "MHD_CID: Offset %d %d %d\n", Offset[0],Offset[1],Offset[2]);
  if(ver==TRUE) fprintf(stderr, "MHD_CID: TempDim %d %d %d\n", 
			TempDim[0], TempDim[1], TempDim[2]);
  if(ver==TRUE) fprintf(stderr, "MHD_CID: Refinement %d %d %d\n", 
			Refinement[0], Refinement[1], Refinement[2]);
  if(ver==TRUE) fprintf(stderr, "MHD_CID: MHDParentTempPermanent %d %d %d\n", 
	  MHDParentTempPermanent[0], MHDParentTempPermanent[1], MHDParentTempPermanent[2]);

  while( OldFineLevelIterator != NULL ){
    
    OldFineGrid = OldFineLevelIterator->GridData;

    if( OldFineGrid == NULL ){
      OldFineLevelIterator = OldFineLevelIterator->NextGridThisLevel;
      continue;
    }
      
    PoutF("Pig: CID: Other Grid LeftEdge ", 
	  16*OldFineGrid->GridLeftEdge[0],16*OldFineGrid->GridLeftEdge[1],
	  16*OldFineGrid->GridLeftEdge[2]);
     Pout("Pig: CID : (OtherProcessor )", OldFineGrid->ProcessorNumber);

    if (MyProcessorNumber != ProcessorNumber && 
	MyProcessorNumber != OldFineGrid->ProcessorNumber){
      Pout("Pig: CID:    Bail.  I have neither.");
      OldFineLevelIterator = OldFineLevelIterator->NextGridThisLevel;
      continue;
    }

    //Before Interpolate Field Values:
    //Run this routine if the procesors arent the same AND this processor has the OldGrid
    //(the negation of the second clause: !(a || b ) <=> !a && !b)
    if (CommunicationDirection == COMMUNICATION_SEND &&
	(ProcessorNumber == OldFineGrid->ProcessorNumber ||
	 MyProcessorNumber != OldFineGrid->ProcessorNumber )){
      
      Pout("Pig: CID: Bail.  Send sweep: No comm necessary, or not the OtherGrids Proc");
      OldFineLevelIterator = OldFineLevelIterator->NextGridThisLevel;
      continue;
    }
    
    //Inside Interpolate Field Values:
    //Only run on the processor with the New Subgrid (this.)
    if (CommunicationDirection == COMMUNICATION_RECEIVE &&
	MyProcessorNumber != ProcessorNumber ){
      Pout("Pig: CID: Bail: Receive sweep, I'm not the New Processor");
      OldFineLevelIterator = OldFineLevelIterator->NextGridThisLevel;
      continue;
    }  

    FLOAT PeriodicOffset[3] = {0,0,0};
    this->MHD_CIDWorker( OldFineGrid,PeriodicOffset);

    //for(int fart=0;fart<MAX_DIMENSION;fart++){
    //TopLeftBoundary[fart] = reflecting;
    //TopRightBoundary[fart] = reflecting;
    //}
    dccCounter9 = 1;
    if( this->CheckForOverlap( OldFineGrid, TopLeftBoundary, TopRightBoundary,
			       &grid::MHD_CIDWorker) == FAIL ){
      fprintf(stderr,"MHD_CID: CheckForOverlap failed.\n");
      return FAIL;
    }
    dccCounter9 = 0;
    OldFineLevelIterator = OldFineLevelIterator->NextGridThisLevel;

  }


  Pout("MHD_CID: Exit\n");

  return SUCCESS;

}

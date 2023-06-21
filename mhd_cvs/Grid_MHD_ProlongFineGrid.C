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

int grid::MHD_ProlongFineGrid(grid * OtherGrid, int Step, char * Label)
{

  int ver = FALSE;
  if( MyProcessorNumber != ProcessorNumber &&
      MyProcessorNumber != OtherGrid->ProcessorNumber )
    return SUCCESS;

  this->DebugCheck("MHD_ProlongFineGrid (before)");

  if(ver==TRUE) fprintf(stderr, "PFG: --------------- Enter\n");
  Pout("PF:    Prolong Fine Grid: enter");

  // variables (OnlyOneFace will be explained in a moment)

  int dim, i, OnlyOneFace = -1;
  
  FLOAT GridLeft[MAX_DIMENSION], GridRight[MAX_DIMENSION], Some;
  FLOAT Left, Right;

  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  int StartOther[MAX_DIMENSION], Dim[MAX_DIMENSION];
  int OtherDim[MAX_DIMENSION];

  // Edges

  for(dim=0;dim<GridRank; dim++){
    if(ver==TRUE) fprintf(stderr, "PFG: dim %d \n",dim);

    GridLeft[dim]  = CellLeftEdge[dim][0];
    GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] +
      CellWidth[dim][GridDimension[dim]-1];

    //Some is to make the ensuing floating point logic more rigorous.
    //Just comparing two floats for equality is a dumb idea, but there's no native
    //global integer position for the grid boundaries.

    Some=0.1*CellWidth[dim][0];

    if( GridLeftEdge[dim] > (OtherGrid->GridRightEdge[dim]+Some) ){
      if(ver==TRUE) fprintf(stderr, "PFG: Left > Right %d %f %f\n"
	      , dim,GridLeft[dim], (OtherGrid->GridRightEdge[dim]+Some));
      return SUCCESS;
    }

    if( GridRightEdge[dim] < (OtherGrid->GridLeftEdge[dim]-Some) ){
      if(ver==TRUE) fprintf(stderr, "PFG: Right < Left\n");
      return SUCCESS;
    }


    //If these two grids share share a face (and not nested) then only one prolongation
    //needs to be called, on that face.  If two faces are shared, then the grids share an
    //edge, and no interesting prolongatin takes place. 
    
    //This abs business is because there isn't a native integer for comparison,
    //only floating point stuff.  

    if( fabs(GridLeftEdge[dim]-OtherGrid->GridRightEdge[dim]) < Some ){
      if( OnlyOneFace == -1 )
	OnlyOneFace = dim;
      else{
	if(ver==TRUE) fprintf(stderr, "PFG: too many faces: %d, %d\n",OnlyOneFace, dim);
	return SUCCESS;
      }
    }

    if( fabs(GridRightEdge[dim]-OtherGrid->GridLeftEdge[dim]) < Some ){
      if( OnlyOneFace == -1 )
	OnlyOneFace = 10+dim;
      else{
	if(ver==TRUE) fprintf(stderr, "PFG: too many faces: %d, %d\n",OnlyOneFace, dim);
	return SUCCESS;
      }
    }
    
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
    Left  = max(GridLeft[dim], OtherGrid->GridLeftEdge[dim]);
    Right = min(GridRight[dim], OtherGrid->GridRightEdge[dim]);
    if(ver==TRUE) fprintf(stderr, "PFG:     dim %d Left %f Right %f \n", dim, Left, Right);
    //convert to GridCoordinates
    Start[dim] = nint((Left  - GridLeft[dim]) / CellWidth[dim][0]);
    End[dim]   = nint((Right - GridLeft[dim]) / CellWidth[dim][0]) - 1;

    if (End[dim] - Start[dim] < 0){
      if(ver==TRUE) fprintf(stderr, "PFG: End < Start\n");
      return SUCCESS;
    }

    StartOther[dim] = nint((Left - OtherGrid->CellLeftEdge[dim][0])/
			   CellWidth[dim][0]);

    //This is more than just simplification: if the two grids aren't on the same processor,
    //only the relevant information is coppied.  OtherDim must reflect that.
    OtherDim[dim] = OtherGrid->GridDimension[dim];


  }//dim

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Dim[dim] = End[dim] - Start[dim] + 1;


  float mult = 16;
  if(ver==TRUE) fprintf(stderr, "PFG: this->GridDimension %d %d %d\n",
	  this->GridDimension[0],
	  this->GridDimension[1],
	  this->GridDimension[2]);
  if(ver==TRUE) fprintf(stderr,"PFG: GridLeft %f %f %f \n PFG: Grid Right %f %f %f\n",
	  mult*GridLeft[0],mult*GridLeft[1],mult*GridLeft[2],
	  mult*GridRight[0],mult*GridRight[1],mult*GridRight[2]);
  if(ver==TRUE) fprintf(stderr, "PFG: OtherGrid->GridDimension %d %d %d\n",
	  OtherGrid->GridDimension[0],
	  OtherGrid->GridDimension[1],
	  OtherGrid->GridDimension[2]);
  if(ver==TRUE) fprintf(stderr,"PFG: Start %d %d %d \nPFG: End %d %d %d\nPFG: Dim %d %d %d",
	  Start[0],Start[1],Start[2],End[0],End[1],End[2],Dim[0],Dim[1],Dim[2]);
  if(ver==TRUE) fprintf(stderr,"PFG: start other %d %d %d\nPFG: Odim %d %d %d\n",
	  StartOther[0],StartOther[1],StartOther[2],
	  OtherDim[0],OtherDim[1],OtherDim[2]);



  if (ProcessorNumber != OtherGrid->ProcessorNumber) {
    OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber, 
				       ALL_FIELDS, NEW_ONLY, StartOther, Dim);
    for (dim = 0; dim < GridRank; dim++) {
      OtherDim[dim] = Dim[dim];
      StartOther[dim] = 0;
    }
  }


  //The processor that doesn't have ThisGrid has nothing else to do.
  //We've gotten what we want from it, now we'll leave it to die alone in the desert.
  //Bwa ha ha ha!

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if(MHD_Verbose == TRUE ) Pout("pfg                  CopyZones (baryon)");
  if(MHD_Verbose == TRUE ) Pout("pfg                     Start:",Start[0],Start[1],Start[2]);
  if(MHD_Verbose == TRUE ) Pout("pfg                     StartOther:"
				,StartOther[0],StartOther[1],StartOther[2]);
  if(MHD_Verbose == TRUE ) Pout("pfg                     Dim:", Dim[0],Dim[1],Dim[2]);


  //
  // here is the prolong itself.
  //
      
  if( this->MHDParentTemp[0] == NULL ){
    if(ver==TRUE) fprintf(stderr, "PFG: What's wrong with the Temp Grid?\n");
    return FAIL;
  }


  if( dccCounter > 0 ){
    //dccCounter12 += 10;
    fprintf(stderr, "PFG FLAG: setting derivatives to %d.  \n", dccCounter12);

  }
  if(this->MHDParentTemp[0] != NULL ){
	
    //doom is a dummy field used for debugging. Feel free to remove it.
    //fprintf(stderr,"kludge:  no prolongation\n");
    int doom[6] = {1,1,1,1,1,1};
    //int doom[6] = {1,0,0,0,0,0};

    int ProlongStart[3], ProlongStartOther[3], ProlongDim[3] = {1,1,1}, Method = -12;    
    int Prolong = FALSE;
    
    for( int Looper=1;Looper<=6;Looper++){
      Prolong = FALSE;
      switch(Looper){
      case 1:

	//This is the important flag.
	if( OnlyOneFace == -1 || OnlyOneFace == 0 ){
	  
	  //These logical flags are probably redundant.  They're in this code because I had them
	  //in the other code, where they were a.) necessary and b.) almost did their job right.
	  if( Start[0] > GridStartIndex[0] && 
	      Start[0] < GridEndIndex[0] &&
	      doom[0] == 1 ){
	    fprintf(stderr, "Prolong x l %s\n", Label);
	    for(i=0;i<3;i++){
	      ProlongStart[i] = Start[i];
	      ProlongDim[i] = End[i] - Start[i]+1;
	      ProlongStartOther[i] = StartOther[i];
	    }
	    
	    ProlongDim[0] = 2;
	    ProlongStart[0] -= 2;
	    Method = -1;
	    
	    Prolong = TRUE;
	  }//Start Logicals.
	}//OnlyOneFace
	
	break;
	
      case 2:
	if( OnlyOneFace == -1 || OnlyOneFace == 10 ){
	  if( End[0] < GridEndIndex[0] &&
	      End[0] > GridStartIndex[0] && 
	      doom[1] == 1 ){
	    fprintf(stderr, "Prolong x r %s\n", Label);
	    
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
	  }
	}//OnlyOneFace
	break;
	
      case 3:
	if( OnlyOneFace == -1 || OnlyOneFace == 1 ){
	  if( Start[1] > GridStartIndex[1] &&
	      Start[1] < GridEndIndex[1] && 
	      doom[2] == 1 ){
	    fprintf(stderr, "Prolong y l %s\n", Label);
	    
	    for(i=0;i<3;i++){
	      ProlongStart[i] = Start[i];
	      ProlongDim[i] = End[i] - Start[i]+1;
	      ProlongStartOther[i] = StartOther[i];
	    }
	    
	    ProlongDim[1] = 2;
	    ProlongStart[1] -= 2;
	    Method = -3;
	      
	    Prolong = TRUE;
	  }
	}//OnlyOneFace
	break;
	
      case 4:
	if(OnlyOneFace == -1 || OnlyOneFace == 11 ){
	  if( End[1] < GridEndIndex[1] &&
	      End[1] > GridStartIndex[1] && 
	      doom[3] == 1 ){
	    fprintf(stderr, "Prolong y r %s\n", Label);
	    
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
	  }
	}//OnlyOneFace
	break;
	
      case 5:
	if(OnlyOneFace == -1 || OnlyOneFace == 2){
	  if( Start[2] > GridStartIndex[2] && 
	      Start[2] < GridEndIndex[2] &&
	      doom[4] == 1 ){
	    fprintf(stderr, "Prolong z l %s\n", Label);
	    
	    for(i=0;i<3;i++){
	      ProlongStart[i] = Start[i];
	      ProlongDim[i] = End[i] - Start[i]+1;
	      ProlongStartOther[i] = StartOther[i];
	    }
	    
	    ProlongDim[2] = 2;
	    ProlongStart[2] -= 2;
	    Method = -5;
	    
	    Prolong = TRUE;
	  }
	}//OnlyOneface
	break;
	
      case 6:
	if(OnlyOneFace == -1 || OnlyOneFace == 12 ){    
	  if( End[2] < GridEndIndex[2] &&
	      End[2] > GridStartIndex[2] &&
	      doom[5] == 1 ){
	    fprintf(stderr, "Prolong z r %s\n", Label);
	    
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
	  }
	}//OnlyOneFace
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
	FORTRAN_NAME(mhd_interpolate)(this->MHDParentTemp[0], 
				      this->MHDParentTemp[1], 
				      this->MHDParentTemp[2], 
				      this->MHDParentTempPermanent,
				      this->MHDRefinementFactors,
				      MagneticField[0],
				      MagneticField[1],
				      MagneticField[2],
				      GridDimension,
				      ProlongStart,
				      ProlongDim,
				      OtherGrid->MagneticField[0],
				      OtherGrid->MagneticField[1],
				      OtherGrid->MagneticField[2],
				      OtherGrid->GridDimension,
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

      }	  
    }//Loop	
  }//if(this->MHDParentTemp[0] != NULL )


  Pout("PF:    Prolong Fine Grid: exit");

  this->DebugCheck("MHD_ProlongFineGrid (after)");

  return SUCCESS;
}

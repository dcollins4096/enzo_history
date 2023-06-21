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
/  GRID CLASS (RECEIVES FROM 'FAKE' GRID TO REAL GRID)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "error.h"
#include "pout.h"
/* function prototypes */

extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest, 
                                   int *sdim1, int *sdim2, int *sdim3, 
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3, 
                                   int *dstart1, int *dstart2, int *dststart3);

extern "C" void FORTRAN_NAME(combine3d)(
               float *source1, float *weight1, float *source2, float *weight2,
	       float *dest, int *sdim1, int *sdim2, int *sdim3, 
	       int *ddim1, int *ddim2, int *ddim3,
	       int *sstart1, int *sstart2, int *sstart3, 
	       int *dstart1, int *dstart2, int *dstart3,
	       int *ivel_flag, int *irefine);

float ReturnCPUTime();

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */


int grid::CommunicationReceiveRegion(grid *FromGrid, int FromProcessor, 
				     int SendField, int NewOrOld, 
				     int RegionStart[], int RegionDim[],
				     int IncludeBoundary)
{

  int i, index=0, field, dim, Zero[] = {0, 0, 0};
  float One = 1.0;
  /* Compute size of region to transfer. */

  //this->ProcessorNumber, FromProcessor, FromGrid->ProcessorNumber, MyProcessorNumber);

  int SendAllBaryonFields = FALSE;
  switch( SendField ){
  case ALL_FIELDS:
  case JUST_BARYONS:
  case BARYONS_ELECTRIC:
  case BARYONS_MAGNETIC:
    SendAllBaryonFields = TRUE;
    break;
  default:
    SendAllBaryonFields = FALSE;
    break;
  }
  if(   dccTempVerbosity == TRUE ){
    Pout("CRR: Just Entered.  SendAllBaryonFields",
	    SendAllBaryonFields);
  }

  //Count the number of floats to send.
  //This logic has gotten a bit more complicated with the advent of MHD.
  /* in the original version, the logic looked like this.  This is only for reference.
  int NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
                       ((NewOrOld == NEW_AND_OLD)? 2 : 1);
  */

  int NumberOfFields = 0;
  if( SendAllBaryonFields == TRUE ) 
    NumberOfFields = NumberOfBaryonFields;

  //if SendField >= 0 then send only that field.
  if( SendField >= 0 )
    NumberOfFields = 1;

  if( NewOrOld == NEW_AND_OLD )
    NumberOfFields *= 2;

  int RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];

  /* Allocate buffer. */
  
  int TransferSize = RegionSize*NumberOfFields;
  float *buffer = NULL;

  int FromDim[MAX_DIMENSION], FromOffset[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    FromOffset[dim] = (dim < GridRank && IncludeBoundary == FALSE)? 
      DEFAULT_GHOST_ZONES : 0;
    FromDim[dim] = RegionDim[dim] + 2*FromOffset[dim];
  }

  //MHD if for the Magnetic Field
  //MHDe is for the ElectricField
  
  int MHDRegionDim[3][3], MHDRegionSize[3]={1,1,1}, MHDFromDim[3][3];
  int MHDeRegionDim[3][3], MHDeRegionSize[3]={1,1,1}, MHDeFromDim[3][3];
  int ThisIsAFaceProjection=FALSE;
  int MHD_SendBFlag[3]={FALSE,FALSE,FALSE};
  int MHD_SendEFlag[3]={FALSE,FALSE,FALSE};
  //This is used for determining the proper Electric dimensions.
  int MHD_BoundaryOnly[3]={FALSE,FALSE,FALSE};

  //This complicated Field Specific Send business is only important in MHD_ProjectFace.
  if( MHD_Used ){
    switch(SendField){
    case ALL_FIELDS:
      for(dim=0;dim<3;dim++){
	MHD_SendBFlag[dim]=TRUE;
	MHD_SendEFlag[dim]=TRUE;  
	MHD_BoundaryOnly[dim]=FALSE;
      }
      break;

    case BARYONS_MAGNETIC:
      for(dim=0;dim<3;dim++){
	MHD_SendBFlag[dim]=TRUE;
        MHD_SendEFlag[dim]=FALSE;
        MHD_BoundaryOnly[dim]=FALSE;
      }
      break;
    case BARYONS_ELECTRIC:
      for(dim=0;dim<3;dim++){
	MHD_SendBFlag[dim]=FALSE;
        MHD_SendEFlag[dim]=TRUE;
        MHD_BoundaryOnly[dim]=FALSE;
      }
      break;

    case JUST_BARYONS:
      for(dim=0;dim<3;dim++){
	MHD_SendBFlag[dim]=FALSE;
	MHD_SendEFlag[dim]=FALSE;
      }
      break;


    case ELECTRIC_FIELD:
      for(dim=0;dim<3;dim++){
	MHD_SendEFlag[dim] = ( MHD_ProjectThisFace[dim] == TRUE ) ? FALSE : TRUE;
	MHD_BoundaryOnly[dim] = MHD_ProjectThisFace[dim]; 

      }//dim
      break;
      
    case MAGNETIC_FIELD:

      for(dim=0;dim<3;dim++){

	if( MHD_ProjectThisFace[dim] == TRUE ){
	  MHD_BoundaryOnly[dim] = TRUE;
	  MHD_SendBFlag[dim] = TRUE;
	}else{
	  MHD_BoundaryOnly[dim] = FALSE;
	  MHD_SendBFlag[dim]=FALSE;
	}

      }//dim

      break;
    default:
      break;
    }

    // For the cell centered magnetic field:
    TransferSize += ((SendField == ALL_FIELDS)? 3*RegionSize : 0 )*
      ((NewOrOld == NEW_AND_OLD)? 2 : 1);
    
    //
    // Face Centered Magnetic Field 
    //

    for(field =0; field<3;field++){

      if( MHD_SendBFlag[field] == TRUE ){
	
	//Calculate size of transfer region
	
	for(dim=0;dim<3;dim++){
	  
	  //Only expand MHD region (in standard way) if 
	  //NOT a face projection.
	  
	  MHDRegionDim[field][dim] = RegionDim[dim]
	    +( (MHD_BoundaryOnly[dim]==TRUE)? 0: ( (field==dim) ? 1:0) );
	
	  MHDRegionSize[field] *= MHDRegionDim[field][dim];
	  MHDFromDim[field][dim] = FromDim[dim] +( (field==dim) ? 1:0 ) ;
	}//dim
      	
	//increase allocation for Face Centered field
   
	TransferSize += MHDRegionSize[field] *
	  ((NewOrOld == NEW_AND_OLD)? 2 : 1);

      }//SendBFlag
    }//field
    
    //
    // Transfer size for electric field
    //

    for(field=0;field<3;field++){

      if(MHD_SendEFlag[field]==TRUE){
	
	for( dim=0; dim<3;dim++){
	  MHDeRegionDim[field][dim] = RegionDim[dim] + 
	    ( (field==dim) ? 0 : ( (MHD_BoundaryOnly[dim]==TRUE) ? 0:1 ) );
	  MHDeRegionSize[field] *= MHDeRegionDim[field][dim];
	  MHDeFromDim[field][dim] = FromDim[dim] +( (field==dim) ? 0 : 1 );
	}
	
	TransferSize += MHDeRegionSize[field];
	
      }//Efield
    }//field
    
  }//MHD_Used
  
 

  if (MyProcessorNumber == FromProcessor || 
      MyProcessorNumber == ProcessorNumber){

    buffer = new float[TransferSize];
  }
  
  /* If this is the from processor, pack fields. */
  if(   dccTempVerbosity == TRUE ) Pout("CRR: Corn");
  if (MyProcessorNumber == FromProcessor) {
    
    index = 0;
    
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < FromGrid->NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS || SendAllBaryonFields == TRUE) {
	  if(   dccTempVerbosity == TRUE ){
	    Pout("CRR: Damnit. Sending BF");
	    Pout("CRR: field SendField ALL_FIELDS SendAll TRUE",
		    field, SendField, ALL_FIELDS, SendAllBaryonFields, TRUE);
	  }
	  FORTRAN_NAME(copy3d)(FromGrid->BaryonField[field], &buffer[index],
			       FromDim, FromDim+1, FromDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       FromOffset, FromOffset+1, FromOffset+2);
	  index += RegionSize;
	}
    
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < FromGrid->NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS || SendAllBaryonFields == TRUE ){
	  FORTRAN_NAME(copy3d)(FromGrid->OldBaryonField[field], &buffer[index],
			       FromDim, FromDim+1, FromDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       FromOffset, FromOffset+1, FromOffset+2);
	  index += RegionSize;
	}
    
    if( MHD_Used ){
      
      /* Send Centered B */
      if( NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY )
	for(field = 0;field<3;field++)
	  if( SendField == ALL_FIELDS){
	    FORTRAN_NAME(copy3d)(FromGrid->CenteredB[field], &buffer[index],
				 FromDim, FromDim+1, FromDim+2,
				 RegionDim,RegionDim+1, RegionDim+2,
				 Zero, Zero+1, Zero+2,
				 FromOffset, FromOffset+1, FromOffset+2);
	    
	    
	    index += RegionSize;
	    
	  }//field
      
      
      if( NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY )
	for(field = 0;field<3;field++)
	  if( SendField==ALL_FIELDS && 0==1){
	    FORTRAN_NAME(copy3d)(FromGrid->OldCenteredB[field], &buffer[index],
				 FromDim, FromDim+1, FromDim+2,
				 RegionDim,RegionDim+1, RegionDim+2,
				 Zero, Zero+1, Zero+2,
				 FromOffset, FromOffset+1, FromOffset+2);
	    index += RegionSize;
	    
	  }//field
      
      /* send Face B */
      if( NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY )
	for( field=0;field<3;field++)
	  if( MHD_SendBFlag[field]==TRUE){

	    FORTRAN_NAME(copy3d)(FromGrid->MagneticField[field], &buffer[index],
				 MHDFromDim[field], MHDFromDim[field]+1, MHDFromDim[field]+2,
				 MHDRegionDim[field], MHDRegionDim[field]+1, MHDRegionDim[field]+2,
				 Zero, Zero+1, Zero+2,
				 FromOffset, FromOffset+1, FromOffset+2);

       
	    index += MHDRegionSize[field];
	  }
      
      if( NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY )
	  for( field=0;field<3;field++)
	  if( MHD_SendBFlag[field]==TRUE){
	    FORTRAN_NAME(copy3d)(FromGrid->OldMagneticField[field], &buffer[index],
				 MHDFromDim[field], MHDFromDim[field]+1, MHDFromDim[field]+2,
				 MHDRegionDim[field], MHDRegionDim[field]+1, MHDRegionDim[field]+2,
				 Zero, Zero+1, Zero+2,
				 FromOffset, FromOffset+1, FromOffset+2);
	    
	    
	    index += MHDRegionSize[field];
	  }
      

      for( field=0; field<3; field++)
	if( MHD_SendEFlag[field]==TRUE){	
	  
	  FORTRAN_NAME(copy3d)(FromGrid->ElectricField[field], &buffer[index],
		       MHDeFromDim[field], MHDeFromDim[field]+1, MHDeFromDim[field]+2,
		       MHDeRegionDim[field], MHDeRegionDim[field]+1, MHDeRegionDim[field]+2, 
		       Zero, Zero+1, Zero+2,
		       FromOffset, FromOffset+1, FromOffset+2);
	  index += MHDeRegionSize[field];
	}//field
      
      
    }//MHD_Used
  }
  /* Send buffer. */
         
  
#ifdef USE_MPI

  /* only send if processor numbers are not identical */

  if (ProcessorNumber != FromProcessor) {
    
    MPI_Status status;
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;

//    MPI_Datatype DataType = MPI_FLOAT;
//    if (sizeof(float) == 8)
//      DataType = MPI_DOUBLE;

    float time1 = ReturnCPUTime();
    ZLAN_START;

    if (MyProcessorNumber == FromProcessor) {
//      fprintf(stderr, "RF: Sending %d floats from %d to %d\n", TransferSize, 
//	      FromProcessor, ProcessorNumber);

      CommunicationBufferedSend(buffer, TransferSize, DataType, 
		       ProcessorNumber, 0, MPI_COMM_WORLD, BUFFER_IN_PLACE);



    }

    if (MyProcessorNumber == ProcessorNumber) {


      JBPERF_START_MPI_RECV("MPI_Recv",TransferSize,DataType);

      CHECK_MPI_ERROR(MPI_Recv(buffer, TransferSize, DataType, FromProcessor, 
			       0, MPI_COMM_WORLD, &status));

      JBPERF_STOP_MPI_RECV("MPI_Recv",TransferSize,DataType);




    }

//    if (MyProcessorNumber == FromProcessor) {
//      fprintf(stderr, "RF: Sending %d floats from %d to %d\n", TransferSize,
//            FromProcessor, ProcessorNumber);
//      CHECK_MPI_ERROR(MPI_Bsend(buffer, TransferSize, DataType, 
//                      ProcessorNumber, 0, MPI_COMM_WORLD));
//    }

//    if (MyProcessorNumber == ProcessorNumber) {
//      fprintf(stderr, "RF: Waiting for %d floats at %d from %d\n",
//            TransferSize, MyProcessorNumber, FromProcessor);
//      CHECK_MPI_ERROR(MPI_Recv(buffer, TransferSize, DataType, FromProcessor,
//                      0,MPI_COMM_WORLD, &status));
//    }

    ZLAN_STOP_RECV(5);

    ZLAN_COUNT(6,TransferSize);

    CommunicationTime += ReturnCPUTime() - time1;

  }


#endif /* USE_MPI */

  /* If this is the to processor, unpack fields. */

  int GridSize = GridDimension[0]*GridDimension[1]*GridDimension[2];

  if (MyProcessorNumber == ProcessorNumber) {


    index = 0;
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS || SendAllBaryonFields == TRUE ){
	  if (BaryonField[field] == NULL) {
	    BaryonField[field] = new float[GridSize];
	    for (i = 0; i < GridSize; i++)
	      BaryonField[field][i] = 0;
          }
	  FORTRAN_NAME(copy3d)(&buffer[index], BaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionStart, RegionStart+1, RegionStart+2,
			       Zero, Zero+1, Zero+2);

	  index += RegionSize;
	}
    
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS || SendAllBaryonFields == TRUE) {
	  if (OldBaryonField[field] == NULL) {
	    OldBaryonField[field] = new float[GridSize];
	    for (i = 0; i < GridSize; i++)
	      BaryonField[field][i] = 0;
          }
	  FORTRAN_NAME(copy3d)(&buffer[index], OldBaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionStart, RegionStart+1, RegionStart+2,
			       Zero, Zero+1, Zero+2);

	  index += RegionSize;
	}

    if( MHD_Used ){     
      /* unpack centeredB */
      if( NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY )
	if( SendField == ALL_FIELDS )
	  for(field = 0; field<3; field++){
	    
	    if(CenteredB[field] == NULL){
	      CenteredB[field] = new float[GridSize];
	      for(i=0;i<GridSize;i++) CenteredB[field][i] = 0.0;
	    }//allocate Bc
	    
	    
	    FORTRAN_NAME(copy3d)(&buffer[index], CenteredB[field],
				 RegionDim, RegionDim+1, RegionDim+2,
				 GridDimension, GridDimension+1, GridDimension+2,
				 RegionStart, RegionStart+1, RegionStart+2,
				 Zero, Zero+1, Zero+2);
	    
	    index += RegionSize;
	    
	  }//new Bc
      
      // I don't think I'm ever going to want the OldCenteredB
      // This is here "just in case"
      /*
	if( NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY )
	if(SendField == ALL_FIELDS && 0==1)
	for(field = 0; field<3; field++){
	
	if(OldCenteredB[field] == NULL){
	OldCenteredB[field] = new float[GridSize];
	for(i=0;i<GridSize;i++) OldCenteredB[field][i] = 0.0;
	}//allocate Bc
	
	FORTRAN_NAME(copy3d)(&buffer[index], OldCenteredB[field],
	RegionDim, RegionDim+1, RegionDim+2,
	GridDimension, GridDimension+1, GridDimension+2,
	RegionStart, RegionStart+1, RegionStart+2,
	Zero, Zero+1, Zero+2);
	
	index += RegionSize;
	
	}//Old Bc
      */

      /* unpack face B */
      if( NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY )
	for(field = 0; field<3; field++)
	  if( MHD_SendBFlag[field]==TRUE){

	    if(MagneticField[field] == NULL){
	      MagneticField[field] = new float[MagneticSize[field] ];
	      for(i=0;i<MagneticSize[field];i++) MagneticField[field][i] = 0.0;
	    }//allocate Bf
	    
	    
	    FORTRAN_NAME(copy3d)(&buffer[index], MagneticField[field],
				 MHDRegionDim[field], MHDRegionDim[field]+1, MHDRegionDim[field]+2,
				 MagneticDims[field], MagneticDims[field]+1, MagneticDims[field]+2,
				 RegionStart, RegionStart+1, RegionStart+2,
				 Zero, Zero+1, Zero+2);
	    

	    index += MHDRegionSize[field];
	    
	  }//unpack new bf
    
      if( NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY )

	for(field = 0; field<3; field++)
	  if( MHD_SendBFlag[field]==TRUE ){
	    
	    if(OldMagneticField[field] == NULL){
	      OldMagneticField[field] = new float[MagneticSize[field] ];
	      for(i=0;i<MagneticSize[field];i++) OldMagneticField[field][i] = 0.0;
	    }//allocate Bf
	    
	    FORTRAN_NAME(copy3d)(&buffer[index], OldMagneticField[field],
				 MHDRegionDim[field], MHDRegionDim[field]+1, MHDRegionDim[field]+2,
				 MagneticDims[field], MagneticDims[field]+1, MagneticDims[field]+2,
				 RegionStart, RegionStart+1, RegionStart+2,
				 Zero, Zero+1, Zero+2);
	    
	    index += MHDRegionSize[field];
	    
	  }//unpack old bf
      
      for(field=0; field<3;field++)
	if(MHD_SendEFlag[field]==TRUE){

	  if( ElectricField[field] == NULL){
	    ElectricField[field] = new float[ ElectricSize[field] ];
	    for(i=0;i<ElectricSize[field];i++)
	      ElectricField[field][i] = 0.0;
	  }

	  FORTRAN_NAME(copy3d)(&buffer[index], ElectricField[field],
		       MHDeRegionDim[field], MHDeRegionDim[field]+1, MHDeRegionDim[field]+2,
		       ElectricDims[field],ElectricDims[field]+1,ElectricDims[field]+2,
		       RegionStart, RegionStart + 1, RegionStart + 2,
		       Zero, Zero+1, Zero+2);
	  index += MHDeRegionSize[field];
	  

	}
            
    }//MHD_Used

  }//unpack

  /* Clean up. */

  if (MyProcessorNumber == ProcessorNumber)
    delete [] buffer;


  return SUCCESS;
}


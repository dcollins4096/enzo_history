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
/  GRID CLASS (SEND FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
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
#include "performance.h"
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "error.h"

/* function prototypes */


extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest, 
                                   int *sdim1, int *sdim2, int *sdim3, 
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3, 
                                   int *dstart1, int *dstart2, int *dststart3);
float ReturnCPUTime();

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */


int grid::CommunicationSendRegion(grid *ToGrid, int ToProcessor,int SendField, 
			      int NewOrOld, int RegionStart[], int RegionDim[])
{

  int index, field, dim, Zero[] = {0, 0, 0};

  /* Compute size of region to transfer. */

  int NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
                       ((NewOrOld == NEW_AND_OLD)? 2 : 1);
  if (SendField == ACCELERATION_FIELDS)
    NumberOfFields = GridRank;
  int RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];
  int TransferSize = RegionSize * NumberOfFields;

  /* MHD Dimension stuff */

  int MHDRegionDim[3][3], MHDRegionSize[3]={1,1,1};
  int MHDeRegionDim[3][3], MHDeRegionSize[3]={1,1,1};

  if( MHD_Used ){
    //Account for face centered field.  Note that I don't want to communicate the OldCenteredField
    TransferSize += ((SendField == ALL_FIELDS)? 3*RegionSize : 0 )*
                     ((NewOrOld == NEW_AND_OLD)? 2 : 1);
   
    if(MHD_SendFace == TRUE){

      for(field =0; field<3;field++){
	for(dim=0;dim<3;dim++){
	  MHDRegionDim[field][dim] = RegionDim[dim]+MHDAdd[field][dim];
	  MHDRegionSize[field] *= MHDRegionDim[field][dim];
	  MHDeRegionDim[field][dim]= RegionDim[dim]+( (field==dim)?0:1);
	  MHDeRegionSize[field] *=MHDeRegionDim[field][dim];
	}
	
	TransferSize += ((SendField == ALL_FIELDS)? MHDRegionSize[field]: 0 )*
	  ((NewOrOld == NEW_AND_OLD)? 2 : 1);
	
      }//field
    }//MHD_SendFace
    
  }//MHD_Used

  
  /* Allocate buffer. */
  
  float *buffer = NULL;
  if (MyProcessorNumber == ProcessorNumber || MyProcessorNumber == ToProcessor)
    buffer = new float[TransferSize];
  
  /* If this is the from processor, pack fields. */
  if (MyProcessorNumber == ProcessorNumber) {
    
    /*    printf("SendRegion: RegionStart = %d %d %d\n", RegionStart[0], 
	  RegionStart[1], RegionStart[2]); */
    
    index = 0;
    
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY){
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(BaryonField[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
	}
    }

    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY){
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(OldBaryonField[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
	}
      

    }

    if( MHD_Used && SendField == ALL_FIELDS ){
      
      if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY ){

	for( field = 0; field<3; field++){
	  if(CenteredB[field] == NULL ){
	    fprintf(stderr, "Severe Error in Grid_CommunicationSendRegion:  CenteredB = NULL.");
	  }
	  FORTRAN_NAME(copy3d)(CenteredB[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
	}

      }

      if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY){
	for( field = 0; field<3; field++){
	  FORTRAN_NAME(copy3d)(OldCenteredB[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
	}
            

      }
      
      if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY )
	for(field=0;field<3;field++){
	  FORTRAN_NAME(copy3d)(MagneticField[field], &buffer[index],
			       &MagneticDims[field][0], 
			       &MagneticDims[field][1], 
			       &MagneticDims[field][2], 
			       &MHDRegionDim[field][0],
			       &MHDRegionDim[field][1],
			       &MHDRegionDim[field][2],
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += MHDRegionSize[field];
	}

      
      if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY){
	for(field=0;field<3;field++){
	  FORTRAN_NAME(copy3d)(OldMagneticField[field], &buffer[index],
			       &MagneticDims[field][0],
			       &MagneticDims[field][1],
			       &MagneticDims[field][2],
			       &MHDRegionDim[field][0],
			       &MHDRegionDim[field][1],
			       &MHDRegionDim[field][2],
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += MHDRegionSize[field];
	  
	  
	}
      }//new and old or old only
    }//MHD
  

    if (SendField == GRAVITATING_MASS_FIELD_PARTICLES)
      FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles, buffer,
			   GravitatingMassFieldParticlesDimension,
			   GravitatingMassFieldParticlesDimension+1,
			   GravitatingMassFieldParticlesDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
    
    if (SendField == GRAVITATING_MASS_FIELD)
      FORTRAN_NAME(copy3d)(GravitatingMassField, buffer,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
    
    if (SendField == POTENTIAL_FIELD)
      FORTRAN_NAME(copy3d)(PotentialField, buffer,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
    
    if (SendField == ACCELERATION_FIELDS)
      for (dim = 0; dim < GridRank; dim++) {
	FORTRAN_NAME(copy3d)(AccelerationField[dim], &buffer[index],
			     GridDimension, GridDimension+1, GridDimension+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     RegionStart, RegionStart+1, RegionStart+2);
	index += RegionSize;
      }
  }
  
  /* Send buffer. */
  
#ifdef USE_MPI
  
  /* only send if processor numbers are not identical */
  
  if (ProcessorNumber != ToProcessor) {
    
    MPI_Status status;
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    
    //    MPI_Datatype DataType = MPI_FLOAT;
    //    if (sizeof(float) == 8)
    //      DataType = MPI_DOUBLE;
    
    float time1 = ReturnCPUTime();
    ZLAN_START;
    
    if (MyProcessorNumber == ProcessorNumber) {

      //      fprintf(stderr, "Sending %d floats from %d to %d\n", TransferSize, 
      //	      MyProcessorNumber, ToProcessor);
      CommunicationBufferedSend(buffer, TransferSize, DataType, ToProcessor, 
				MPI_SENDREGION_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);

    }
    
    
    if (MyProcessorNumber == ToProcessor) {
      //      cerr << "\a";

      //      fprintf(stderr, "Waiting for %d floats at %d from %d\n", TransferSize, 
      //	      MyProcessorNumber, ProcessorNumber);
      JBPERF_START_MPI_RECV("MPI_Recv",TransferSize, DataType);
      CHECK_MPI_ERROR(MPI_Recv(buffer, TransferSize, DataType, ProcessorNumber,
			       MPI_SENDREGION_TAG, MPI_COMM_WORLD, &status));
      JBPERF_STOP_MPI_RECV("MPI_Recv",TransferSize, DataType);

      
    }

//    if (MyProcessorNumber == ProcessorNumber) {
//      fprintf(stderr, "Sending %d floats from %d to %d\n", TransferSize,
//            MyProcessorNumber, ToProcessor);
//      CHECK_MPI_ERROR(MPI_Bsend(buffer, TransferSize, DataType, ToProcessor,
//               MPI_SENDREGION_TAG, MPI_COMM_WORLD));
//    }

//    if (MyProcessorNumber == ToProcessor) {
//      fprintf(stderr, "Waiting for %d floats at %d from %d\n", TransferSize,
//            MyProcessorNumber, ProcessorNumber);
//      CHECK_MPI_ERROR(MPI_Recv(buffer, TransferSize, DataType, ProcessorNumber,
//               MPI_SENDREGION_TAG, MPI_COMM_WORLD, &status));
//    }

    ZLAN_STOP_RECV(5);

    ZLAN_COUNT(6,TransferSize);

    CommunicationTime += ReturnCPUTime() - time1;

  }
  

#endif /* USE_MPI */

  /* If this is the to processor, unpack fields. */

  if (MyProcessorNumber == ToProcessor) {

    
    index = 0;
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY){
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  delete ToGrid->BaryonField[field];
	  ToGrid->BaryonField[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->BaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	}


    }
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY){
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  delete ToGrid->OldBaryonField[field];
	  ToGrid->OldBaryonField[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->OldBaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	}

    }
    if( MHD_Used && SendField == ALL_FIELDS ){
      if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY){
	for(field=0;field<3;field++){
	  delete ToGrid->CenteredB[field];
	  ToGrid->CenteredB[field]=new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], CenteredB[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	  
	}

      }


      if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY){

	
	for(field=0;field<3;field++){
	  delete ToGrid->OldCenteredB[field];
	  ToGrid->OldCenteredB[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], OldCenteredB[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	  
	}

      }



      /* send Bf */
      if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY){
	if( MHD_SendFace == TRUE ){
	  for(field=0;field<3;field++){
	    delete ToGrid->MagneticField[field];
	    ToGrid->MagneticField[field] = new float[MHDRegionSize[field] ];
	    FORTRAN_NAME(copy3d)(&buffer[index], MagneticField[field],
				 &MHDRegionDim[field][0],
				 &MHDRegionDim[field][1],
				 &MHDRegionDim[field][2],
				 &MHDRegionDim[field][0],
				 &MHDRegionDim[field][1],
				 &MHDRegionDim[field][2],
				 Zero, Zero+1, Zero+2,
				 Zero, Zero+1, Zero+2);
	    index += MHDRegionSize[field];

	  }

	}
      }

      if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
	for(field=0;field<3;field++){
	  
	  if( OldMagneticField[field] != NULL )
	    fprintf(stderr,"shit: CommSendRegion: OldMagneticField != NULL" );

	  delete ToGrid->OldMagneticField[field];
	  ToGrid->OldMagneticField[field] = new float[MHDRegionSize[field]];
	  
	  FORTRAN_NAME(copy3d)(&buffer[index], OldMagneticField[field],
			       &MHDRegionDim[field][0],
			       &MHDRegionDim[field][1],
			       &MHDRegionDim[field][2],
			       &MHDRegionDim[field][0],
			       &MHDRegionDim[field][1],
			       &MHDRegionDim[field][2],
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += MHDRegionSize[field];
	}

      
      //Allocate E: This is really only usefull for the serial Root Grid IO case,
      //where PartitionGrid is used.  
      if(0==1)
      for(field=0;field<3;field++){
	if(ToGrid->ElectricField[field] != NULL )
	  delete ToGrid->ElectricField[field];
	ToGrid->ElectricField[field] = new float[MHDeRegionSize[field]];
      }
    }//MHD

  
    if (SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
      delete ToGrid->GravitatingMassFieldParticles;
      ToGrid->GravitatingMassFieldParticles = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassFieldParticles,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
    
    if (SendField == GRAVITATING_MASS_FIELD) {
      delete ToGrid->GravitatingMassField;
      ToGrid->GravitatingMassField = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassField,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }


    if (SendField == POTENTIAL_FIELD) {
      delete ToGrid->PotentialField;
      ToGrid->PotentialField = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->PotentialField,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
    
    if (SendField == ACCELERATION_FIELDS)
      for (dim = 0; dim < GridRank; dim++) {
	delete ToGrid->AccelerationField[dim];
	ToGrid->AccelerationField[dim] = new float[RegionSize];
	FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->AccelerationField[dim],
			     RegionDim, RegionDim+1, RegionDim+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     Zero, Zero+1, Zero+2);
	index += RegionSize;
      }
			  
  }

  if (MyProcessorNumber == ToProcessor)
    delete [] buffer;


  return SUCCESS;
}


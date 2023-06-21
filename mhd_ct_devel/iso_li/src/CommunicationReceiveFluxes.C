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
/  COMMUNICATION ROUTINE: RECEIVE FLUXES FROM ANOTHER PROCESSOR
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
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "error.h"

/* function prototypes */

float ReturnCPUTime();


int CommunicationReceiveFluxes(fluxes *Fluxes, int FromProc, 
			       int NumberOfFields, int Rank)
{

  /* Count space and allocate buffer. */

  int dim1, dim2, field, i, TotalSize = 0, Sizes[MAX_DIMENSION], TempDim;
  int ESizes[3][3];
  for (dim1 = 0; dim1 < Rank; dim1++) {
    int size = 1;
    for (dim2 = 0; dim2 < Rank; dim2++) {
      TempDim = (Fluxes->LeftFluxEndGlobalIndex[dim1][dim2] -
	         Fluxes->LeftFluxStartGlobalIndex[dim1][dim2]) + 1;
      if (dim2 == dim1)
	TempDim = 1;
      size *= TempDim;
    }
    Sizes[dim1] = size;
    TotalSize += 2*size;
  }

  TotalSize *= NumberOfFields;

  /* Space counting for Electric Fields */

  if( MHD_Used == TRUE && MHD_FluxCorrection == TRUE ){
    for(dim1=0;dim1<Rank; dim1++){
      for(field=0;field<3;field++)
	if(field != dim1){
	  ESizes[field][dim1]=1;
	  for( dim2=0; dim2<Rank; dim2++){
	    TempDim = (Fluxes->LeftFluxEndGlobalIndex[dim1][dim2] -
		       Fluxes->LeftFluxStartGlobalIndex[dim1][dim2]) + 1;
	    if( dim2 != dim1 && dim2 != field ) 
	      TempDim++;
	    
	    ESizes[field][dim1]*= TempDim;
	  }
	  TotalSize += 2*ESizes[field][dim1];

	}else{
	  ESizes[field][dim1] = 0;
	}//field
    }//dim1
    
  
  }//MHD 

  float *buffer = new float[TotalSize];

  /* receive into buffer. */
  
#ifdef USE_MPI

  MPI_Status status;
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;

//  MPI_Datatype DataType = MPI_FLOAT;
//  if (sizeof(float) == 8)
//    DataType = MPI_DOUBLE;

  float time1 = ReturnCPUTime();

  ZLAN_START;

  JBPERF_START_MPI_RECV("MPI_Recv",TotalSize,DataType);

  CHECK_MPI_ERROR(MPI_Recv(buffer, TotalSize, DataType, FromProc, 
			   MPI_FLUX_TAG, MPI_COMM_WORLD, &status));

  JBPERF_STOP_MPI_RECV("MPI_Recv",TotalSize,DataType);


  ZLAN_STOP_RECV(12);

  CommunicationTime += ReturnCPUTime() - time1;

#endif /* USE_MPI */

  /* Unpack buffer */

  int index = 0;
  for (dim1 = 0; dim1 < Rank; dim1++)
    for (field = 0; field < NumberOfFields; field++) {
      for (i = 0; i < Sizes[dim1]; i++)
	Fluxes->LeftFluxes[field][dim1][i] = buffer[index++];
      for (i = 0; i < Sizes[dim1]; i++)
	Fluxes->RightFluxes[field][dim1][i] = buffer[index++];
    }

  if( MHD_Used && MHD_FluxCorrection ){
    for(dim1=0;dim1<Rank; dim1++){
      for(field=0;field<3;field++)
	if(field != dim1) {

	  for(i=0;i<ESizes[field][dim1];i++){
	    Fluxes->LeftElectric[field][dim1][i]=buffer[index++];
	  }

	  for(i=0;i<ESizes[field][dim1];i++){
	    Fluxes->RightElectric[field][dim1][i]=buffer[index++];
	  }

	}//field
    }//dim1
  }//mhd
  
  delete buffer;
  
  return SUCCESS;
}

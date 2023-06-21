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
/  GRID CLASS (CHECKS FOR NANS)
/
/  written by: Greg Bryan
/  date:       1999
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
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
#include "fortran.def"
#include "Grid.h"

// Turn the first on to check the gas and dark matter data for nan's.
// Turn the second on just to print the message (i.e. which routine
//     has called it.

#define DEBUG_CHECK_OFF
#define TRACE_OFF

/* function prototypes */

int  ReportMemoryUsage(char *header = NULL);

int grid::DebugCheck(char *message)
{

#ifdef TRACE_ON

  if (ProcessorNumber == MyProcessorNumber)
    fprintf(stderr, "P(%d): %s\n", MyProcessorNumber, message);
    //    ReportMemoryUsage(message);

#endif /* TRACE_ON */

#ifdef DEBUG_CHECK_ON

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Set this to zero but don't the compiler so it doesn't optimize everything
     away. */

  int i,j,k,field,index, Shit=FALSE, FailHard=FALSE;
  
  int ThisIsZero = GridStartIndex[0] - DEFAULT_GHOST_ZONES, size = 1, 
      dim, k1, k2;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (field = 0; field < NumberOfBaryonFields; field++)
    for(k=0;k<GridDimension[2];k++)
      for(j=0;j<GridDimension[1];j++)
	for(i=0;i<GridDimension[0];i++){
	  index=index0(i,j,k);

	  if( (FailHard == TRUE && Shit == FALSE ) || FailHard == FALSE ) {
	    if (BaryonField[field][index+ThisIsZero] != BaryonField[field][index]) {
	      fprintf(stderr, "DebugCheck[%s](Proc %d): BaryonField %d %d %d %d %g\n", message, 
		      MyProcessorNumber, field, i,j,k, BaryonField[field][index]);
	      Shit=TRUE;
	    }
	  }
	}//i,j,k,field

  for (k1 = 0; k1 < GridRank; k1++)
    for (k2 = 0; k2 < NumberOfParticles; k2++)
      if (ParticlePosition[k1][k2] != ParticlePosition[k1][k2+ThisIsZero] ||
	  ParticleVelocity[k1][k2] != ParticleVelocity[k1][k2+ThisIsZero]  ) {
	fprintf(stderr, "DebugCheck[%s](Proc %d): dm %d (%d/%d) %g %g\n", 
		message, MyProcessorNumber, k1, k2, NumberOfParticles,
		ParticlePosition[k1][k2], ParticleVelocity[k1][k2]);
	exit(EXIT_FAILURE);
      }

  if( MHD_Used ){

    
    for(field=0;field<3;field++){
      
      
      for(k=0;k<MagneticDims[field][2];k++)
	for(j=0;j<MagneticDims[field][1];j++)
	  for(i=0;i<MagneticDims[field][0];i++){

	    index=indexba(i,j,k,field);
	  if( (FailHard == TRUE && Shit == FALSE ) || FailHard == FALSE ) {
	    if (MagneticField[field][index+ThisIsZero] != MagneticField[field][index]) {
	      fprintf(stderr, "DebugCheck[%s](Proc %d): MagneticField[%d][%d,%d,%d]= %g\n", message, 
		      MyProcessorNumber, field, i,j,k, MagneticField[field][index]);
	      Shit=TRUE;
	    }
	  }
	  }//i,j,k
      
      if( ElectricField[field] != NULL )
	for(k=0;k<ElectricDims[field][2];k++)
	  for(j=0;j<ElectricDims[field][1];j++)
	    for(i=0;i<ElectricDims[field][0];i++){
	      index=i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);
	  if( (FailHard == TRUE && Shit == FALSE ) || FailHard == FALSE ) {
	      if (ElectricField[field][index+ThisIsZero] != ElectricField[field][index]) {
		fprintf(stderr, "DebugCheck[%s](Proc %d): ElectricField[%d][%d,%d,%d]= %g\n", message, 
			MyProcessorNumber, field, i,j,k, ElectricField[field][index]);
		Shit=TRUE;
	      }
	  }
	    }//i,j,k
      
    }//field

  }//MHD

  if( Shit == TRUE ){
    FILE * Shit = fopen("dataShit","a");
    this->WriteGrid(Shit,"dataShit",666);
    exit(EXIT_FAILURE);
  }

#if 0		
  if (NumberOfBaryonFields > 0)
    for (k1 = 0; k1 < 2+DualEnergyFormalism; k1++)
      for (k2 = 0; k2 < size; k2++)
	if (BaryonField[k1][k2+ThisIsZero] <= 0) {
	  fprintf(stderr, "DebugCheck[%s](Proc %d): <0 %d %d %g\n", message, 
		  MyProcessorNumber, k1, k2, BaryonField[k1][k2]);
	  exit(EXIT_FAILURE);
	}
#endif

#endif /* DEBUG_CHECK_ON */
		
  return SUCCESS;
}

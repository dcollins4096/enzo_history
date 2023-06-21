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
/  GRID CLASS (INTERPOLATE AND OUTPUT GRID TO ARBITRARY TIME)
/
/  written by: Greg Bryan
/  date:       April, 2000
/  modified1:
/
/  PURPOSE:  This routine interpolates grid data to the time passed
/            in and then call the regular grid i/o routine.  It is
/            intended for outputs at arbitrary times, rather than just
/            at the end of a step, as the regular write grid routine
/            assumes.
/
/  RETURNS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hdf4.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "StarParticleData.h"

/* function prototypes */
void my_exit(int status);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/* This macro converts a float and writes it to the local buffer, which,
   when full is written to the file pointed to by fptr. */


int grid::WriteGridInterpolate(FLOAT WriteTime, FILE *fptr, char *base_name, 
			       int grid_id)
{
  
  int dim, i, field, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  /* Compute coefficient factors for linear interpolation in time.
     Note: interp = coef1*Old + coef2*New. */
  
  float coef1 = 0, coef2 = 1;
  if (Time != WriteTime) {
    if (Time <= OldTime) {
      fprintf(stderr, "WGI: fields are at the same time or worse.\n");
      return FAIL;
    } else {
      coef1 = max((Time - WriteTime)/
		  (Time - OldTime), 0.0);
      coef2  = (1.0 - coef1);
    }
  }
  
  /* Interpolate grid to given time and save old variables. */
  
  
  
  float *SavedBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];
  float *SavedMagneticField[3];
  float *SavedCentered[3];


  //float *SavedElectricField[3];
  if (coef2 != 1 && MyProcessorNumber == ProcessorNumber){
    for (field = 0; field < NumberOfBaryonFields; field++) {
      SavedBaryonField[field] = BaryonField[field];
      BaryonField[field] = new float[size];
      if( BaryonField[field] == NULL ){
	fprintf(stderr,"BF = Null!\n");
	my_exit(EXIT_FAILURE);
      }

      //<dbg>

      if( OldBaryonField[field] != NULL){
	for (i = 0; i < size; i++)
	  BaryonField[field][i] = coef1*OldBaryonField[field][i] +
	    coef2*SavedBaryonField[field][i];
      }else{
	for (i = 0; i < size; i++)
	  BaryonField[field][i] =SavedBaryonField[field][i];
      }
      //for (i = 0; i < size; i++)
      //BaryonField[field][i] = coef1*OldBaryonField[field][i] +
      //coef2*SavedBaryonField[field][i];
      //</dbg>
    }
    if(MHD_Used){    
    for(field = 0; field<3;field++){
      
      SavedMagneticField[field] = MagneticField[field];
      SavedCentered[field] = CenteredB[field];
      
      //SavedElectricField[field] = ElectricField[field];
      
      MagneticField[field] = new float[MagneticSize[field]];
      CenteredB[field] = new float[size];
      //ElectricField[field] = new float[ElectricSize[field]];
      
      for(i=0;i<MagneticSize[field]; i++)
	MagneticField[field][i] = coef1*OldMagneticField[field][i] +
	  coef2*SavedMagneticField[field][i];
      for(i=0;i<size;i++)
	CenteredB[field][i] = coef1*OldCenteredB[field][i]+
	  coef2*SavedCentered[field][i];
      /*
	for(i=0;i<ElectricSize[field];i++)
	ElectricField[field][i] = coef1*OldElectricField[field][i] +
	coef2*SavedElectricField[field][i];
      */
    }
  }					
  
}
/* Move particles to given time. */

  float TimeDifference = WriteTime - Time;
  if (this->UpdateParticlePosition(TimeDifference) == FAIL) {
    fprintf(stderr, "Error in grid->UpdateParticlePosition.\n");
    return FAIL;
  }

  /* Write grid (temporarily replace Time with WriteTime). */

  FLOAT SavedTime = Time;
  Time = WriteTime;
  if (this->WriteGrid(fptr, base_name, grid_id) == FAIL) {
    fprintf(stderr, "Error in grid->WriteGrid.\n");
    return FAIL;
  }
  Time = SavedTime;

  /* Move particles back. */

  this->UpdateParticlePosition(-TimeDifference);

  /* Copy saved baryon fields back. */

  if (coef2 != 1 && MyProcessorNumber == ProcessorNumber){
    for (field = 0; field < NumberOfBaryonFields; field++) {
      delete [] BaryonField[field];
      BaryonField[field] = SavedBaryonField[field];
    }
    
    if(MHD_Used){
      for(field=0; field<3; field++){
	delete [] MagneticField[field];
	MagneticField[field] = SavedMagneticField[field];
	delete [] CenteredB[field];
	CenteredB[field] = SavedCentered[field];
	//delete [] ElectricField[field];
	//ElectricField[field] = SavedElectricField[field];
      }
    }
  }
  return SUCCESS;
}


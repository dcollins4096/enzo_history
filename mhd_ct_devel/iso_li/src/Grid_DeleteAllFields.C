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
/  GRID CLASS (REMOVE ALL FIELDS)
/
/  written by: Greg Bryan
/  date:       April, 1996
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

void grid::DeleteAllFields()
{

  int i;

  this->DeleteParticles();

  for (i = 0; i < MAX_DIMENSION; i++) {
    delete [] ParticleAcceleration[i];
    delete [] AccelerationField[i];

    ParticleAcceleration[i]      = NULL;
   AccelerationField[i]         = NULL;
  }
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ParticleAcceleration[MAX_DIMENSION] = NULL;

  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    delete [] BaryonField[i];
    delete [] OldBaryonField[i];
    BaryonField[i]    = NULL;
    OldBaryonField[i] = NULL;
  }


  
  for(i=0; i<3; i++)
    if(OldAccelerationField[i] != NULL ){
      delete OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }

  for(i=0;i<3;i++){
    if(MagneticField[i] != NULL){
      delete [] MagneticField[i];
      MagneticField[i] = NULL;
    }
    if( ElectricField[i] != NULL ){
      delete [] ElectricField[i];
      ElectricField[i] = NULL;
    }
    if( CenteredB[i] != NULL ){
      delete [] CenteredB[i];
      CenteredB[i]     = NULL;
    }
    if(OldMagneticField[i] != NULL){
      delete [] OldMagneticField[i];
      OldMagneticField[i] = NULL;
    }
    if(OldElectricField[i] != NULL){
      delete [] OldElectricField[i];
      OldElectricField[i] = NULL;
    }
    if(OldCenteredB[i] != NULL){
      delete [] OldCenteredB[i];
      OldCenteredB[i]     = NULL;
    }
    if( Current[i] != NULL ){
      delete[] Current[i];
      Current[i] = NULL;
    }
    if( AvgElectricField[i] != NULL ){
      delete[] AvgElectricField[i];
      AvgElectricField[i] = NULL;
    }
  }

  
  delete [] PotentialField;
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;

  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;

}
